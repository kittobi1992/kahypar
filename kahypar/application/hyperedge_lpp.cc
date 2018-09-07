/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>

#include "kahypar/application/command_line_options.h"
#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/io/partitioning_output.h"
#include "kahypar/io/sql_plottools_serializer.h"
#include "kahypar/macros.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/utils/math.h"
#include "kahypar/utils/randomize.h"

using kahypar::HighResClockTimepoint;
using kahypar::Context;
using kahypar::HypernodeID;
using kahypar::PartitionID;
using kahypar::HyperedgeWeight;
using kahypar::HypernodeWeight;

int main(int argc, char* argv[]) {
  Context context;

  kahypar::processCommandLineInput(context, argc, argv);
  kahypar::Randomize::instance().setSeed(context.partition.seed);

  kahypar::Hypergraph hypergraph(
      kahypar::io::createHypergraphFromFile(context.partition.graph_filename,
                                            context.partition.k));
  context.setupPartWeights(hypergraph.totalWeight());
  LOG << V(hypergraph.totalWeight());
  std::chrono::duration<double> elapsed_time(0);

  const HighResClockTimepoint complete_start = std::chrono::high_resolution_clock::now();


  // Label Propagation Partitioning (LPP) algorithm of Jiang et al.
  // W. Jiang, J. Qi, J. X. Yu, J. Huang and R. Zhang,
  // "HyperX: A Scalable Hypergraph Framework,"
  // in IEEE Transactions on Knowledge and Data Engineering.
  // doi: 10.1109/TKDE.2018.2848257
  // URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8391739&isnumber=4358933

  // + 1 for initial sentinel
  // used to compute L(h)
  std::vector<HypernodeID> part_arities(context.partition.k + 1, 0.0);

  std::vector<HypernodeID> pins_in_part(context.partition.k, 0);

  // used to compute L(v)
  std::vector<double> he_scores(context.partition.k, 0.0);

  // sentinel initialization
  for (const auto he : hypergraph.edges()) {
    hypergraph.hyperedgeData(he).label = context.partition.k;
    part_arities[context.partition.k] += hypergraph.edgeSize(he);
  }

  // initialization
  for (const auto hn : hypergraph.nodes()) {
    hypergraph.hypernodeData(hn).label = kahypar::Randomize::instance().getRandomInt(0, context.partition.k - 1);
  }

  // $\overline{A}$ for L(v) computation
  const double avg_part_arity = hypergraph.initialNumPins() /static_cast<double>(context.partition.k);
  const double avg_part_arity_squared = std::pow(static_cast<double>(avg_part_arity), 2);

  std::vector<PartitionID> max_pins_parts;
  std::vector<PartitionID> max_scores;

  int iteration  = 0;
  HyperedgeWeight current_cut = std::numeric_limits<HyperedgeWeight>::max();
  HyperedgeWeight best_cut = std::numeric_limits<HyperedgeWeight>::max();

  // TODO: Try iterating in increasing order of degree/net size

  do {
    best_cut = current_cut;
    ASSERT(iteration == 0 || part_arities[context.partition.k] == 0);
    /////////////////////////////////////////////////////////////////
    // STEP 1: Update hyperedge labels
    /////////////////////////////////////////////////////////////////
    // Goal:
    // Minimize the number of replicas:
    // Nets adapt the label which is most common among their pins.
    // Ties are broken randomly.
    for (const auto he : hypergraph.edges()) {
      for (const auto pin : hypergraph.pins(he)) {
        ++pins_in_part[hypergraph.hypernodeData(pin).label];
      }

      HypernodeID max_pins = std::numeric_limits<HypernodeID>::min();
      max_pins_parts.clear();
      for (PartitionID part = 0; part != context.partition.k; ++part) {
        if (pins_in_part[part] > max_pins) {
          max_pins = pins_in_part[part];
          max_pins_parts.clear();
          max_pins_parts.push_back(part);
        } else if (pins_in_part[part] == max_pins) {
          max_pins_parts.push_back(part);
        }
      }
      ASSERT(!max_pins_parts.empty(),"");
      ASSERT(std::all_of(max_pins_parts.cbegin(), max_pins_parts.cend(),
                         [&](PartitionID i) {return pins_in_part[i] == max_pins;}));
      std::fill(pins_in_part.begin(), pins_in_part.end(), 0);

      const PartitionID max_part = max_pins_parts[kahypar::Randomize::instance().getRandomInt(0, max_pins_parts.size() - 1)];

      part_arities[hypergraph.hyperedgeData(he).label] -= hypergraph.edgeSize(he);
      hypergraph.hyperedgeData(he).label = max_part;
      part_arities[max_part] += hypergraph.edgeSize(he);
    }
    ASSERT((hypergraph.initialNumPins()) ==
           std::accumulate(part_arities.begin(), part_arities.end(), 0.0));

    /////////////////////////////////////////////////////////////////
    // Step 2: Update vertex labels
    /////////////////////////////////////////////////////////////////
    // Goal:
    // - balance net sizes across blocks -> assign vertex to blocks with smaller sum of pins
    // - minimize the number of replicas -> assign vertex to block containing most incident nets
    for (const auto hn : hypergraph.nodes()) {
      ASSERT(std::all_of(he_scores.begin(), he_scores.end(),
                         [](const double score) { return score == 0.0;}));
      for (const auto incident_he : hypergraph.incidentEdges(hn)) {
        ++he_scores[hypergraph.hyperedgeData(incident_he).label];
      }

      for (PartitionID part = 0; part != context.partition.k; ++part) {
        he_scores[part] =  he_scores[part] *
                           std::exp((avg_part_arity_squared - std::pow(part_arities[part], 2))
                                    / avg_part_arity_squared);
      }

      double  max_score = std::numeric_limits<double>::lowest();
      max_scores.clear();
      for (PartitionID part = 0; part != context.partition.k; ++part) {
        if (he_scores[part] > max_score) {
          max_score = he_scores[part];
          max_scores.clear();
          max_scores.push_back(part);
        } else if (he_scores[part] == max_score) {
          max_scores.push_back(part);
        }
      }
      ASSERT(!he_scores.empty(),"");
      ASSERT(std::all_of(max_scores.cbegin(), max_scores.cend(),
                         [&](PartitionID i) {return he_scores[i] == max_score;}));
      std::fill(he_scores.begin(), he_scores.end(), 0.0);

      const PartitionID max_part = max_scores[kahypar::Randomize::instance().getRandomInt(0, max_scores.size() - 1)];

      hypergraph.hypernodeData(hn).label = max_part;
    }

    if (iteration == 0) {
      for (const auto hn : hypergraph.nodes()) {
        ASSERT(hypergraph.hypernodeData(hn).label < context.partition.k,"");
        hypergraph.setNodePart(hn, hypergraph.hypernodeData(hn).label);
      }
      hypergraph.initializeNumCutHyperedges();
    } else{
      for (const auto hn : hypergraph.nodes()) {
        ASSERT(hypergraph.hypernodeData(hn).label < context.partition.k,"");
        if (hypergraph.partID(hn) != hypergraph.hypernodeData(hn).label) {
          hypergraph.changeNodePart(hn,hypergraph.partID(hn),hypergraph.hypernodeData(hn).label);
        }
      }
    }
    LOG << V(iteration);
    kahypar::io::printObjectives(hypergraph, context);
    current_cut = kahypar::metrics::hyperedgeCut(hypergraph);
    ++iteration;
  } while (current_cut != best_cut);

  const HighResClockTimepoint complete_end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_seconds = complete_end - complete_start;

  LOG << "********************************************************************************";
  LOG << "*                          FINAL Partitioning Result                           *";
  LOG << "********************************************************************************";
  for (const auto hn : hypergraph.nodes()) {
    ASSERT(std::any_of(hypergraph.incidentEdges(hn).first,
                       hypergraph.incidentEdges(hn).second,
                       [&](const kahypar::HyperedgeID he) {
                         return hypergraph.hyperedgeData(he).label == hypergraph.partID(hn);
                       }), "");
  }

  LOG << "#########################################################################";
  LOG << "Evaluation based on labels";
  LOG << "#########################################################################";
  std::vector<HyperedgeWeight> he_label_weights(context.partition.k, 0);
  for (const auto he : hypergraph.edges()) {
    ++he_label_weights[hypergraph.hyperedgeData(he).label];
  }
  for (PartitionID part = 0; part != context.partition.k; ++part){
    LOG << "Hyperedgeload of" << part << "=" <<he_label_weights[part];
  }

  std::vector<HypernodeWeight> hn_label_weights(context.partition.k, 0);
  for (const auto hn : hypergraph.nodes()) {
    ++hn_label_weights[hypergraph.hypernodeData(hn).label];
  }
  for (PartitionID part = 0; part != context.partition.k; ++part){
    LOG << "Vertexload of" << part << "=" <<hn_label_weights[part];
  }

  LOG << "#########################################################################";
  LOG << "Evaluation based on vertex assignment";
  LOG << "#########################################################################";
  std::vector<HyperedgeWeight> he_part_weights(context.partition.k, 0);
  for (const auto he : hypergraph.edges()) {
    for (const auto part : hypergraph.connectivitySet(he)) {
      he_part_weights[part] += hypergraph.edgeWeight(he);
    }
  }

  for (PartitionID part = 0; part != context.partition.k; ++part){
    LOG << "Hyperedgeload of" << part << "=" <<he_part_weights[part];
  }

  std::vector<HypernodeWeight> hn_part_weights(context.partition.k, 0);
  for (const auto hn : hypergraph.nodes()) {
    ++hn_part_weights[hypergraph.partID(hn)];
  }
  for (PartitionID part = 0; part != context.partition.k; ++part){
    LOG << "Vertexload of" << part << "=" <<hn_part_weights[part];
  }

  LOG << "#########################################################################";
  LOG << "Evaluation based on hypergaph";
  LOG << "#########################################################################";
  kahypar::io::printPartSizesAndWeights(hypergraph);



  kahypar::io::writePartitionFile(hypergraph,
                                  context.partition.graph_partition_filename);

  // In case a time limit is used, the last partitioning step is already serialized
  if (context.partition.sp_process_output && context.partition.time_limit == 0) {
    kahypar::io::serializer::serialize(context, hypergraph, elapsed_seconds, iteration);
  }

  return 0;
}
