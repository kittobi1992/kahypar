/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#pragma once

#include <algorithm>
#include <limits>
#include <queue>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/utils/randomize.h"

namespace kahypar {
class FlowCutter : public IInitialPartitioner,
                   private InitialPartitionerBase<FlowCutter >{
  using Base = InitialPartitionerBase<FlowCutter >;
  friend Base;

 public:
  FlowCutter(Hypergraph& hypergraph, Context& context) :
      Base(hypergraph, context) {
    ALWAYS_ASSERT(context.partition.k == 2,context.partition.k);
  }

  ~FlowCutter() override = default;

  FlowCutter(const FlowCutter&) = delete;
  FlowCutter& operator= (const FlowCutter&) = delete;

  FlowCutter(FlowCutter&&) = delete;
  FlowCutter& operator= (FlowCutter&&) = delete;

 private:
  void partitionImpl() override final {
    Base::multipleRunsInitialPartitioning();
  }

  void initialPartition() {
    // your implementation here
    PartitionID part = 0;
    for (const HypernodeID hn : _hg.nodes()){
      _hg.setNodePart(hn, part);
      part = (part + 1) % _context.partition.k;
    }
    Base::performFMRefinement();
  }
  using Base::_hg;
  using Base::_context;
  using Base::kInvalidNode;
};
}  // namespace kahypar
