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

#include <limits>
#include <string>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/macros.h"
#include "kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "kahypar/partition/coarsening/policies/rating_community_policy.h"
#include "kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "kahypar/partition/coarsening/policies/rating_tie_breaking_policy.h"
#include "kahypar/partition/coarsening/vertex_pair_rater.h"

namespace kahypar {
template <class ScorePolicy = HeavyEdgeScore,
          class HeavyNodePenaltyPolicy = NoWeightPenalty,
          class CommunityPolicy = UseCommunityStructure,
          class AcceptancePolicy = BestRatingPreferringUnmatched<>,
          typename RatingType = RatingType>
class AlgebraicDistanceCoarsener final : public ICoarsener,
                                         private VertexPairCoarsenerBase<>{
 private:
  static constexpr bool debug = false;

  using Base = VertexPairCoarsenerBase;
  using Rater = VertexPairRater<ScorePolicy,
                                HeavyNodePenaltyPolicy,
                                CommunityPolicy,
                                AcceptancePolicy,
                                RatingType>;
  using Rating = typename Rater::Rating;

 public:
  AlgebraicDistanceCoarsener(Hypergraph& hypergraph, const Context& context,
                             const HypernodeWeight weight_of_heaviest_node) :
    Base(hypergraph, context, weight_of_heaviest_node),
    _rater(_hg, _context) { }

  ~AlgebraicDistanceCoarsener() override = default;

  AlgebraicDistanceCoarsener(const AlgebraicDistanceCoarsener&) = delete;
  AlgebraicDistanceCoarsener& operator= (const AlgebraicDistanceCoarsener&) = delete;

  AlgebraicDistanceCoarsener(AlgebraicDistanceCoarsener&&) = delete;
  AlgebraicDistanceCoarsener& operator= (AlgebraicDistanceCoarsener&&) = delete;

 private:
  void coarsenImpl(const HypernodeID limit) override final {
    // TODO: implement
    LOG << "Hello World";
    // To contract vertex pair (u,v) call:
    // performContraction(u, v);
  }

  bool uncoarsenImpl(IRefiner& refiner) override final {
    return doUncoarsen(refiner);
  }

  using Base::_pq;
  using Base::_hg;
  using Base::_context;
  using Base::_history;
  Rater _rater;
};
}  // namespace kahypar
