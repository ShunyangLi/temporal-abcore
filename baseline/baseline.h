//
// Created by lsy on 2022/3/24.
//

#ifndef TEMPORAL_ABCORE_BASELINE_H
#define TEMPORAL_ABCORE_BASELINE_H

#include "../bigraph/bigraph.h"

auto index_baseline(BiGraph& g) -> void ;
auto baseline_query(const int& alpha, const int& beta, const int& ts, const int& te, BiGraph& g, vector<bool>& node_u,
                    vector<bool>& node_v) -> void ;


#endif //TEMPORAL_ABCORE_BASELINE_H
