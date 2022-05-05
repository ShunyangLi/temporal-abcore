//
// Created by lsy on 2022/5/2.
//

#ifndef TEMPORAL_ABCORE_TABCORE_BASELINE_H
#define TEMPORAL_ABCORE_TABCORE_BASELINE_H

#include "../bigraph/bigraph.h"
#include "../config/config.h"

auto tabcore_baseline(BiGraph& g) -> void;
auto query(const int& alpha, const int& beta, const int& ts, const int& te, BiGraph& g, vector<bool>& node_u,
           vector<bool>& node_v) -> void;

#endif //TEMPORAL_ABCORE_TABCORE_BASELINE_H
