//
// Created by lsy on 5/5/22.
//

#ifndef TEMPORAL_ABCORE_ONLINE_H
#define TEMPORAL_ABCORE_ONLINE_H

#include "../bigraph/bigraph.h"

auto online_peeling (const int& alpha, const int& beta, const int& ts, const int& te,
                     BiGraph& g, vector<bool>& node_u, vector<bool>& node_v) -> void;

#endif //TEMPORAL_ABCORE_ONLINE_H
