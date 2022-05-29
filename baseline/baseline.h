//
// Created by lsy on 2022/3/24.
//

#ifndef TEMPORAL_ABCORE_BASELINE_H
#define TEMPORAL_ABCORE_BASELINE_H

#include "../bigraph/bigraph.h"

auto index_baseline(BiGraph& g) -> void ;
auto baseline_query(const int& alpha, const int& beta, const int& ts, const int& te, BiGraph& g, vector<bool>& node_u,vector<bool>& node_v) -> void ;

auto back_del_edges(BiGraph& g, BiGraph& tg, const int& ts, const int& tmax) -> void;
auto compute_core_neighbor(const vid_t& ts, vector<unordered_map<int, int>>& cn, const num_t& num, vector<vector<pair<vid_t,vid_t>>>& nu) -> void;

#endif //TEMPORAL_ABCORE_BASELINE_H
