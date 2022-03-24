//
// Created by lsy on 2022/3/24.
// used to create temporal abcore
//

#include "tabcore.h"
#include "bigraph.h"
#include "abcore.h"

// computer core neighbor in ts to te
auto compute_core_neighbor(const vid_t& ts, vector<unordered_map<int, int>>& cn,
                           const num_t& num, vector<vector<pair<vid_t,vid_t>>> nu) -> void {

    for (auto u = 0; u < num; u ++) {
        cn[u].clear();
        for (auto index = nu[u].size() - 1; ts >= 0; index --) {
            auto v = nu[u][index].first;
            auto t = nu[u][index].second;

            // because time from largest to smallest
            if (t < ts) break;

            cn[u][v] += 1;

        }
    }
}


/**
 * baseline for temporal abcore
 * @param g graph object
 */
auto index_baseline(BiGraph& g) -> void  {

    cout << "starting abcore decomposition" << endl;
    // ab core decomposition
    coreIndexKCore(g);

    // init the index with vertex
    g.u_index.resize(g.num_v1);
    g.v_index.resize(g.num_v2);

    // initial the upper index
    for (auto u = 0; u < g.num_v1; u ++) {
        g.u_index[u].resize(g.left_index[u].size());
        for (auto alpha = 1; alpha < g.left_index[u].size(); alpha ++) {
            g.u_index[u][alpha].resize(g.left_index[u][alpha]);
        }
    }

    // initial the lower index
    for (auto v = 0; v < g.num_v2; v ++) {
        g.v_index.resize(g.right_index[v].size());
        for (auto beta = 1; beta < g.right_index[v].size(); beta ++) {
            g.v_index[v][beta].resize(g.right_index[v][beta]);
        }
    }

    // init the core neighbor number
    compute_core_neighbor(0, g.ucn, g.num_v1, g.tnu);
    compute_core_neighbor(0, g.vcn, g.num_v2, g.tnv);

    for (auto ts = 0; ts < g.tmax; ++ ts) {
        // then working here

        if (ts == g.tmax - 1) break;
        compute_core_neighbor(0, g.ucn, g.num_v1, g.tnu);
        compute_core_neighbor(0, g.vcn, g.num_v2, g.tnv);

    }

}