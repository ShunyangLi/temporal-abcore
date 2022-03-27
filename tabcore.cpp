//
// Created by lsy on 2022/3/24.
// used to create temporal abcore
//

#include "tabcore.h"
#include "bigraph.h"
#include "abcore.h"

// computer core neighbor in ts to te
// the core neighbor is the number of neighbors that in ts and te
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
 * delete edges for peeling ts to te
 * @param g
 */
auto compute_del_edges(BiGraph& g, const vid_t& ts, const vid_t& te, vector<bool>& vu, vector<bool>& vv) -> void {
    auto tg = g;
    auto q = queue<pair<vid_t, vid_t>>();

    // to do the loop delete edges
    for (auto _te = tg.tmax; _te >= ts; --_te) {

        // because there are may one more than one edge that between ts and te
        for (auto index = tg.edges_idx[_te]; index < tg.edges_idx[_te + 1]; ++index) {
            auto u = tg.edges[index].first;
            auto v = tg.edges[index].second;

            // then we process u and v
            -- tg.ucn[u][v];
            -- tg.vcn[v][u];

            // it is one edge, so just detect only once
            if (tg.ucn[u][v] < tg.left_index.size() || tg.vcn[v][u] < tg.right_index.size()) {
                q.push(std::make_pair(u, v));
                vu[u] = true;
                vv[v] = true;
            }
        }

        // then try to delete the edge, and check whether the core number changed
        while (!q.empty()) {
            auto const u_index = g.left_index[u];
            auto const v_index = g.right_index[v];
            auto v_size = g.neighbor_v2[v].size();
            auto u_size = g.neighbor_v1[u].size();

            auto afftect_u = vector<vector<vid_t>>(g.left_index.size());
            auto afftect_v = vector<vector<vid_t>>(g.right_index.size());

            update_bicore_index(g,u,v,addition, afftect_u, afftect_v);


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
    auto vu = vector<bool>(g.num_v1, false);
    auto vv = vector<bool>(g.num_v2, false);

    for (auto ts = 0; ts < g.tmax; ++ ts) {
        // then working here

        if (ts == g.tmax - 1) break;
        compute_core_neighbor(0, g.ucn, g.num_v1, g.tnu);
        compute_core_neighbor(0, g.vcn, g.num_v2, g.tnv);



    }

}