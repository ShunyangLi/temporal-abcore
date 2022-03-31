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
                           const num_t& num, vector<vector<pair<vid_t,vid_t>>>& nu) -> void {

    for (auto u = 0; u < num; u ++) {
        cn[u].clear();
        for (int index = nu[u].size() - 1; index >= 0; index --) {
            auto v = nu[u][index].first;
            auto t = nu[u][index].second;
            // because time from largest to smallest
            if (t < ts) break;

            cn[u][v] += 1;

        }
    }
}

// just set the index value based on ts and te, alpha, and beta
auto update_index(vector<vector<vector<vector<pair<vid_t,vid_t>>>>>& index,
                  const int& ts, const int& te,
                  const vid_t& u, const int& alpha, const int& beta, BiGraph& g) -> void {

    if (index[u][alpha][beta].empty()) {
        index[u][alpha][beta].push_back(make_pair(g.time_new_to_old[ts], g.time_new_to_old[te]));
    } else {
        auto times = index[u][alpha][beta].back();
        if (te > times.second) {
            index[u][alpha][beta].push_back(make_pair(g.time_new_to_old[ts], g.time_new_to_old[te]));
        }
    }
}

/**
 * delete edges for peeling ts to te
 * @param g
 */
auto compute_del_edges(BiGraph& g, BiGraph& tg,
                       const int& ts, vector<bool>& vu, vector<bool>& vv, ) -> void {
//    auto tg = g;
    auto q = queue<pair<vid_t, vid_t>>();

    // to do the loop delete edges
    for (auto _te = tg.tmax - 1; _te >= ts; _te --) {

        // because there are may one more than one edge that between ts and te
        for (auto index = tg.edges_idx[_te]; index < tg.edges_idx[_te + 1]; ++index) {
            auto u = tg.edges[index].first;
            auto v = tg.edges[index].second;

            // then we process u and v
            -- tg.ucn[u][v];
            -- tg.vcn[v][u];

            if (tg.ucn[u][v] == 0) tg.ucn[u].erase(v);
            if (tg.vcn[v][u] == 0) tg.vcn[v].erase(u);

            // when the neighbors of u is less then a, then update
            // it is one edge, so just detect only once
            if (tg.ucn[u].size() < tg.left_index.size() || tg.vcn[v].size() < tg.right_index.size()) {
                if (vu[u] && vv[v]) continue;
                q.push(std::make_pair(u, v));
                vu[u] = true;
                vv[v] = true;
            }
        }

        // then try to delete the edge, and check whether the core number changed
        while (!q.empty()) {
            auto edge = q.front();
            q.pop();
            auto tu = edge.first;
            auto tv = edge.second;

            auto const u_alpha_offset = tg.left_index;
            auto const v_beta_offset = tg.right_index;

            auto au = unordered_map<vid_t, vector<vid_t>>();
            auto av = unordered_map<vid_t, vector<vid_t>>();

            update_bicore_index(tg, tu, tv, 0, au, av);

            // then we just check whether the the core number is changed
            for (const pair<vid_t, vector<vid_t>>& it : au) {
                auto u = it.first;
                auto changed_alpha = it.second;

                // the alpha value of u becomes smaller
                if (u_alpha_offset[u].size() > tg.left_index[u].size()) {
                    // then record it.
                    for (auto alpha =  u_alpha_offset[u].size() - 1; alpha > tg.left_index[u].size() - 1; --alpha) {
                        auto beta = u_alpha_offset[u][alpha];
                        update_index(g.u_index, ts, _te, u, alpha, beta, g);
                    }
                }

                for (auto alpha  = tg.left_index[u].size() - 1; alpha >= 1; --alpha) {
                    if (tg.left_index[u][alpha] != u_alpha_offset[u][alpha]) {
                        auto beta = u_alpha_offset[u][alpha];
                        update_index(g.u_index, ts, _te, u, alpha, beta, g);
                    }
                }
            }

            for (const pair<vid_t, vector<vid_t>>& it : av) {
                auto v = it.first;
                auto changed_alpha = it.second;

                // the alpha value of u becomes smaller
                if (v_beta_offset[v].size() > tg.right_index[v].size()) {
                    // then record it.
                    for (auto beta =  v_beta_offset[v].size() - 1; beta > tg.right_index[v].size() - 1; --beta) {
                        auto alpha = v_beta_offset[v][beta];
                        update_index(g.v_index, ts, _te, v, beta, alpha, g);
                    }
                }

                for (int beta  = tg.right_index[v].size() - 1; beta >= 1; beta --) {
                    if (tg.right_index[v][beta] != v_beta_offset[v][beta]) {
                        auto alpha = v_beta_offset[v][beta];
                        update_index(g.v_index, ts, _te, v, beta, alpha, g);
                    }
                }
            }
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
        g.u_index[u].resize(g.left_index[u].size() + 1);
        for (auto alpha = 1; alpha < g.left_index[u].size(); alpha ++) {
            g.u_index[u][alpha].resize(g.left_index[u][alpha] + 1);
        }
    }

    // initial the lower index
    for (auto v = 0; v < g.num_v2; v ++) {
        g.v_index[v].resize(g.right_index[v].size() + 1);
        for (auto beta = 1; beta < g.right_index[v].size(); beta ++) {
            g.v_index[v][beta].resize(g.right_index[v][beta] + 1);
        }
    }

    // init the core neighbor number
    // compute_core_neighbor(0, g.ucn, g.num_v1, g.tnu);
    // compute_core_neighbor(0, g.vcn, g.num_v2, g.tnv);
    auto vu = vector<bool>(g.num_v1, false);
    auto vv = vector<bool>(g.num_v2, false);

    for (auto ts = 0; ts < g.tmax; ++ ts) {
        // then working here

        if (ts == g.tmax - 1) break;
        compute_core_neighbor(0, g.ucn, g.num_v1, g.tnu);
        compute_core_neighbor(0, g.vcn, g.num_v2, g.tnv);

        auto tg = g;
        compute_del_edges(g, tg, ts, vu, vv);

        // delete the visited edges
        if (ts < g.tmax) compute_del_edges(g, g, ts, vu, vv);

    }

    cout << "finished baseline" << endl;
}