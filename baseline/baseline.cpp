//
// Created by lsy on 2022/3/24.
// used to create temporal abcore
//

#include "baseline.h"
#include "../bigraph/bigraph.h"
#include "../abcore/abcore.h"
#include "../config/config.h"

/**
 * compute the vertex storage index size
 * @param g graph objext
 */
auto vertex_index_size(BiGraph& g) -> void {
    double idx_size = 0;

    // add the vertices size
    idx_size += sizeof(int) * g.num_v1;

    // then for each vertex compute the real usage
    for (auto u = 0; u < g.num_v1; u ++) {
        for (auto alpha = 1; alpha < g.u_index[u].size(); alpha ++) {
            for (auto beta = 1; beta < g.u_index[u][alpha].size(); beta ++) {
                idx_size += g.u_index[u][alpha][beta].size() * 2 * sizeof(int);
            }
        }
    }

    idx_size += sizeof(int) * g.num_v2;
    for (auto v = 0; v < g.num_v2; v ++) {
        for (auto beta = 1; beta < g.v_index[v].size(); beta ++ ) {
            for (auto alpha = 1; alpha < g.v_index[v][beta].size(); alpha ++) {
                idx_size += g.v_index[v][beta][alpha].size() * 2 * sizeof(int );
            }
        }
    }

    double graph_size = g.num_edges * 3 * sizeof(int) ;

#ifdef MBS
    cout << "Graph size: " << graph_size / 1024 / 1024 << " MB." << endl;
    cout << "Index size: " << double (idx_size / 1024 / 1024) << " MB." << endl;
#else
    cout << "Graph size: " << graph_size << " bytes." << endl;
    cout << "Index size: " << idx_size  << " bytes." << endl;
#endif
}


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
auto update_index(vector<vector<vector<vector<pair<int,int>>>>>& index,
                  const int& ts, const int& te,
                  const vid_t& u, const int& alpha, const int& beta, BiGraph& g) -> void {

    if (index[u][alpha][beta].empty()) {
//        index[u][alpha][beta].push_back(make_pair(g.time_new_to_old[ts], g.time_new_to_old[te]));
        index[u][alpha][beta].push_back(make_pair(ts, te));
    } else {
        auto times = index[u][alpha][beta].back();
        if (te > times.second) {
            index[u][alpha][beta].push_back(make_pair(ts, te));
//            index[u][alpha][beta].push_back(make_pair(g.time_new_to_old[ts], g.time_new_to_old[te]));
        }
    }
}

/**
 * delete edges for peeling ts to te
 * @param g
 */
auto back_del_edges(BiGraph& g, BiGraph& tg,
                    const int& ts, const int& tmax) -> void {
//    auto tg = g;
    auto q = queue<pair<vid_t, vid_t>>();

    // to do the loop delete edges
    for (auto _te = tmax; _te >= ts; _te --) {

        // because there are may one more than one edge that between ts and te
        for (auto index = tg.edges_idx[_te]; index < tg.edges_idx[_te + 1]; ++index) {
            auto u = tg.edges[index - 1].first;
            auto v = tg.edges[index - 1].second;

            // then we process u and v
            -- tg.ucn[u][v];
            -- tg.vcn[v][u];

            if (tg.ucn[u][v] == 0) tg.ucn[u].erase(v);
            if (tg.vcn[v][u] == 0) tg.vcn[v].erase(u);

            // when the neighbors of u is less then a, then update
            // it is one edge, so just detect only once
            if (tg.ucn[u].size() < tg.left_index[u].size() - 1 || tg.vcn[v].size() < tg.right_index[v].size() - 1) {
                q.push(std::make_pair(u, v));
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

            auto au = vector<vid_t>();
            auto av = vector<vid_t>();


            update_bicore_index(tg, tu, tv, DELETION, au, av);

            // then we just check whether the the core number is changed
            for (auto const & u : set<vid_t>(au.begin(), au.end())) {
                for (auto alpha = u_alpha_offset[u].size() - 1; alpha >= 1; --alpha) {
                    auto beta = u_alpha_offset[u][alpha];

                    for (; beta >= 1; beta --) {
                        if (alpha > tg.left_index[u].size() - 1) update_index(g.u_index, ts, _te, u, alpha, beta, g);
                        else {
                            if (tg.left_index[u][alpha] != u_alpha_offset[u][alpha])
                                update_index(g.u_index, ts, _te, u, alpha, beta, g);
                            else break;
                        }
                    }
                }
            }

            for (auto const & v: set<vid_t>(av.begin(), av.end())) {
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
 * delete the used edge from front to end
 */
auto advance_del_edge(BiGraph& g, const int& ts, const int& _te) -> void {
//    auto tg = g;
    auto q = queue<pair<vid_t, vid_t>>();

    auto index = g.edges_idx[_te];
    if (_te - 1 < 0) index = 0;
    else index = g.edges_idx[_te - 1];

    // because there are may one more than one edge that between ts and te
    for (; index < g.edges_idx[_te + 1]; ++index) {
        auto u = g.edges[index].first;
        auto v = g.edges[index].second;

        for (auto it = g.tnu[u].begin(); it != g.tnu[u].end(); ++it) {
            if (it->first == v && it->second == _te) {
                g.tnu[u].erase(it);
                break;
            }
        }

        for (auto it = g.tnv[v].begin(); it != g.tnv[v].end(); ++it) {
            if (it->first == u && it->second == _te) {
                g.tnv[v].erase(it);
                break;
            }
        }


        // when the neighbors of u is less then a, then update
        // it is one edge, so just detect only once
        if (g.tnu[u].size() < g.left_index[u].size() - 1 || g.tnv[v].size() < g.right_index[v].size() - 1) {
            q.push(std::make_pair(u, v));
        }
    }

    // then try to delete the edge, and check whether the core number changed
    while (!q.empty()) {
        auto edge = q.front();
        q.pop();
        auto tu = edge.first;
        auto tv = edge.second;

        auto const u_alpha_offset = g.left_index;
        auto const v_beta_offset = g.right_index;

        auto au = vector<vid_t>();
        auto av = vector<vid_t>();

        update_bicore_index(g, tu, tv, DELETION, au, av);

        // then we just check whether the the core number is changed
        for (auto const& u : au) {
            // the alpha value of u becomes smaller
            if (u_alpha_offset[u].size() > g.left_index[u].size()) {
                // then record it.
                for (auto alpha =  u_alpha_offset[u].size() - 1; alpha > g.left_index[u].size() - 1; --alpha) {
                    auto beta = u_alpha_offset[u][alpha];
                    if (!g.u_index[u][alpha][beta].empty() && g.u_index[u][alpha][beta].back().second == END) continue;
                    g.u_index[u][alpha][beta].push_back(make_pair(g.time_new_to_old[ts + 1], END));
                }
            }

            for (auto alpha  = g.left_index[u].size() - 1; alpha >= 1; --alpha) {
                if (g.left_index[u][alpha] != u_alpha_offset[u][alpha]) {
                    auto beta = u_alpha_offset[u][alpha];
                    if (!g.u_index[u][alpha][beta].empty() && g.u_index[u][alpha][beta].back().second == END) continue;
                    g.u_index[u][alpha][beta].push_back(make_pair(g.time_new_to_old[ts + 1], END));
                }
            }
        }

        for (auto const& v : av) {
            // the alpha value of u becomes smaller
            if (v_beta_offset[v].size() > g.right_index[v].size()) {
                // then record it.
                for (auto beta =  v_beta_offset[v].size() - 1; beta > g.right_index[v].size() - 1; --beta) {
                    auto alpha = v_beta_offset[v][beta];
                    if (!g.v_index[v][beta][alpha].empty() && g.v_index[v][beta][alpha].back().second == END) continue;
                    g.v_index[v][beta][alpha].push_back(make_pair(g.time_new_to_old[ts + 1], END));

                }
            }

            for (int beta  = g.right_index[v].size() - 1; beta >= 1; beta --) {
                if (g.right_index[v][beta] != v_beta_offset[v][beta]) {
                    auto alpha = v_beta_offset[v][beta];
                    if (!g.v_index[v][beta][alpha].empty() && g.v_index[v][beta][alpha].back().second == END) continue;
                    g.v_index[v][beta][alpha].push_back(make_pair(g.time_new_to_old[ts + 1], END));
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
#ifdef TIME
    auto start = chrono::system_clock::now();
#endif

    for (auto ts = 0; ts < g.tmax; ++ ts) {
        // then working here
        if (ts == g.tmax - 1) break;

        // count the neighbor in the time interval ts to te
        compute_core_neighbor(ts, g.ucn, g.num_v1, g.tnu);
        compute_core_neighbor(ts, g.vcn, g.num_v2, g.tnv);

        auto tg = g;
        back_del_edges(g, tg, ts, g.tmax - 1);

        // delete the visited edges
        if (ts < g.tmax) advance_del_edge(g, ts, ts);

    }

#ifdef TIME
    cout << "finished baseline" << endl;
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "construction: " << elapsed_seconds.count() << endl;
#endif

#ifdef INDEX_SIZE
    vertex_index_size(g);
#endif

}

auto baseline_query(const int& alpha, const int& beta, const int& ts, const int& te, BiGraph& g, vector<bool>& node_u,
                    vector<bool>& node_v) -> void {

#ifdef TIME
    auto start = chrono::system_clock::now();
#endif

    node_u = vector<bool>(g.num_v1, false);
    node_v = vector<bool>(g.num_v2, false);

    for (auto u = 0; u < g.u_index.size(); u ++) {
        if (g.u_index[u].size() < alpha) continue;
        if (g.u_index[u][alpha].size() < beta) continue;
        if (g.u_index[u][alpha][beta].empty()) continue;

        for (auto const& ti : g.u_index[u][alpha][beta]) {
            if (ti.first <= ts) {
                if (ti.second == END) {
                    node_u[u] = false;
                    break;
                }
                if (ti.second <= te) node_u[u] = true;
            } else {
                break;
            }
        }
    }

    for (auto v = 0; v < g.v_index.size(); v++) {
        if (g.v_index[v].size() < beta) continue;
        if (g.v_index[v][beta].size() < alpha) continue;
        if (g.v_index[v][beta][alpha].empty()) continue;

        for (auto const& ti : g.v_index[v][beta][alpha]) {
            if (ti.first <= ts) {
                if (ti.second == END) {
                    node_v[v] = false;
                    break;
                }
                if (ti.second <= te) node_v[v] = true;
            } else {
                break;
            }
        }
    }

#ifdef TIME
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "query for baseline: " << elapsed_seconds.count() << endl;
#endif

#ifdef PRINTABCORE
    cout << "upper vertices: " << endl;
    for (auto u = 0; u < node_u.size(); u ++) {
        if (node_u[u]) cout << " " << u;
    }
    cout << endl;

    cout << "lower vertices: " << endl;
    for (auto v = 0; v < node_v.size(); v ++) {
        if (node_v[v]) cout << " " << v;
    }
    cout << endl;
#endif
}
