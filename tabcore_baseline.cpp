//
// Created by lsy on 2022/5/2.
//

#include "tabcore_baseline.h"
#include "abcore.h"
#include "config.h"

/**
 * compute the vertex storage index size
 * @param g graph objext
 */
auto tab_vertex_index_size(BiGraph& g) -> void {
    double idx_size = 0;

    // compute the index size for upper index and lower index
    idx_size += sizeof(int) * g.tbcore_uindex.size();

    for (auto alpha = 1; alpha < g.tbcore_uindex.size(); alpha ++) {
        idx_size += sizeof(int) * g.tbcore_uindex[alpha].size() - 1;
        for (auto beta = 1; beta < g.tbcore_uindex[alpha].size(); beta ++) {
            idx_size += sizeof(int) * g.tbcore_uindex[alpha][beta].size() - 1;

            for (auto ts = 0; ts < g.tbcore_uindex[alpha][beta].size(); ts ++) {
                if (g.tbcore_uindex[alpha][beta][ts].empty()) {
                    idx_size -= sizeof(int);
                    break;
                }

                auto block = g.tbcore_uindex[alpha][beta][ts].front();

                while (block != nullptr) {
                    idx_size += sizeof(int);
                    idx_size += sizeof(vertex_block) * 2;
                    idx_size += block->nodeset.size() * sizeof(int );
                    block = block->child;
                }
            }
        }
    }

    idx_size += sizeof(int) * g.tbcore_vindex.size();
    for (auto beta = 1; beta < g.tbcore_vindex.size(); beta ++) {
        idx_size += sizeof(int) * g.tbcore_vindex[beta].size() - 1;
        for (auto alpha = 1; alpha < g.tbcore_vindex[beta].size(); alpha ++) {
            idx_size += sizeof(int) * g.tbcore_vindex[beta][alpha].size();

            for (auto ts = 0; ts < g.tbcore_vindex[beta][alpha].size(); ts ++) {
                if (g.tbcore_vindex[beta][alpha][ts].empty()) {
                    idx_size -= sizeof(int);
                    break;
                }
                auto block = g.tbcore_vindex[beta][alpha][ts].front();
                while (block != nullptr) {
                    idx_size += sizeof(int);
                    idx_size += sizeof(vertex_block) * 2;
                    idx_size += block->nodeset.size() * sizeof(int );
                    block = block->child;
                }
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
auto tab_compute_core_neighbor(const vid_t& ts, vector<unordered_map<int, int>>& cn,
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
auto tab_update_index(vector<vector<vector<vector<vertex_block*>>>>& index,
                  const int& ts, const int& te,
                  const vid_t& u, const int& alpha, const int& beta) -> void {

    if (index.size() < alpha + 1) index.resize(alpha + 1);
    if (index[alpha].size() < beta + 1) index[alpha].resize(beta + 1);
    if (index[alpha][beta].size() < ts + 1) index[alpha][beta].resize(ts + 1);

    // if not then just create a new one
    if (index[alpha][beta][ts].empty()) {
        auto block = new vertex_block;
        block->nodeset.push_back(u);
        block->te = te;
        index[alpha][beta][ts].push_back(block);
        return;
    }


    for (int i = index[alpha][beta][ts].size() - 1; i >=0 ; --i) {
        if (index[alpha][beta][ts][i]->te == te) {
            index[alpha][beta][ts][i]->nodeset.push_back(te);
            return;
        }
    }

    // insert it in the front
    auto block = new vertex_block;
    block->nodeset.push_back(u);
    block->te = te;
    block->child = index[alpha][beta][ts].front();
    index[alpha][beta][ts].front()->parent = block;
    index[alpha][beta][ts].insert(index[alpha][beta][ts].begin(), block);

}

/**
 * delete edges for peeling ts to te
 * @param g
 */
auto tab_back_del_edges(BiGraph& g, BiGraph& tg,
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
                        if (alpha > tg.left_index[u].size() - 1) tab_update_index(g.tbcore_uindex, ts, _te, u, alpha, beta);
                        else {
                            if (tg.left_index[u][alpha] != u_alpha_offset[u][alpha])
                                tab_update_index(g.tbcore_uindex, ts, _te, u, alpha, beta);
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
                        tab_update_index(g.tbcore_vindex, ts, _te, v, beta, alpha);
                    }
                }

                for (int beta  = tg.right_index[v].size() - 1; beta >= 1; beta --) {
                    if (tg.right_index[v][beta] != v_beta_offset[v][beta]) {
                        auto alpha = v_beta_offset[v][beta];
                        tab_update_index(g.tbcore_vindex, ts, _te, v, beta, alpha);
                    }
                }
            }
        }
    }
}

/**
 * delete the used edge from front to end
 */
auto tab_advance_del_edge(BiGraph& g, const int& ts, const int& _te) -> void {
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
    }
}

/**
 * baseline for computing abcore
 * @param g
 */
auto tabcore_baseline(BiGraph& g) -> void {
    coreIndexKCore(g);

    g.tbcore_uindex.resize(2);
    g.tbcore_vindex.resize(2);

    auto start = chrono::system_clock::now();
    // then start peeling
    for (auto ts = 0; ts < g.tmax; ++ ts) {
        // then working here
        if (ts == g.tmax - 1) break;

        // count the neighbor in the time interval ts to te
        tab_compute_core_neighbor(ts, g.ucn, g.num_v1, g.tnu);
        tab_compute_core_neighbor(ts, g.vcn, g.num_v2, g.tnv);

        auto tg = g;
        tab_back_del_edges(g, tg, ts, g.tmax - 1);

        // delete the visited edges
        if (ts < g.tmax) tab_advance_del_edge(g, ts, ts);
//        cout << "ts: " << ts << endl;
    }

    cout << "finished temporal abcore baseline" << endl;

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "construction: " << elapsed_seconds.count() << endl;


#ifdef INDEX_SIZE
    tab_vertex_index_size(g);
#endif
}


/**
 * Given a ts, te and alpha beta value, we aim to find the (a,b)-core vertices
 * in G[ts,te]
 */
auto query(const int& ts, const int& te, const int& alpha, const int& beta, vector<bool>& node_u,
           vector<bool>& node_v, BiGraph& g) -> void {

}