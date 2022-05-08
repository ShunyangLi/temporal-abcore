//
// Created by lsy on 2022/5/2.
//

#include "tabcore_baseline.h"
#include "../abcore/abcore.h"
#include "../config/config.h"

/**
 * compute the vertex storage index size
 * @param g graph objext
 */
auto tab_vertex_index_size(BiGraph& g) -> void {
    double idx_size = 0;

    // compute the index size for upper index and lower index
    idx_size += double (sizeof(int) * g.tbcore_uindex.size());

    for (auto alpha = 1; alpha < g.tbcore_uindex.size(); alpha ++) {
        idx_size += double (sizeof(int) * g.tbcore_uindex[alpha].size() - 1);
        for (auto beta = 1; beta < g.tbcore_uindex[alpha].size(); beta ++) {
            idx_size += double (sizeof(int) * g.tbcore_uindex[alpha][beta].size() - 1);

            for (auto ts = 0; ts < g.tbcore_uindex[alpha][beta].size(); ts ++) {
                if (g.tbcore_uindex[alpha][beta][ts].empty()) {
                    idx_size -= sizeof(int);
                    break;
                }

                auto block = g.tbcore_uindex[alpha][beta][ts].front();

                while (block != nullptr) {
                    idx_size += sizeof(int);
                    idx_size += sizeof(vertex_block) * 2;
                    idx_size += double (block->nodeset.size() * sizeof(int ));
                    block = block->child;
                }
            }
        }
    }

    idx_size += double (sizeof(int) * g.tbcore_vindex.size());
    for (auto beta = 1; beta < g.tbcore_vindex.size(); beta ++) {
        idx_size += double (sizeof(int) * g.tbcore_vindex[beta].size() - 1);
        for (auto alpha = 1; alpha < g.tbcore_vindex[beta].size(); alpha ++) {
            idx_size += double (sizeof(int) * g.tbcore_vindex[beta][alpha].size());

            for (auto ts = 0; ts < g.tbcore_vindex[beta][alpha].size(); ts ++) {
                if (g.tbcore_vindex[beta][alpha][ts].empty()) {
                    idx_size -= sizeof(int);
                    break;
                }
                auto block = g.tbcore_vindex[beta][alpha][ts].front();
                while (block != nullptr) {
                    idx_size += sizeof(int);
                    idx_size += sizeof(vertex_block) * 2;
                    idx_size += double (block->nodeset.size() * sizeof(int ));
                    block = block->child;
                }
            }
        }
    }

    auto graph_size = double (g.num_edges * 3 * sizeof(int));

#ifdef MBS
    cout << "Graph size: " << graph_size / 1000000 << " MB." << endl;
    cout << "Index size: " << double (idx_size / 1000000) << " MB." << endl;
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
        for (int index = int(nu[u].size() - 1); index >= 0; index --) {
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

    for (int i = int(index[alpha][beta][ts].size() - 1); i >=0 ; --i) {
        if (index[alpha][beta][ts][i]->te == te) {
            index[alpha][beta][ts][i]->nodeset.push_back(u);
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
            auto u = tg.edges[index].first;
            auto v = tg.edges[index].second;

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

            auto vu = vector<bool>(g.num_v1, false);
            auto vv = vector<bool>(g.num_v2, false);
            // then we just check whether the the core number is changed
            for (auto const & u : au)  {

                if (vu[u]) continue;
                else vu[u] = true;

                if (u_alpha_offset[u].size() > tg.left_index[u].size()) {
                    // then record it.
                    for (auto alpha =  u_alpha_offset[u].size() - 1; alpha > tg.left_index[u].size() - 1; --alpha) {
                        for (auto beta = 1; beta <= u_alpha_offset[u][alpha]; beta ++)
                            tab_update_index(g.tbcore_uindex, ts, _te, u, alpha, beta);
                    }
                }

                for (auto alpha  = int(tg.left_index[u].size() - 1); alpha >= 1; --alpha) {
                    if (tg.left_index[u][alpha] != u_alpha_offset[u][alpha]) {
                        for (auto beta = int(u_alpha_offset[u][alpha]); beta > int(tg.left_index[u][alpha]); beta --)
                            tab_update_index(g.tbcore_uindex, ts, _te, u, alpha, beta);
                    }
                }
            }

            for (auto const & v: av) {
                if (vv[v]) continue;
                else vv[v] = true;

                // the alpha value of u becomes smaller
                if (v_beta_offset[v].size() > tg.right_index[v].size()) {
                    // then record it.
                    for (auto beta =  v_beta_offset[v].size() - 1; beta > tg.right_index[v].size() - 1; --beta) {
                        for (auto alpha = 1; alpha <= v_beta_offset[v][beta]; alpha ++)
                            tab_update_index(g.tbcore_vindex, ts, _te, v, beta, alpha);
                    }
                }

                for (auto beta  = int(tg.right_index[v].size() - 1); beta >= 1; beta --) {
                    if (tg.right_index[v][beta] != v_beta_offset[v][beta]) {
                        for (auto alpha = int(v_beta_offset[v][beta]); alpha > tg.right_index[v][beta]; alpha --)
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
auto tab_advance_del_edge(BiGraph& g, const int& _te) -> void {
//    auto tg = g;
    auto q = queue<pair<vid_t, vid_t>>();

    // because there are may one more than one edge that between ts and te
    for (auto index = g.edges_idx[_te]; index < g.edges_idx[_te + 1]; ++index) {
        auto u = g.edges[index].first;
        auto v = g.edges[index].second;

        // then we process u and v
        -- g.ucn[u][v];
        -- g.vcn[v][u];

        if (g.ucn[u][v] == 0) g.ucn[u].erase(v);
        if (g.vcn[v][u] == 0) g.vcn[v].erase(u);

        // when the neighbors of u is less then a, then update
        // it is one edge, so just detect only once
        if (g.ucn[u].size() < g.left_index[u].size() - 1 || g.vcn[v].size() < g.right_index[v].size() - 1) {
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
    auto start = chrono::system_clock::now();

    coreIndexKCore(g);

    g.tbcore_uindex.resize(2);
    g.tbcore_vindex.resize(2);

    // count the neighbor in the time interval ts to te
    tab_compute_core_neighbor(0, g.ucn, g.num_v1, g.tnu);
    tab_compute_core_neighbor(0, g.vcn, g.num_v2, g.tnv);

    // then start peeling
    for (auto ts = 0; ts < g.tmax; ++ ts) {
        // then working here
        if (ts == g.tmax - 1) break;

        auto tg = g;
        tab_back_del_edges(g, tg, ts, g.tmax - 1);

        // count the neighbor in the time interval ts to te
        tab_compute_core_neighbor(ts, g.ucn, g.num_v1, g.tnu);
        tab_compute_core_neighbor(ts, g.vcn, g.num_v2, g.tnv);

        // delete the visited edges
        if (ts < g.tmax) tab_advance_del_edge(g, ts);
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
auto query(const int& alpha, const int& beta, const int& ts, const int& te, BiGraph& g, vector<bool>& node_u,
           vector<bool>& node_v) -> void {
#ifdef TIME
    auto start = chrono::system_clock::now();
#endif
    node_u = vector<bool>(g.num_v1, false);
    node_v = vector<bool>(g.num_v2, false);

    if (alpha > g.tbcore_uindex.size()) return;
    if (beta > g.tbcore_uindex[alpha].size()) return;
    if (ts > g.tbcore_uindex[alpha][beta].size()) return;
    if (g.tbcore_uindex[alpha][beta][ts].empty()) return;

    auto block = g.tbcore_uindex[alpha][beta][ts].front();
    while (block != nullptr && block->te <= te) {
        for (auto const& u: block->nodeset) node_u[u] = true;
        block = block->child;
    }

    block = g.tbcore_vindex[beta][alpha][ts].front();
    while (block != nullptr && block->te <= te) {
        for (auto const& v: block->nodeset) node_v[v] = true;
        block = block->child;
    }

#ifdef TIME
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "query for tabcore: " << elapsed_seconds.count() << endl;
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


/**
 * implement the skip storage of baseline index
 * @param g
 */
auto adv_tabcore_baseline(BiGraph& g) -> void {
    auto start = chrono::system_clock::now();

    coreIndexKCore(g);

    g.tbcore_uindex.resize(2);
    g.tbcore_vindex.resize(2);

    tab_compute_core_neighbor(0, g.ucn, g.num_v1, g.tnu);
    tab_compute_core_neighbor(0, g.vcn, g.num_v2, g.tnv);

    // then start peeling
    for (auto ts = 0; ts < g.tmax; ts += 3) {
        // then working here
        if (ts == g.tmax - 1) break;

        auto tg = g;
        tab_back_del_edges(g, tg, ts, g.tmax - 1);

        // count the neighbor in the time interval ts to te
        tab_compute_core_neighbor(ts, g.ucn, g.num_v1, g.tnu);
        tab_compute_core_neighbor(ts, g.vcn, g.num_v2, g.tnv);

        // delete the visited edges
        for (auto tmp = ts; tmp < ts + 3; tmp ++) {
            if (ts < g.tmax) tab_advance_del_edge(g, tmp);
        }
    }

    cout << "finished skip structure baseline" << endl;

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "construction: " << elapsed_seconds.count() << endl;


#ifdef INDEX_SIZE
    tab_vertex_index_size(g);
#endif
}


/**
 * query for skip structure index
 */
auto query_skip(const int& alpha, const int& beta, const int& ts, const int& te, BiGraph& g, vector<bool>& node_u,
           vector<bool>& node_v) -> void {
#ifdef TIME
    auto start = chrono::system_clock::now();
#endif

    node_u = vector<bool>(g.num_v1, false);
    node_v = vector<bool>(g.num_v2, false);

    if (alpha > g.tbcore_uindex.size()) return;
    if (beta > g.tbcore_uindex[alpha].size()) return;
    if (ts > g.tbcore_uindex[alpha][beta].size()) return;

    // select the last tts that tts <= ts
    auto tts = ts;
    if (ts % 3 != 0) {
        tts = int(ts / 3) * 3 + 3;
    }

    auto block = g.tbcore_uindex[alpha][beta][tts].front();
    while (block != nullptr && block->te <= te) {
        for (auto const& u: block->nodeset) {
            node_u[u] = true;
        }
        block = block->child;
    }

    block = g.tbcore_vindex[beta][alpha][tts].front();
    while (block != nullptr && block->te <= te) {
        for (auto const& v: block->nodeset) node_v[v] = true;
        block = block->child;
    }

    if (ts != tts) {
        // we add the edge to check whether work.
        auto uDegree = map<int,int>();
        auto vDegree = map<int,int>();

        for (auto index = g.edges_idx[ts]; index < g.edges_idx[tts]; index++) {
            auto u = g.edges[index].first;
            auto v = g.edges[index].second;

            if (!node_u[u]) {
                uDegree[u] ++;
                if (uDegree[u] >= alpha) node_u[u]= true;
            }
            if (!node_v[v]) {
                vDegree[v] ++;
                if (vDegree[v] >= beta) node_v[v] = true;
            }
        }
    }


#ifdef TIME
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "query for tabcore: " << elapsed_seconds.count() << endl;
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