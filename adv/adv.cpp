//
// Created by lsy on 24/5/22.
//

#include "adv.h"
#include "../config/config.h"


/**
 * compute the core time neighbor,
 * the CT neighbors of u is a hash table where each key is a vertex v
 * (u,v) in G[ts,te], each value of v is the number of edges (u,v,t)
 * with ts <= t <= te, for u, core(v) >= core(u)
 * @param ts strating time
 * @param cn neighbor
 * @param num end u id
 * @param ctn ctn
 */
auto compute_ctn(const vid_t& ts,  const int& te,
                 const int& alpha, const int& beta, BiGraph& g) -> void {
    for (auto index = g.edges_idx[ts]; index < g.edges_idx[te]; index++) {
        auto u = g.edges[index].first;
        auto v = g.edges[index].second;

        // firstly check u is in abcore, and v is in abcore
        if (g.left_index[u].size() - 1 < alpha) continue;
        if (g.left_index[u][alpha] < beta) continue;
        if (g.right_index[v].size() - 1 < beta) continue;
        if (g.right_index[v][beta] < alpha) continue;

        // then search the neighbor of u and v
        if (g.right_index[v][beta] >= alpha) {
            g.ucn[u][v] += 1;
        }

        if (g.left_index[u][alpha] >= beta) {
            g.vcn[v][u] += 1;
        }
    }
}

auto adv_del_edges(BiGraph& g, const int& ts, const int& te,
                   const int& alpha, const int& beta) -> void {

    auto qu = queue<vid_t>();
    auto qv = queue<vid_t>();

    for (auto _te = te; _te > ts; _te --) {
        // try to remove the current edges, front to end
        auto visited_u = vector<bool>(g.num_v1);
        auto visited_v = vector<bool>(g.num_v2);

        for (auto index = g.edges_idx[_te]; index < g.edges_idx[_te + 1]; ++index) {
            auto u = g.edges[index].first;
            auto v = g.edges[index].second;

            if (g.left_index[u].size() - 1 < alpha) continue;
            if (g.left_index[u][alpha] < beta) continue;
            if (g.right_index[v].size() - 1 < beta) continue;
            if (g.right_index[v][beta] < alpha) continue;

            // then we process u and v
            -- g.ucn[u][v];
            -- g.vcn[v][u];

            if (g.ucn[u][v] == 0) g.ucn[u].erase(v);
            if (g.vcn[v][u] == 0) g.vcn[v].erase(u);

            // then add it into the quenu to compute again
            if (g.ucn[u].size() < g.left_index[u].size() - 1 ) {
                if (!visited_u[u]) {
                    visited_u[u] = true;
                    qu.push(u);
                }

                for (auto const& it : g.ucn[u]) {
                    if (it.second > te) break;
                    auto tv = it.first;
                    if (g.right_index[tv].size() - 1 >= beta) {
                        if (g.right_index[tv][beta] >= alpha) {
                            if (visited_v[tv]) continue;
                            visited_v[tv] = true;
                            qv.push(tv);
                        }
                    }
                }
            }

            // then process lower vertices
            if (g.vcn[v].size() < g.right_index[v].size() - 1) {
                if(!visited_v[v]) {
                    visited_v[v] = true;
                    qv.push(v);
                }
                for (auto const& it: g.vcn[v]) {
                    if (it.second > te) break;
                    auto tu = it.first;

                    if (g.left_index[tu].size() - 1 >= alpha) {
                        if (g.left_index[tu][alpha] >= beta) {
                            if (visited_u[tu]) continue;
                            visited_u[tu] = true;
                            qu.push(tu);
                        }
                    }
                }
            }
        }

        // for upper just ensure alphe
        while (!qu.empty()) {
            auto u = qu.front();
            qu.pop();
            visited_u[u] = false;

            auto nbr_t = vector<vid_t >();
            auto visited = vector<bool>(g.num_v2);

            auto core_time = 0;
            for (auto const& it : g.ucn[u]) {
                if (it.second < ts) continue;
                if (nbr_t.size() >= alpha && it.second > core_time) break;

                auto v = it.first;

                if (g.right_index[v].size() - 1 < beta ||
                    g.right_index[v][beta] < alpha || visited[v]) continue;

                visited[v] = true;

                auto v_ct = g.v_index[v][beta][alpha].back().second;
                nbr_t.push_back(max(it.second, v_ct));

                if (nbr_t.size() <= alpha) core_time = max(core_time, v_ct);

                // no more process
                if (nbr_t.size() < alpha) {
                    visited_u[u] = true;
                } else {
                    nth_element(nbr_t.begin(),nbr_t.begin()+alpha-1,nbr_t.end());
                    cout << "ts: " << ts << " te: " << nbr_t[alpha-1] << endl;
                }


                // then process the neighbot

            }

        }

        // for lower just ensure beta
        while (!qv.empty()) {

        }

    }

}

/**
 * implement the front to end algorithms
 * based on the core values
 * @param g
 */
auto tabcore_adv(BiGraph& g) -> void {
    auto start = chrono::system_clock::now();

    coreIndexKCore(g);

    // firstly into the core number from 1 - tmax
    compute_core_neighbor(0, g.ucn, g.num_v1, g.tnu);
    compute_core_neighbor(0, g.vcn, g.num_v2, g.tnv);
    // compute the core number from 1-tmax
    auto tg = g;
    back_del_edges(g, tg, 0, g.tmax);

    for (auto alpha = 2; alpha < g.a_to_b.size(); alpha ++) {
        if (g.a_to_b[alpha] == 0) continue;

        for (auto beta = 1; beta < g.a_to_b[alpha]; beta ++) {
            compute_ctn(0, g.tmax, alpha, beta, g);
            // TODO then delelte edges
            /*
             * the delete edge is similary with baseline
             * need to implement CTN(u)k, need to consider how to implement
             * and how to manage the core numbers
             */

            // then ts from 1 to tmax
            for (auto ts = 1; ts < g.tmax; ++ ts) {

            }
        }
    }

#ifdef TIME
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "construction: " << elapsed_seconds.count() << endl;
#endif

#ifdef INDEX_SIZE
    tab_vertex_index_size(g);
#endif

}