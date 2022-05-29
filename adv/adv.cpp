//
// Created by lsy on 24/5/22.
//

#include "adv.h"
#include "../config/config.h"


/**
 * compute the core time neighbor,
 * the CT neighbors of u is a hash table where each key is a vertex v
 * (u,v) in G[ts,te], each value of v is the number of edges (u,v,t)
 * with ts <= t <= te
 * @param ts strating time
 * @param cn neighbor
 * @param num end u id
 * @param ctn ctn
 */
auto compute_ctn(const vid_t& ts,  const int& te,
                 const int& alpha, const int& beta, BiGraph& g) -> void {
    // TODO te should be the end time
    for (auto index = g.edges_idx[ts]; index < g.edges_idx[te]; index++) {
        auto u = g.edges[index].first;
        auto v = g.edges[index].second;

        // firstly check u is in abcore, and v is in abcore
        if (g.left_index[u].size() - 1 < alpha) continue;
        if (g.left_index[u][alpha] < beta) continue;
        if (g.right_index[v].size() - 1 < beta) continue;
        if (g.right_index[v][beta] < alpha) continue;

        // then search the neighbor of u and v
        
    }

//    for (auto u = 0; u < u_max; u ++) {
//        ctn[u].clear();
//
//        for (auto index = int(nu[u].size() - 1); index >= 0; index --) {
//            auto v = nu[u][index].first;
//            auto t = nu[u][index].second;
//
//            // because time from largest to smallest
//            if (t < ts) break;
//            ctn[u][v] += 1;
//        }
//    }
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
    // TODO

    for (auto alpha = 2; alpha < g.a_to_b.size(); alpha ++) {
        if (g.a_to_b[alpha] == 0) continue;
        for (auto beta = 1; beta < g.a_to_b[alpha]; beta ++) {
            // then ts from 1 to tmax
            for (auto ts = 0; ts < g.tmax; ++ ts) {
                // TODO then delelte edges
                /*
                 * the delete edge is similary with baseline
                 * need to implement CTN(u)k, need to consider how to implement
                 * and how to manage the core numbers
                 */

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