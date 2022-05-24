//
// Created by lsy on 24/5/22.
//

#include "adv.h"
#include "../config/config.h"

/**
 * implement the front to end algorithms
 * based on the core values
 * @param g
 */
auto tabcore_adv(BiGraph& g) -> void {
    auto start = chrono::system_clock::now();

    coreIndexKCore(g);

    for (auto alpha = 1; alpha < g.a_to_b.size(); alpha ++) {
        if (g.a_to_b[alpha] == 0) continue;
        for (auto beta = 1; beta < g.a_to_b[alpha]; beta ++) {
            // then ts from 1 to tmax

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