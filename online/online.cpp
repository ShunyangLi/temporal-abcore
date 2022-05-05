//
// Created by lsy on 5/5/22.
//

#include "online.h"


auto compute_a_b_core_nodes(vector<vector<vid_t>>& nu, vector<vector<vid_t>>& nv,
                            int alpha, int beta,
                            vector<bool>& nodes_u,
                            vector<bool>& nodes_v) -> void {

    auto num_v1 = nu.size();
    auto num_v2 = nv.size();
    vector<bool > left_node;
    vector<bool> right_node;
    left_node.resize(num_v1);
    right_node.resize(num_v2);
    fill_n(left_node.begin(), left_node.size(), true);
    fill_n(right_node.begin(), right_node.size(), true);
    vector<vid_t> valid_left_vertex;
    vector<vid_t> valid_right_vertex;
    auto total_visited_edge = 0;
    int left_degree_max = 0;
    for (int i = 0; i < num_v1; i++) {
        if (left_degree_max < nu[i].size()) left_degree_max = nu[i].size();
    }
    int right_degree_max = 0;
    for (int i = 0; i < num_v2; i++) {
        if (right_degree_max < nv[i].size()) right_degree_max = nv[i].size();
    }
    if (left_degree_max < alpha || right_degree_max < beta) return;

    int visited_num = 0;
    int visited_num_bound = 100;

    auto start1 = chrono::steady_clock::now();

    auto degree_v1 = vector<int>(num_v1);
    auto degree_v2 = vector<int>(num_v2);

    for (auto i = 0; i < nu.size(); i++) degree_v1[i] = nu[i].size();
    for (auto i = 0; i < nv.size(); i++) degree_v2[i] = nv[i].size();

    valid_left_vertex.clear();
    valid_right_vertex.clear();
    fill_n(left_node.begin(), left_node.size(), true);
    fill_n(right_node.begin(), right_node.size(), true);

    queue <vid_t> connect_v1;
    queue <vid_t> connect_v2;
    for (int i = 0; i < degree_v1.size(); i++) {
        if (degree_v1[i] < alpha) {
            connect_v1.push(i);
            left_node[i] = false;
        }
    }
    for (int i = 0; i < degree_v2.size(); i++) {
        if (degree_v2[i] < beta) {
            connect_v2.push(i);
            right_node[i] = false;
        }
    }
    while (!connect_v1.empty() || !connect_v2.empty()) {
        while (!connect_v1.empty()) {
            vid_t qq = connect_v1.front();
            connect_v1.pop();
            total_visited_edge += nu[qq].size();
            for (auto &v : nu[qq]) {
                if (right_node[v]) {
                    degree_v2[v]--;
                    if (degree_v2[v] < beta) {
                        connect_v2.push(v);
                        right_node[v] = false;
                    }
                }
            }
        }
        while (!connect_v2.empty()) {
            vid_t qq = connect_v2.front();
            //right_node_list.push_back(qq);
            connect_v2.pop();
            for (auto & u : nv[qq]) {
                if (left_node[u]) {
                    degree_v1[u]--;
                    if (degree_v1[u] < alpha) {
                        connect_v1.push(u);
                        left_node[u] = false;
                    }
                }
            }
        }
    }

    for (int i = 0; i < left_node.size(); i++) {
        if (left_node[i]) nodes_u[i] = true;
    }
    for (int i = 0; i < right_node.size(); i++) {
        if (right_node[i]) nodes_v[i] = true;
    }
}


auto online_peeling (const int& alpha, const int& beta, const int& ts, const int& te,
                     BiGraph& g, vector<bool>& node_u, vector<bool>& node_v) -> void {
#ifdef TIME
    auto start = chrono::system_clock::now();
#endif

    auto nu = vector<vector<vid_t>>(g.num_v1);
    auto nv = vector<vector<vid_t>>(g.num_v2);

    for (auto u = 0; u < g.tnu.size(); u++) {
        for (auto const& e : g.tnu[u]) {
            auto v = e.first;
            auto t = e.second;
            if (ts <= t && te >= t) {
                nu[u].push_back(v);
                break;
            }
        }
    }

    for (auto v = 0; v < g.tnv.size(); v++) {
        for (auto const& e : g.tnv[v]) {
            auto u = e.first;
            auto t = e.second;
            if (ts <= t && te >= t) {
                nv[v].push_back(u);
                break;
            }
        }
    }

    auto nodes_u = vector<bool>(nu.size());
    auto nodes_v = vector<bool>(nv.size());
    compute_a_b_core_nodes(nu, nv, alpha, beta, nodes_u, nodes_v);

#ifdef TIME
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;

    cout << "online query: " << elapsed_seconds.count() << endl;
#endif
}