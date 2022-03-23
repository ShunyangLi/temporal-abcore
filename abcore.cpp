#include "abcore.h"

using namespace std;

void crossUpdate_for_kcore(BiGraph& g, int alpha, int k_x, vid_t v) {
    for (int beta = k_x; beta > 0; beta--) {
        if (g.right_index[v][beta] < alpha) {
            g.right_index[v][beta] = alpha;
        }
        else {
            break;
        }
    }
}

void alphaCopyPeel_for_kcore(int left_k, BiGraph& g) {
    int dd_;
    int pre_left_k_ = left_k - 1;
    vector<bool> left_deletion_next_round;
    vector<bool> right_deletion_next_round;
    vector<int> left_degree_next_round;
    vector<int> right_degree_next_round;
    vector<vid_t> left_vertices_to_be_peeled;
    vector<vid_t> right_vertices_to_be_peeled;
    for (vid_t u = 0; u < g.getV1Num(); u++) {
        if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
            left_vertices_to_be_peeled.push_back(u);
        }
    }
    int right_remain_nodes_num = g.num_v2;
    vector<vid_t> right_remain_nodes; right_remain_nodes.resize(g.num_v2);
    for (int i = 0; i < right_remain_nodes.size(); i++) {
        right_remain_nodes[i] = i;
    }
    int right_remain_nodes_tmp_num = 0;
    vector<vid_t> right_remain_nodes_tmp; right_remain_nodes_tmp.resize(g.num_v2);
    bool update_flag = false;
    for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
        if (right_k - 1 > 0) {
            update_flag = true;
        }
        int pre_ = right_k - 1;
        bool stop = true;
        right_remain_nodes_tmp_num = 0;
        for (int i = 0; i < right_remain_nodes_num; i++) {
            vid_t v = right_remain_nodes[i];
            if (!g.right_delete[v]) {
                stop = false;
                right_remain_nodes_tmp[right_remain_nodes_tmp_num] = v;
                right_remain_nodes_tmp_num++;
                if (g.degree_v2[v] < right_k) {
                    right_vertices_to_be_peeled.push_back(v);
                }
            }
        }
        swap(right_remain_nodes, right_remain_nodes_tmp);
        right_remain_nodes_num = right_remain_nodes_tmp_num;
        if (stop) break;
        while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
            // peel left
            int oo_ = left_vertices_to_be_peeled.size();
            for (int j = 0; j < oo_; j++) {
                vid_t u = left_vertices_to_be_peeled[j];
                if (g.left_delete[u]) continue;
                vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
                int ss = tmp_neigh_.size();
                for (int k = 0; k < ss; k++) {
                    vid_t v = tmp_neigh_[k];
                    if (g.right_delete[v]) continue;
                    dd_ = --g.degree_v2[v];
                    if (update_flag && dd_ == 0) {
                        crossUpdate_for_kcore(g, left_k, pre_, v);
                        g.right_delete[v] = true;
                    }
                    if (dd_ == pre_) {
                        right_vertices_to_be_peeled.push_back(v);
                    }
                }
                g.degree_v1[u] = 0;
                g.left_delete[u] = true;
                if (update_flag) {
                    g.left_index[u][left_k] = pre_;
                }
            }
            left_vertices_to_be_peeled.clear();
            // peel right
            oo_ = right_vertices_to_be_peeled.size();
            for (int j = 0; j < oo_; j++) {
                vid_t v = right_vertices_to_be_peeled[j];
                if (g.right_delete[v]) continue;
                vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v];
                int ss = tmp_neigh_.size();
                for (int k = 0; k < ss; k++) {
                    vid_t u = tmp_neigh_[k];
                    if (g.left_delete[u]) continue;
                    dd_ = --g.degree_v1[u];
                    if (update_flag && dd_ == 0) {
                        g.left_index[u][left_k] = pre_;
                        g.left_delete[u] = true;
                    }
                    if (dd_ == pre_left_k_) {
                        left_vertices_to_be_peeled.push_back(u);
                    }
                }
                g.degree_v2[v] = 0;
                g.right_delete[v] = true;
                if (update_flag) {
                    crossUpdate_for_kcore(g, left_k, pre_, v);
                }
            }
            right_vertices_to_be_peeled.clear();
        }
        if (right_k == 1) {
            left_degree_next_round = g.degree_v1;
            right_degree_next_round = g.degree_v2;
            left_deletion_next_round = g.left_delete;
            right_deletion_next_round = g.right_delete;
        }
    }
    g.degree_v1 = left_degree_next_round;
    g.degree_v2 = right_degree_next_round;
    g.left_delete = left_deletion_next_round;
    g.right_delete = right_deletion_next_round;
    g.v1_max_degree = 0;
    g.v2_max_degree = 0;
    for (vid_t u = 0; u < g.degree_v1.size(); u++) {
        if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
    }
    for (vid_t v = 0; v < g.degree_v2.size(); v++) {
        if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
    }
}

int coreIndexKCore(BiGraph& g) {
    int left_degree_max = 0;
    for (int i = 0; i < g.getV1Num(); i++) {
        if (left_degree_max < g.getV1Degree(i)) left_degree_max = g.getV1Degree(i);
    }
    int right_degree_max = 0;
    for (int i = 0; i < g.getV2Num(); i++) {
        if (right_degree_max < g.getV2Degree(i)) right_degree_max = g.getV2Degree(i);
    }
    // init g's max degree and index
    g.v1_max_degree = left_degree_max;
    g.v2_max_degree = right_degree_max;
    g.left_index.resize(g.getV1Num());
    g.right_index.resize(g.getV2Num());
    g.left_delete.resize(g.getV1Num());
    g.right_delete.resize(g.getV2Num());
    fill_n(g.left_delete.begin(), g.left_delete.size(), false);
    fill_n(g.right_delete.begin(), g.right_delete.size(), false);
    for (int i = 0; i < g.getV1Num(); i++) {
        g.left_index[i].resize(g.getV1Degree(i) + 1);
        fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
    }
    for (int i = 0; i < g.getV2Num(); i++) {
        g.right_index[i].resize(g.getV2Degree(i) + 1);
        fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
    }
    int beta_s = 0;
    for (int left_k = 1; left_k <= g.v1_max_degree; left_k++) {
        alphaCopyPeel_for_kcore(left_k, g);
        beta_s = 0;
        for (vid_t u = 0; u < g.num_v1; u++) {
            if (g.degree_v1[u] <= left_k) continue;
            int right_k = g.left_index[u][left_k];
            if (beta_s < right_k) beta_s = right_k;
        }
        if (beta_s <= left_k) break;
    }
    // restore g
    fill_n(g.left_delete.begin(), g.left_delete.size(), false);
    g.v1_max_degree = left_degree_max;
    for (vid_t u = 0; u < g.num_v1; u++) {
        g.degree_v1[u] = g.neighbor_v1[u].size();
    }
    fill_n(g.right_delete.begin(), g.right_delete.size(), false);
    g.v2_max_degree = right_degree_max;
    for (vid_t v = 0; v < g.num_v2; v++) {
        g.degree_v2[v] = g.neighbor_v2[v].size();
    }
    inv(g);
    for (int left_k = 1; left_k <= beta_s; left_k++) {
        alphaCopyPeel_for_kcore(left_k, g);
    }
    inv(g);
    // restore g
    fill_n(g.left_delete.begin(), g.left_delete.size(), false);
    g.v1_max_degree = left_degree_max;
    for (vid_t u = 0; u < g.num_v1; u++) {
        g.degree_v1[u] = g.neighbor_v1[u].size();
    }
    fill_n(g.right_delete.begin(), g.right_delete.size(), false);
    g.v2_max_degree = right_degree_max;
    for (vid_t v = 0; v < g.num_v2; v++) {
        g.degree_v2[v] = g.neighbor_v2[v].size();
    }
    return beta_s;
}

///////////////////////////
/// dynamic maintenance ///
///////////////////////////

void dyn_crossUpdate_addition(BiGraph& g, int alpha, int k_x, vid_t v) {
    for (int beta = k_x; beta > 0; beta--) {
        int oldalpha = g.right_index[v][beta];
        if (oldalpha < alpha) {
            g.right_index[v][beta] = alpha;
        }
        else {
            break;
        }
    }
}
void dyn_crossUpdate_deletion(BiGraph& g, int alpha, int k_x, vid_t v) {
    int newalpha = alpha - 1;
    int truedegree = g.neighbor_v2[v].size();
    for (int i = k_x; i <= truedegree; i++) {
        int oldalpha = g.right_index[v][i];
        if (oldalpha > newalpha) {
            g.right_index[v][i] = newalpha;
        }
        else {
            break;
        }
    }
}

void compute_a_b_core(BiGraph& g, int alpha, int beta) {
    vector<vid_t> left_vertices_to_be_peeled;
    vector<vid_t> right_vertices_to_be_peeled;
    bool stop = true;
    for (vid_t u = 0; u < g.num_v1; u++) {
        if (!g.left_delete[u]) {
            if (g.degree_v1[u] < alpha) {
                left_vertices_to_be_peeled.push_back(u);
            }
        }
    }
    for (vid_t v = 0; v < g.num_v2; v++) {
        if (!g.right_delete[v]) {
            if (g.degree_v2[v] < beta) {
                right_vertices_to_be_peeled.push_back(v);
            }
        }
    }
    while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
        // peel left
        int oo_ = left_vertices_to_be_peeled.size();
        for (int j = 0; j < oo_; j++) {
            vid_t u = left_vertices_to_be_peeled[j];
            if (g.left_delete[u]) continue;
            vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
            int ss = tmp_neigh_.size();
            for (int k = 0; k < ss; k++) {
                vid_t v = tmp_neigh_[k];
                if (g.right_delete[v]) continue;
                int dd_ = --g.degree_v2[v];
                if (dd_ == 0) {
                    // core part
                    g.right_delete[v] = true;
                }
                if (dd_ == beta - 1) {
                    right_vertices_to_be_peeled.push_back(v);
                }
            }
            g.degree_v1[u] = 0;
            g.left_delete[u] = true;
        }
        left_vertices_to_be_peeled.clear();
        // peel right
        oo_ = right_vertices_to_be_peeled.size();
        for (int j = 0; j < oo_; j++) {
            vid_t v = right_vertices_to_be_peeled[j];
            if (g.right_delete[v]) continue;
            vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v];
            int ss = tmp_neigh_.size();
            for (int k = 0; k < ss; k++) {
                vid_t u = tmp_neigh_[k];
                if (g.left_delete[u]) continue;
                int dd_ = --g.degree_v1[u];
                if (dd_ == 0) {
                    g.left_delete[u] = true;
                }
                if (dd_ == alpha - 1) {
                    left_vertices_to_be_peeled.push_back(u);
                }
            }
            g.degree_v2[v] = 0;
            g.right_delete[v] = true;
        }
        right_vertices_to_be_peeled.clear();
    }
}

void update_index_with_fixed_left_k_deletion_with_limit_swap(BiGraph& g, int alpha, int tau_alpha, int start_bound, vid_t u, vid_t v) {
    int beta = tau_alpha;

    compute_a_b_core(g, alpha, beta);

    for (vid_t u = 0; u < g.num_v1; u++) {
        if (g.left_delete[u] && g.left_index[u].size() >= alpha + 1) {
            int oldbeta = g.left_index[u][alpha];
            if (oldbeta == beta) {
                g.left_index[u][alpha] = beta - 1;
            }
        }
    }
    for (vid_t v = 0; v < g.num_v2; v++) {
        if (g.right_delete[v] && g.right_index[v].size() >= beta + 1) {
            dyn_crossUpdate_deletion(g, alpha, beta, v);
        }
    }
    // contrary to insertion
    int bound = start_bound;
    if (g.left_index[u][alpha] > bound) {
        g.left_index[u][alpha] = bound;
    }
}

pair<int, int> calculate_Delta(BiGraph& g) {
    int delta = 0;
    for (vid_t u = 0; u < g.num_v1; u++) {
        if (g.left_index[u].size() > delta + 1) {
            for (int alpha = delta + 1; alpha < g.left_index[u].size(); alpha++) {
                if (g.left_index[u][alpha] >= alpha) {
                    delta = alpha;
                }
                else {
                    break;
                }
            }
        }
    }
    if (delta == 0) {
        return make_pair(1, 0);
    }
    else {
        return make_pair(delta, delta);
    }
}

int calculate_vbeta_with_fixed_alpha(BiGraph& g, vid_t v, int alpha) {
    int ss = g.right_index[v].size();
    for (int beta = 1; beta < ss; beta++) {
        if (g.right_index[v][beta] < alpha) {
            return beta - 1;
        }
    }
    return ss - 1;
}

int calculate_bound_for_left_node(BiGraph& g, vid_t u, int alpha) {
    // deletion operation when removing level
    if (alpha > g.degree_v1[u]) {
        return 0;
    }
    vector<int> k_index;
    vector<vid_t>& neigh = g.neighbor_v1[u];
    int ss = neigh.size();
    for (int i = 0; i < ss; i++) {
        int vbeta = calculate_vbeta_with_fixed_alpha(g, neigh[i], alpha);
        k_index.push_back(vbeta);
    }
    sort(k_index.begin(), k_index.end());
    return k_index[k_index.size() - alpha];
}

int calculate_ualpha_with_fixed_beta(BiGraph& g, vid_t u, int beta) {
    int ss = g.left_index[u].size();
    for (int alpha = 1; alpha < ss; alpha++) {
        if (g.left_index[u][alpha] < beta) {
            return alpha - 1;
        }
    }
    return ss - 1;
}

int calculate_bound_for_right_node(BiGraph& g, vid_t v, int beta) {
    // deletion operation when removing level
    if (beta > g.degree_v2[v]) {
        return 0;
    }
    vector<int> k_index;
    vector<vid_t>& neigh = g.neighbor_v2[v];
    int ss = neigh.size();
    for (int i = 0; i < ss; i++) {
        int ualpha = calculate_ualpha_with_fixed_beta(g, neigh[i], beta);
        k_index.push_back(ualpha);
    }
    sort(k_index.begin(), k_index.end());
    return k_index[k_index.size() - beta];
}

void update_index_with_fixed_left_k_addition_with_limit_swap(BiGraph& g, int alpha, int tau_alpha, vid_t u, vid_t v) {
    if (tau_alpha > g.left_index[u][alpha]) {
        g.left_index[u][alpha] = tau_alpha;
    }
    int beta = tau_alpha;
    // compute alpha-beta+1-core
    compute_a_b_core(g, alpha, beta + 1);
    // update core value
    for (vid_t u = 0; u < g.num_v1; u++) {
        if (!g.left_delete[u]) {
            int oldbeta = g.left_index[u][alpha];
            if (oldbeta < beta + 1) {
                g.left_index[u][alpha] = beta + 1;
            }
        }
    }
    for (vid_t v = 0; v < g.num_v2; v++) {
        if (!g.right_delete[v]) {
            dyn_crossUpdate_addition(g, alpha, beta + 1, v);
        }
    }
}

double update_bicore_index(BiGraph& g, vid_t u, vid_t v, bool addition)
{
    auto start = chrono::system_clock::now();
    // calculate swap threshold
    pair<int, int> r = calculate_Delta(g);
    // check whether operation is legal and insert/remove the edge
    if (addition) {
        if (u < g.num_v1 && v < g.num_v2) {
            vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
            int ss = tmp_neigh_.size();
            for (int i = 0; i < ss; i++) {
                if (tmp_neigh_[i] == v) {
                    cout << "illegal insertion 1" << endl;
                    return 0;
                }
            }
            g.addEdge(u, v);
        }
        else {
            cout << "illegal insertion 2" << endl;
            return 0;
        }
    }
    else {
        if (u < g.num_v1 && v < g.num_v2) {
            vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
            int ss = tmp_neigh_.size();
            bool deleted = false;
            for (int i = 0; i < ss; i++) {
                if (tmp_neigh_[i] == v) {
                    deleted = true;
                    g.deleteEdge(u, v);
                    break;
                }
            }
            if (!deleted) {
                cout << "illegal deletion 1" << endl;
                return 0;
            }
        }
        else {
            cout << "illegal deletion 2" << endl;
            return 0;
        }
    }
    // addition
    if (addition) {
        // extend core index and bicore index and set threshold
        g.left_index[u].push_back(0);
        g.right_index[v].push_back(0);
        int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
        int th_alpha, th_beta;
        if (r.first + 1 + r.second + 1 >= max_alpha) {
            th_alpha = max_alpha; th_beta = 0;
        }
        else if (r.first + 1 + r.second + 1 >= max_beta) {
            th_alpha = 0; th_beta = max_beta;
        }
        else {
            th_alpha = r.first + 1; th_beta = r.second + 1;
        }
        // compute alpha beta condition in advance
        vector<int> tau_alpha; vector<int> tau_beta;
        tau_alpha.push_back(-1); tau_beta.push_back(-1);
        for (int alpha = 1; alpha <= th_alpha; alpha++) {
            int ubeta = g.left_index[u][alpha];
            int bound = calculate_bound_for_left_node(g, u, alpha);
            int vbeta = -1;
            bool changed = false;
            for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
                if (g.right_index[v][b] < alpha) {
                    vbeta = b - 1;
                    changed = true;
                    break;
                }
            }
            if (!changed) {
                vbeta = g.degree_v2[v];
            }
            if (vbeta < 0) {
                cout << "error: vbeta" << endl;
            }
            int beta = bound > vbeta ? vbeta : bound;
            tau_alpha.push_back(beta);
        }
        for (int beta = 1; beta <= th_beta; beta++) {
            int valpha = g.right_index[v][beta];
            int bound = calculate_bound_for_right_node(g, v, beta);
            int ualpha = -1;
            bool changed = false;
            for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
                if (g.left_index[u][a] < beta) {
                    ualpha = a - 1;
                    changed = true;
                    break;
                }
            }
            if (!changed) {
                ualpha = g.degree_v1[u];
            }
            if (ualpha < 0) {
                cout << "error: ualpha" << endl;
            }
            int alpha = bound > ualpha ? ualpha : bound;
            tau_beta.push_back(alpha);
        }
        for (int alpha = 1; alpha <= th_alpha; alpha++) {
            update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha, tau_alpha[alpha], u, v);
            // restore graph
            int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
            fill_n(g.left_delete.begin(), g.left_delete.size(), false);
            for (vid_t u = 0; u < g.num_v1; u++) {
                g.degree_v1[u] = g.neighbor_v1[u].size();
                left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
            }
            g.v1_max_degree = left_degree_max;
            fill_n(g.right_delete.begin(), g.right_delete.size(), false);
            for (vid_t v = 0; v < g.num_v2; v++) {
                g.degree_v2[v] = g.neighbor_v2[v].size();
                right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
            }
            g.v2_max_degree = right_degree_max;
        }
        inv(g);
        for (int alpha = 1; alpha <= th_beta; alpha++) {
            update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha, tau_beta[alpha], v, u);
            // restore graph
            int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
            fill_n(g.left_delete.begin(), g.left_delete.size(), false);
            for (vid_t u = 0; u < g.num_v1; u++) {
                g.degree_v1[u] = g.neighbor_v1[u].size();
                left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
            }
            g.v1_max_degree = left_degree_max;
            fill_n(g.right_delete.begin(), g.right_delete.size(), false);
            for (vid_t v = 0; v < g.num_v2; v++) {
                g.degree_v2[v] = g.neighbor_v2[v].size();
                right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
            }
            g.v2_max_degree = right_degree_max;
        }
        inv(g);
    }
        // deletion
    else {
        int oldusat = g.left_index[u].back();
        int oldvsat = g.right_index[v].back();
        int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
        int th_alpha, th_beta;
        if (r.first + r.second >= max_alpha) {
            th_alpha = max_alpha; th_beta = 0;
        }
        else if (r.first + r.second >= max_beta) {
            th_alpha = 0; th_beta = max_beta;
        }
        else {
            th_alpha = r.first; th_beta = r.second;
        }
        // compute alpha beta condition in advance
        vector<int> alpha_start; vector<int> beta_start; vector<int> alpha_start_bound; vector<int> beta_start_bound;
        alpha_start.push_back(-1); beta_start.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
        for (int alpha = 1; alpha <= th_alpha; alpha++) {
            int ubeta = g.left_index[u][alpha];
            int bound = calculate_bound_for_left_node(g, u, alpha);
            int vbeta = -1;
            bool changed = false;
            // g.degree_v2[v] has already been decreased
            for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
                if (g.right_index[v][b] < alpha) {
                    vbeta = b - 1;
                    changed = true;
                    break;
                }
            }
            if (!changed) {
                vbeta = g.degree_v2[v] + 1;
            }
            if (vbeta == -1 || vbeta == 0) {
                cout << "error: vbeta" << endl;
            }
            int beta = ubeta > vbeta ? vbeta : ubeta;
            alpha_start.push_back(beta);
            alpha_start_bound.push_back(bound);
        }
        for (int beta = 1; beta <= th_beta; beta++) {
            int valpha = g.right_index[v][beta];
            int bound = calculate_bound_for_right_node(g, v, beta);
            int ualpha = -1;
            bool changed = false;
            for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
                if (g.left_index[u][a] < beta) {
                    ualpha = a - 1;
                    changed = true;
                    break;
                }
            }
            if (!changed) {
                ualpha = g.degree_v1[u] + 1;
            }
            if (ualpha == -1 || ualpha == 0) {
                cout << "error: ualpha" << endl;
            }
            int alpha = valpha > ualpha ? ualpha : valpha;
            beta_start.push_back(alpha);
            beta_start_bound.push_back(bound);
        }
        for (int alpha = 1; alpha <= th_alpha; alpha++) {
            update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha, alpha_start[alpha], alpha_start_bound[alpha], u, v);
            // restore graph
            int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
            fill_n(g.left_delete.begin(), g.left_delete.size(), false);
            for (vid_t u = 0; u < g.num_v1; u++) {
                g.degree_v1[u] = g.neighbor_v1[u].size();
                left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
            }
            g.v1_max_degree = left_degree_max;
            fill_n(g.right_delete.begin(), g.right_delete.size(), false);
            for (vid_t v = 0; v < g.num_v2; v++) {
                g.degree_v2[v] = g.neighbor_v2[v].size();
                right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
            }
            g.v2_max_degree = right_degree_max;
        }
        inv(g);
        for (int alpha = 1; alpha <= th_beta; alpha++) {
            update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha, beta_start[alpha], beta_start_bound[alpha], v, u);
            // restore graph
            int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
            fill_n(g.left_delete.begin(), g.left_delete.size(), false);
            for (vid_t u = 0; u < g.num_v1; u++) {
                g.degree_v1[u] = g.neighbor_v1[u].size();
                left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
            }
            g.v1_max_degree = left_degree_max;
            fill_n(g.right_delete.begin(), g.right_delete.size(), false);
            for (vid_t v = 0; v < g.num_v2; v++) {
                g.degree_v2[v] = g.neighbor_v2[v].size();
                right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
            }
            g.v2_max_degree = right_degree_max;
        }
        inv(g);
        g.left_index[u].pop_back();
        g.right_index[v].pop_back();
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    return elapsed_seconds.count();
}