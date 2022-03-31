#include "bigraph.h"

using namespace std;

BiGraph::BiGraph(const string& dir) {
    num_v1 = 0;
    num_v2 = 0;
    num_edges = 0;

    neighbor_v1.clear();
    neighbor_v2.clear();


    degree_v1.clear();
    degree_v2.clear();

    core_v1.clear();
    core_v2.clear();

    //KKCore index left (x,*) right (*,x)
    left_index.clear();
    right_index.clear();
    v1_max_degree = 0;
    v2_max_degree = 0;
    delta = -1;
    this->dir = dir;
    loadGraph(dir);
}

BiGraph::BiGraph() {
    dir = "";
    num_v1 = 0;
    num_v2 = 0;
    num_edges = 0;

    neighbor_v1.clear();
    neighbor_v2.clear();


    degree_v1.clear();
    degree_v2.clear();

    core_v1.clear();
    core_v2.clear();

    //KKCore index left (x,*) right (*,x)
    left_index.clear();
    right_index.clear();
    v1_max_degree = 0;
    v2_max_degree = 0;
    delta = -1;
}

void BiGraph::init(unsigned int num1, unsigned int num2) {
    num_v1 = num1;
    num_v2 = num2;
    num_edges = 0;

    neighbor_v1.resize(num_v1);
    neighbor_v2.resize(num_v2);
    neighborHash_v1.resize(num1);
    neighborHash_v2.resize(num2);

    degree_v1.resize(num_v1);
    degree_v2.resize(num_v2);

    fill_n(degree_v1.begin(), num_v1, 0);
    fill_n(degree_v2.begin(), num_v2, 0);

    left_delete.resize(num_v1);
    right_delete.resize(num_v2);

    tnu.resize(num1);
    tnv.resize(num2);
}

void BiGraph::loadGraph(const string& data_graph) {
    int u, v, r, t;

    FILE * edgeGraph = fopen(data_graph.c_str(), "r");

    while ((r = fscanf(edgeGraph, "%d %d %d", &u, &v, &t)) != EOF) {
        if (r != 3) {
            fprintf(stderr, "Bad file format: u v t incorrect\n");
            exit(1);
        }

        if (num_edges == 0) {
            fprintf(stdout, "upper vertices: %d, lower vertices: %d, edges: %ld\n", u, v, t);
            init(u, v);
            num_edges = t;
            continue;
        }
        addEdge(u, v);
        addEdgeT(u, v, t);
    }

    edges_idx.emplace_back(edges.size());

    fclose(edgeGraph);

    for (int i = 0; i < num_v1; ++i) {
        neighbor_v1[i].clear();
        neighbor_v1[i] = vector<vid_t>{neighborHash_v1[i].begin(), neighborHash_v1[i].end()};
    }

    for (int i = 0; i < num_v2; ++i) {
        neighbor_v2[i].clear();
        neighbor_v2[i] = vector<vid_t>{neighborHash_v2[i].begin(), neighborHash_v2[i].end()};
    }

    neighborHash_v1.clear();
    neighborHash_v2.clear();
    tmax = time_new_to_old.size() - 1;
    ucn.resize(num_v1);
    vcn.resize(num_v2);
}

void BiGraph::addEdge(vid_t u, vid_t v) {

    // timestamp only use delete
//    neighbor_v1[u].push_back(v);
//    ++degree_v1[u];
//    if (degree_v1[u] > v1_max_degree) v1_max_degree = degree_v1[u];
//    neighbor_v2[v].push_back(u);
//    ++degree_v2[v];
//    if (degree_v2[v] > v2_max_degree) v2_max_degree = degree_v2[v];
//    num_edges++;

    neighborHash_v1[u].insert(v);
    degree_v1[u] = neighborHash_v1[u].size();
    if (degree_v1[u] > v1_max_degree) v1_max_degree = degree_v1[u];
    neighborHash_v2[v].insert(u);
    degree_v2[v] = neighborHash_v2[v].size();
    if (degree_v2[v] > v2_max_degree) v2_max_degree = degree_v2[v];
}

void BiGraph::addEdgeT(vid_t u, vid_t v, vid_t t) {
    edges.emplace_back(std::make_pair(u,v));
    if (time_new_to_old.empty()) {
        time_new_to_old.push_back(t);
    } else {
        if (t != time_new_to_old.back()) {
            time_new_to_old.push_back(t);
            edges_idx.emplace_back(edges.size() - 1);
        }
    }

    int format_t = int(time_new_to_old.size() - 1);
    tnu[u].emplace_back(make_pair(v, format_t));
    tnv[v].emplace_back(make_pair(u, format_t));
}

// not change max_degree
void BiGraph::deleteEdge(vid_t u, vid_t v)
{
    for (int i = 0; i < degree_v1[u]; ++i)
    {
        int vv = neighbor_v1[u][i];
        if (vv == v)
        {
            swap(neighbor_v1[u][i], neighbor_v1[u][degree_v1[u] - 1]);
            --degree_v1[u];
            neighbor_v1[u].pop_back();
            num_edges--;//only once!!!
            break;
        }
    }

    if (degree_v1[u] + 1 == v1_max_degree) {
        v1_max_degree = 0;
        for (auto d : degree_v1) {
            v1_max_degree = v1_max_degree < d ? d : v1_max_degree;
        }
    }

    for (int i = 0; i < degree_v2[v]; ++i)
    {
        int uu = neighbor_v2[v][i];
        if (uu == u)
        {
            swap(neighbor_v2[v][i], neighbor_v2[v][degree_v2[v] - 1]);
            --degree_v2[v];
            neighbor_v2[v].pop_back();
            break;
        }
    }

    if (degree_v2[v] + 1 == v2_max_degree) {
        v2_max_degree = 0;
        for (auto d : degree_v2) {
            v2_max_degree = v2_max_degree < d ? d : v2_max_degree;
        }
    }
}

bool BiGraph::isEdge(vid_t u, vid_t v)
{
//    for (unsigned int & it : neighbor_v1[u]) {
//        if (it == v) return true;
//    }
//    return false;

    return any_of(neighbor_v1[u].begin(), neighbor_v1[u].end(), [v](unsigned int & it){
        return it == v;
    });
}

int BiGraph::get_left_index_with_fixed_left_k(vid_t u, int left_k) {
    if (left_index[u].size() > left_k) return left_index[u][left_k];
    else return 0;
}

bool isNodesetEqual(unordered_set<vid_t>& set1, unordered_set<vid_t>& set2) {
    if (set1.size() != set2.size()) {
        return false;
    }
    for (unsigned int it : set1) {
        auto got = set2.find(it);
        if (got == set2.end()) {
            return false;
        }
    }
    return true;
}


void inv(BiGraph& g) {
    swap(g.degree_v1, g.degree_v2);
    swap(g.left_index, g.right_index);
    swap(g.num_v1, g.num_v2);
    swap(g.neighbor_v1, g.neighbor_v2);
    swap(g.neighborHash_v1, g.neighborHash_v2);
    swap(g.v1_max_degree, g.v2_max_degree);
    swap(g.left_delete, g.right_delete);
}