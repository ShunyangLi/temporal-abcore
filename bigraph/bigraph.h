#pragma once
#pragma once
#ifndef __BIGRAPH_H
#define __BIGRAPH_H


#include "../utility/utility.h"
#include <set>


// use vertex_block to record the vertices and ending time
struct vertex_block {
    int te;
    std::vector<vid_t> nodeset;
    // to record the max node that combine the componet
    vertex_block* parent = nullptr;
    vertex_block* child = nullptr;
};

using namespace std;
class BiGraph
{

public:

    explicit BiGraph(const std::string& dir);
    BiGraph();
    ~BiGraph() = default;

    void addEdge(vid_t u, vid_t v);
    void addEdgeT(vid_t u, vid_t v, vid_t t);
    void deleteEdge(vid_t u, vid_t v);
    bool isEdge(vid_t u, vid_t v);
    [[nodiscard]] num_t getV1Num() const { return num_v1; }
    [[nodiscard]] num_t getV2Num() const { return num_v2; }
    num_t getV1Degree(vid_t u) { return degree_v1[u]; }
    num_t getV2Degree(vid_t u) { return degree_v2[u]; }
    std::vector<vid_t> & getV2Neighbors(vid_t u) { return neighbor_v2[u]; }
    std::vector<vid_t> & getV1Neighbors(vid_t u) { return neighbor_v1[u]; }

public:

    void init(unsigned int num_v1, unsigned int num_v2);
    void loadGraph(const std::string& path);

    std::string dir;
    num_t num_v1;
    num_t num_v2;
    num_t num_edges;

    std::vector<std::vector<vid_t>> neighbor_v1;
    std::vector<std::vector<vid_t>> neighbor_v2;

    // this part for timestamp
    std::vector<vid_t> time_new_to_old;
    // edge index
    std::vector<vid_t> edges_idx;
    std::vector<std::pair<vid_t,vid_t>> edges;
    int tmax = 0;

    std::vector<std::vector<std::pair<vid_t,vid_t>>> tnu;
    std::vector<std::vector<std::pair<vid_t,vid_t>>> tnv;

    // u, alpha, beta, timestamp
    vector<vector<vector<vector<pair<int,int>>>>> u_index;
    std::vector<std::vector<std::vector<std::vector<std::pair<int,int>>>>> v_index;

    // record abcore subgraph, alpha, beta, ts -> records
    vector<vector<vector<vector<vertex_block*>>>> tbcore_uindex;
    vector<vector<vector<vector<vertex_block*>>>> tbcore_vindex;

    vector<unordered_map<int, int>> ucn;
    vector<unordered_map<int, int>> vcn;


    std::vector<std::set<vid_t>> neighborHash_v1;
    std::vector<std::set<vid_t>> neighborHash_v2;

    std::vector<int> degree_v1;
    std::vector<int> degree_v2;

    std::vector<num_t> core_v1;
    std::vector<num_t> core_v2;

public:

    //KKCore index left (x,*) right (*,x)
    std::vector<std::vector<int>> left_index;
    std::vector<std::vector<int>> right_index;
    int v1_max_degree;
    int v2_max_degree;
    std::vector<bool> left_delete;
    std::vector<bool> right_delete;
    // for dynamic update
    std::vector<std::vector<int>> left_index_old;
    std::vector<std::vector<int>> right_index_old;
    //BiGraph operator=(const BiGraph& g);
    int delta = -1;

public:
    int get_left_index_with_fixed_left_k(vid_t u, int left_k);
    //BiGraph& operator=(const BiGraph& g_);
};

extern void inv(BiGraph& g);

#endif  /* __BIGRAPH_H */