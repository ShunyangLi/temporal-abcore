#pragma once
#ifndef ABCORE_H
#define ABCORE_H
#include "bigraph.h"

using namespace std;

int coreIndexKCore(BiGraph& g);

double update_bicore_index(BiGraph& g, vid_t u, vid_t v, bool addition,
                           vector<vid_t>& au, vector<vid_t>& av);

#endif