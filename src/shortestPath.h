#ifndef shortestPath_h
#define shortestPath_h

#include "DataStructures.h"
#include "pathUtils.h"

#include <queue>
#include <map>
#include <set>
#include <vector>

using std::map;
using std::queue;
using std::set;
using std::vector;

PATH findShortestPath(const RXNSPACE &rxnspace, const METSPACE &metspace, const METSPACE &inputs, const METABOLITE &output, set<BADIDSTORE> &badIds);
void tracePath(const METSPACE &metspace, const RXNSPACE &rxnspace, map<int, int>  &precursorRxnIds, map<int, bool> &metsExplored, vector<int> &rxnIds, vector<int> &inputIds,
               vector<int> &rxnDirections, queue<int> &nodeList);
vector<int> opposite(const REACTION* rxn, int metId, int sgn);

#endif
