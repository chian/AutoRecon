#ifndef PATHS2MODEL_H
#define PATHS2MODEL_H

#include "shortestPath.h"
#include "Grow.h"
#include "genericLinprog.h"
#include "pathUtils.h"
#include "visual01.h"
#include "DataStructures.h"
#include "RunK.h"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>

/* Unsynonymize */
void addR(PATHSUMMARY &psum, PROBLEM &ProblemSpace);
void addR(vector<vector<vector<PATHSUMMARY> > > &psum, PROBLEM &ProblemSpace);

vector<int> findRxnsToIncludeFromSyn(const REACTION &TMPSYN, const RXNSPACE &rxnspace, int dir);
PATHSUMMARY replaceWithRealRxnIds(const PATHSUMMARY &psum, const vector<PATHSUMMARY> &pList, PROBLEM &ProblemSpace);
vector<int> fillMagicBridges(const PATHSUMMARY &psum, const vector<PATHSUMMARY> &pList, const PROBLEM &ProblemSpace, int id);

/* Random utility functions */
void checkExchangesAndTransporters(PROBLEM &working, const PROBLEM &ProblemSpace, const vector<int> &metIds, const vector<int> &dirs);
void makeSimulatableModel(const vector<PATHSUMMARY> &pList, const PROBLEM &ProblemSpace, const REACTION &biomass, const vector<int> &magicIds, const vector<int> &magicDirs,
                          PROBLEM &working, int &baseNum);


#endif
