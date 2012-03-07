#ifndef ETC_H
#define ETC_H

#include "DataStructures.h"
#include "shortestPath.h"
#include "pathUtils.h"
#include <vector>

/* Main functions */
void ETC(const PROBLEM &ProblemSpace, vector<NETREACTION> &bigout, double cutoff);
void identifyETCrxns(const PROBLEM &ProblemSpace, vector<int> &ETClist);
void identifyETCrxns(const PROBLEM &ProblemSpace, vector<int> &ETClist, double cutoff);
vector<NETREACTION> ETC_dir_check(const RXNSPACE &fullrxnspace, const vector<NETREACTION> &bigout);

/* Utilities */
int Name2Ids(const vector<METABOLITE> &metabolite, const char *met_name);
void decompose(const vector<int> &rxnDirIds, vector<int> &rxnIds,  vector<int> &rxnDirs);
void flip(vector<int> &rxnDirs);
NETREACTION flip(const NETREACTION &one);
void convert(NETREACTION &net, REACTION add, int rxnDirId);
vector<vector<int> > PossiblePairs(const vector<STOICH> &fromnet,  const METSPACE &cofactors);
int netrxn(NETREACTION &net, const PROBLEM &ProblemSpace, int rxnDirId);

#endif
