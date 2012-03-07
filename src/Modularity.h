#ifndef MODULARITY_H_
#define MODULARITY_H_

#include<vector>
#include"pathUtils.h"
#include"DataStructures.h"

void AllHardIncludes(const vector<REACTION> &biglist, vector<REACTION> &smalllist);
void SynIncludes(PROBLEM &Model);

#endif
