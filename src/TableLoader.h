#ifndef TABLELOADER_H
#define TABLELOADER_H

#include <cstdio>
#include <vector>
#include "DataStructures.h"

PROBLEM readReactionTable(const char* filename);
void readLikelihoodTable(PROBLEM &ProblemSpace, const char* filename);
int rxnByName(const RXNSPACE &rxnspace, const char* name);

#endif
