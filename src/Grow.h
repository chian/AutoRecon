#ifndef GROW_H
#define GROW_H

/******************** GROW.h  ****************
    Contains functions for the inner loop
    (gap finding, gap filling, and choosing of
     gap fill solutions based on pre-defined
     scoring criteria)
****************************************/

#include "shortestPath.h"
#include "genericLinprog.h"
#include "pathUtils.h"
#include "visual01.h"
#include "MersenneTwister.h"
#include "DataStructures.h"
#include "RunK.h"

#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Makes a simulatable model out of a vector of PATHSUMMARY, including magic exits for everything
   and entrances for secondary lone and pair metabolites. */
PROBLEM setUpGapfill(const PROBLEM &problemSpace, const vector<PATHSUMMARY> &pList);

/* Primary function for calling gapfind/gapfill and taking care of other stuff 
   like identifying essential magic entrances/exits,
   suggesting gapfill solutions for each other exit as deemed necessary,
   and connecting ETC reactions. */
ANSWER gapfillWrapper(const PROBLEM &problemSpace, const vector<PATHSUMMARY> &pList, const GROWTH &growth);

/* Note - workingRxns and workingMets are purposely passed by value for now...if I can get the add/subtracting exactly right I may be able to avoid it  */
vector<GAPFILLRESULT> gapFindGapFill(PROBLEM &model, const PROBLEM &problemSpace, int gapfillK);
vector<vector<int> > fillGapWithDijkstras(RXNSPACE &workingRxns, METSPACE &workingMets, 
					  PROBLEM &wholeProblem, int toFix, 
					  int direction, int gapfillK);
void addGapfillResultToProblem(PROBLEM &model, const PROBLEM &problemSpace, 
			       const GAPFILLRESULT &gapfillResult, int whichK);
vector<int> minimizeExits(PROBLEM &model);

/* [inner loop] Genetic algorithm helper functions */
vector<int> findSolutionsMinimizingExits(const PROBLEM &baseModel, const PROBLEM &problemSpace, const vector<GAPFILLRESULT> &res, vector<int> &essential);
double innerScore(PROBLEM &model, const PROBLEM &problemSpace, vector<int> &essentialExits);
vector<int> randomK(const vector<GAPFILLRESULT> &res, MTRand &rng);
int getRandomK(const GAPFILLRESULT &res, MTRand &rng);
INNERPOPSTORE crossOver(const vector<INNERPOPSTORE> &parents, int mostFitCutoff, MTRand &rng);
INNERPOPSTORE mutate(const vector<INNERPOPSTORE> &parents, const vector<GAPFILLRESULT> &res, int nToMutate, int mostFitCutoff, MTRand &rng);

/* Switch from one GROWTH condition to another */
void setSpecificGrowthConditions(PROBLEM &model, const GROWTH &growth);

/* Test if knockout of a particular gene is predicted to be lethal in the given PROBLEM... */
bool knockoutLethality(PROBLEM &modified, const map<string, vector<VALUESTORE> > &geneRxnMap, const string &genename);

/* Modify the GAM or NGAM values */
void addATPM(PROBLEM &A, ANSWER &B);

#endif
