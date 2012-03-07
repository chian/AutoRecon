#ifndef RUNK_H
#define RUNK_H

#include "DataStructures.h"
#include "shortestPath.h"
#include "genericLinprog.h"
#include "pathUtils.h"
#include "visual01.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>

/* Class for results from RandomMagic.cpp - needed as a container
   to shuffle around results */
class MAGICRESULT{
 public:
  vector<int> metIds;
  vector<int> dir; /* 1 means exit, -1 means entrance, -2 means entrance for a cofactor */
  int numActiveRxns;
  int numTotalRxns;
  int numExits; /* For purposes of scoring, cofactor pairs are both considered an "exit" */
  int numEntrances;
  int numCofEntrances;
  int outputId; /* OutputID for the current path */
  
  MAGICRESULT() {
    numActiveRxns = -1;
    numTotalRxns = -1;
    outputId = -1;
    numExits = -1;
    numEntrances = -1;
    numCofEntrances = -1;
  }

  void reset() {
    metIds.clear();
    dir.clear();
    numActiveRxns = -1;
    numTotalRxns = -1;
    numExits = -1;
    numEntrances = -1;
    numCofEntrances = -1;
    outputId = -1;
  }
};

void InputSetup(int argc, char *argv[], PROBLEM &ProblemSpace);
int inOutPair(int metId, const METSPACE &metspace);
REACTION MagicExit(const vector<REACTION> &reaction, int met_id, const char* name);
REACTION MagicExit(const vector<REACTION> &reaction, int met_id, const char* name, int R);
void Run_K(PROBLEM &ProblemSpace, vector<MEDIA> &media, RXNSPACE &rxnspace, int outputId, int K, 
	   vector<PATHSUMMARY> &result, int direction, int growthIdx);
void Run_K2(PROBLEM &ProblemSpace, vector<MEDIA> &media, int outputId, int K, int startingK,
	   vector<PATHSUMMARY> &result, int direction, int growthIdx);
void FirstKPass(PROBLEM &ProblemSpace, int K, vector<vector<vector<PATHSUMMARY> > > &psum);
void SecondKPass(PROBLEM &ProblemSpace, int K, vector<vector<vector<PATHSUMMARY> > > &psum);

void GoodReversible(PROBLEM &ProblemSpace);
void ReverseReversible(vector<REACTION> &reaction);

void AddNeededReversibles(vector<REACTION> &reaction, const PATH &path, map<int,int> rxnIds2Idx);

void Load_Forwards_Reactions(vector<REACTION> &reaction, const vector<int> &known_rxn_ids);
vector<int> Load_Outputs_From_Growth(const PROBLEM &parentSpace, int growthIdx);
vector<int> Load_Inputs_From_Growth(const GROWTH &growth);

REACTION MakeObjRxn(const vector<STOICH> &stoich);
REACTION MakeBiomassFromGrowth(const GROWTH &growth);

REACTION GrowthExit(const vector<REACTION> &reaction, int met_id, int reversible, 
		    double fluxBound, const char* name);
void adjustLikelihoods(vector<REACTION> &rxnList, double spont_likely, double black_magic_likely, double hard_include_likely, double no_likely, bool adjustNonspecial);
void calcMetRxnRelations(const RXNSPACE &rxnspace, METSPACE &mets);
void calcMetRxnRelations_nosec(const RXNSPACE &rxns, METSPACE &mets);

#endif
