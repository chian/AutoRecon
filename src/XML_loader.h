#ifndef _XML_LOADER
#define _XML_LOADER

/* C++ library */
#include "DataStructures.h"

using std::vector;
using std::map;

/* From XML_loader.c - parseALL is the wrapper for the other smaller functions */
void parseALL(char* file1, char* file2, PROBLEM &ProblemSpace);
void parseData(char *docname, vector<GROWTH> &growth);
void setUpMaintenanceReactions(PROBLEM &ProblemSpace);
void Load_Stoic_Part(vector<REACTION> &reaction, const vector<METABOLITE> &metabolite);
void addMetNameToStoichs(vector<REACTION> &reaction, METSPACE &metspace);
void addMetNameToStoichs(vector<GROWTH> &growth, METSPACE &metspace);
void loadInputsOutputs(vector<GROWTH> &growth, int growthIdx, METSPACE &metspace, vector<int> &inputIds, vector<int> &outputIds);
void makeMagicBridges(PROBLEM &ProblemSpace);
void identifyFreeReactions(vector<REACTION> &reactions);
bool secondaryPairPresent(const vector<int> &secondaryList1, const vector<int> &secondaryList2, const REACTION &reaction, int metId, map<int, int> &metId2StoichIdx);
void checkConsistency(const PROBLEM &ProblemSpace);
bool pairDone(const vector<vector<int> > &pairList, const vector<int> &onePair);

#endif
