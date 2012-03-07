#ifndef _PATHUTILS
#define _PATHUTILS

#include "DataStructures.h"

bool isExternalMet(const char *name, const char* E_tag);

/* Vector utilities */
bool anyVec(const vector<bool> &vec);
void setdiff(vector<int> &A, vector<int> &B, vector<int> &diff);
void setdiff1(vector<int> &A, const vector<int> &B);
vector<int> setdiff2(vector<int> A, vector<int> B);
vector<int> custom_intersect(vector<int> A, vector<int> B);
void catVectors1(vector<int> &vector1, const vector<int> &vector2);
template <class T> void custom_unique(vector<T> &inVector);
template <class T> vector<T> catVectors(vector<T> vector1, vector<T> vector2);

/* Get products and reactants lists on the OPPOSITE side of a reaction from a given metabolite */
void getProducts(const REACTION &rxn, int reactantId, vector<int> &res);
void pullOutRxnsbyIds(const RXNSPACE &rxnspace, const vector<int> &rxnsToSearchIds, RXNSPACE &result);
void pullOutRxnsbyIds(const PROBLEM &ProblemSpace, const vector<int> &rxnsToSearchIds, RXNSPACE &result);
void pullOutMets(const RXNSPACE &rxnspace, const METSPACE &metspace, METSPACE &metabolite);

/* Syn list making */
int instoich(int id, const vector<STOICH> &list);
int diff2rxns(const REACTION &one, const REACTION &two);
void MakeSynList(PROBLEM &ProblemSpace);

/* Get a reaction coefficient for a reaction "rxnId" and metabolite "metId" */
double rxnCoeff(const RXNSPACE &rxnspace, int rxnId, int metId);

/* Utility functions */
int Name2Ids(const vector<METABOLITE> &metabolite, const char *met_name);
int Name2Ids(const vector<REACTION> &reaction, const char *rxn_name);

vector<PATHSUMMARY> flattenPsum(const vector<vector<vector<PATHSUMMARY> > > &from);
vector<PATHSUMMARY> uniquePsum(const vector<PATHSUMMARY> &from);

int rougheq(const double one, const double two);
int rougheq(const double one, const double two, const double constant);

vector<int> getAllPathRxns(const vector<PATHSUMMARY> &pList);
vector<int> getAllPathMets(const vector<PATHSUMMARY> &pList, const PROBLEM &problemSpace, int whichOne);
vector<int> getAllSecondaryIds(const PROBLEM &ProblemSpace);

void turnOffMagicExits(RXNSPACE &rxnspace, const vector<int> &exitsToKeep);

/* Magic bridge-related functions */
bool fixedBridge(const REACTION &rxn, const REACTION &bridge);
bool isDanglingSecondary(const REACTION &reaction, const METSPACE &metspace, int met_id);
void biomassMB(PROBLEM &problemSpace);
void changeMetToMB(REACTION &rxn, int met_id, const char* metName, bool isDangling);
void modifyRxnsForMagicBridge(PROBLEM &problemSpace, int bridgeId, vector<int> &rxnIds, vector<bool> &toErase);

void modifyGAM(ANSWER &model, double amountToChange);
void modifyNGAM(ANSWER &model, double newAtpm);

/* For template classes to work for anything but built-in types (int, double, etc), 
   the definition of the function has to be in the same file as the prototype.
   So I put them here... */
template <class T>
void custom_unique(std::vector<T> &inVector) {
  typename std::vector<T>::iterator iter;
  sort(inVector.begin(), inVector.end());
  iter = unique(inVector.begin(), inVector.end());
  inVector.resize(iter-inVector.begin());
}

template <class T>
vector<T> catVectors(vector<T> vector1, vector<T> vector2) {
  vector<T> newVector = vector1;
  newVector.reserve(vector1.size() + vector2.size());
  while(!vector2.empty()) {
    newVector.push_back(vector2.back());
    vector2.pop_back();
  }
  return newVector;
}

#endif
