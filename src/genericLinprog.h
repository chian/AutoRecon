#include <glpk.h>
#include <vector>
#include "DataStructures.h"
#include "RunK.h"

#ifndef GENERICLINPROG_H
#define GENERICLINPROG_H

/* This is for reverse-compatibility */
vector<double> FBA_SOLVE(const vector<REACTION> &rxns, const METSPACE &metspace);
vector<double> FBA_SOLVE(const RXNSPACE &rxnspace, const METSPACE &metspace);
void FVA_SOLVE(const RXNSPACE &rxnspace, const METSPACE &metspace, double optPct, vector<double> &minflux, 
	       vector<double> &maxflux);

/* The best way to ensure consistency is perhaps to force the user to make a new one of these if the problem changes - 
   I try to enforce that with private and public here */
class GLPKDATA {
 public:
  GLPKDATA();
  GLPKDATA(const RXNSPACE &rxns, const METSPACE &mets, const vector<int> &objId, const vector<double> &objCoeff, int sense);
  ~GLPKDATA();
  /* TODO: need copy constructor */

  /* Solver routines */
  int gapFindLinprog(vector<int> &usedExits);
  vector<double> FBA_SOLVE();
  void FVA_SOLVE(vector<double> &minFlux, vector<double> &maxFlux, double optPercentage);
  int FastFVA(vector<double> &minFlux, vector<double> &maxFlux, double optPercentage);
  void printPrivateStuff();

 private:
  int numrows;
  int numcols;
  double* lb;
  double* ub;
  int totalDataSize;
  int* ia;
  int* ja; 
  double* ar;

  /* This is only temporary - I need it to help me debug this damned thing */
  RXNSPACE rxnsUsed;
  METSPACE metsUsed;

  vector<int> objIdx;
  vector<double> objCoef;

  int objSense; /* -1 = MIN, 1 = MAX */

  glp_prob* problem;

  /* I intentionally did not make a public version of this function - it is much less confusing to only have ONE place to change LB and UB
     and that is in the RXNSPACE itself. Just make a new GLPKDATA if you need to! */
  void changeLb(double newLb, int idx);
  void changeUb(double newUb, int idx);
  void changeObjective(vector<int> newIdx, vector<double> newCoeffs);
  void validateSense(int sense);
  void initialize(const METSPACE &metspace, const RXNSPACE &rxnspace, vector<int> objId, vector<double> objCoeff, int sense);
  void setUpProblem();

  void printGlpkError(int errorCode);

  /* This has to be declared as static because we're making a function pointer to it and we can't have the address changing on us mid-program, can we? */
  static int suppressGLPKOutput(void *info, const char *s);
};

#endif
