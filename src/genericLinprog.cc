#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <glpk.h>
#include <vector>

#include "DataStructures.h"
#include "Exchanges.h"
#include "genericLinprog.h"
#include "MyConstants.h"
#include "pathUtils.h"
#include "Printers.h"
#include "Exchanges.h"
#include "RunK.h"
#include "Grow.h"

/* Useful plug-in functions */
/* For backwards compatibility - assumes that you want db.BIOMASS as the objective
 with a coefficient of 1 and a sense of 1 (MAX) */
vector<double> FBA_SOLVE(const RXNSPACE &rxnspace, const METSPACE &metspace) {
  vector<int> obj(1, _db.BIOMASS);
  vector<double> coeff(1, 1.0f);
  GLPKDATA prob(rxnspace, metspace, obj, coeff, 1);
  vector<double> result = prob.FBA_SOLVE();

  return result;
}

vector<double> FBA_SOLVE(const vector<REACTION> &rxnList, const METSPACE &metspace) {
  RXNSPACE rxnspace(rxnList);
  vector<double> result = FBA_SOLVE(rxnspace, metspace);
  return result;
}

/* Assumes you want db.BIOMASS as the objective to compare against when running FVA, and a coefficient of 1 for that objective */
void FVA_SOLVE(const RXNSPACE &rxnspace, const METSPACE &metspace, double optPct, vector<double> &minflux, 
	       vector<double> &maxflux) {
  minflux.clear(); maxflux.clear();
  vector<int> obj(1, _db.BIOMASS);
  vector<double> coeff(1, 1.0f);
  GLPKDATA prob(rxnspace, metspace, obj, coeff, 1);
  prob.FVA_SOLVE(minflux, maxflux, optPct);
  return;
}

/**************** Public Methods ****************/

GLPKDATA::~GLPKDATA() {
  free(lb);  free(ub);  free(ar);  free(ia); free(ja);
  glp_delete_prob(problem);
}

GLPKDATA::GLPKDATA() {
  printf("ERROR: Not allowed to use the default constructor for GLPKDATA. Must use one of the other constructors provided below!\n");
  assert(false);
}

/* Objectives are all assumed to be zero except for IDs listed in objId */
GLPKDATA::GLPKDATA(const RXNSPACE &rxnspace, const METSPACE &metspace, const vector<int> &objId, const vector<double> &objCoeff, int sense) {
  initialize(metspace, rxnspace, objId, objCoeff, sense);
}

/* Solve the simplex problem vanilla */
vector<double> GLPKDATA::FBA_SOLVE() {
  setUpProblem();

  /* Turn pre-solving on... 
     without doing this I obtain a basic solution that is not feasible...that isn't useful */
  glp_smcp param;
  glp_init_smcp(&param);

  /* We only want optimial solutions */
  param.presolve=GLP_ON;

  /* Use dual and then switch to primal if dual fails 
   I had more success with the numerical stability of this method
   than with just primal simplex (since the dual was having issues being feasible) */
  param.meth = GLP_DUALP;

  glp_std_basis(problem);
  //  glp_adv_basis(problem, 0);

  /*Scale problem with equilibrium scaling to improve numerical stability 
    (this is the default behavior in GLPKMEX) */
  glp_scale_prob(problem, GLP_SF_EQ);

  int res = glp_simplex(problem, &param);

  /* Check solution quality */
  if(_db.DEBUGFBA) {
    LPXKKT lpxStruct;
    lpx_check_kkt(problem, 1, &lpxStruct);
    printf("QUALITY REPORT: Primal solution quality: %c\n", lpxStruct.pe_quality);
    printf("Primal solution feasibility: %c\n", lpxStruct.pb_quality);
    printf("Dual solution quality: %c\n", lpxStruct.de_quality);
    printf("Dual solution feasibility: %c\n", lpxStruct.de_quality);
  }

  vector<double> fluxResult(numcols, 0);
  if(res!=0) { printf("ERROR: Numerical issues with solution... such cases are treated the same as NO GROWTH cases \n");
    printGlpkError(res);
    glp_write_lp(problem, NULL, "PROBLEMATIC_PROBLEM");
    return fluxResult;
  }

  fluxResult.clear();
  /* Return solution (GLPK doesn't rearrange columns) */
  for(int i=0; i<numcols; i++) {
     fluxResult.push_back(glp_get_col_prim(problem, i+1));
  }
  if(_db.DEBUGFBA) {
    printDoubleVector_rxns(rxnsUsed, fluxResult);
  }
  return fluxResult;
}

void GLPKDATA::FVA_SOLVE(vector<double> &minFlux, vector<double> &maxFlux, double optPercentage) {

  /* You can't have a negative percent or get an objective value more than 100% of the maximum */
  assert(optPercentage >= 0.0f & optPercentage <= 100.0f);

  setUpProblem();
  int status = FastFVA(minFlux, maxFlux, optPercentage);
  assert(status == 0);

  for(int i=0; i<rxnsUsed.rxns.size(); i++) { 
    if(_db.DEBUGFVA) { printf("FVA rxn %s min = %4.5f max = %4.3f lb = %4.3f ub = %4.3f\n", rxnsUsed.rxns[i].name, minFlux[i], maxFlux[i], rxnsUsed.rxns[i].lb, rxnsUsed.rxns[i].ub); }
  }

  return;
}

/* Returns a list of reaction IDs for unnecessary magic exits *
 Requires the metspace and rxnspace in GLPKDATA 

 NOTE: This irreversibly changes the GLPKDATA structure, so I need to write a fucntion to re-initialize it afterwards... 

 To get this, solves the problem:

 MAXIMIZE c1*(NMEF) - c1*(NMER) - c2*(MEF) + c2*(MER) 

 where: c1 is the "normal bonus" parameter
        c2 is the "exit penalty" parameter
        NMEF = Non-magic exits with always-positive flux 
        NMER = Non-magic exits with always-negative flux
        MEF = Magic exits with always-positive flux
        MER = Magic exits with always-negative flux

 This formulation is not perfect but it at least gets us close...and guarantees that we obtain growth from the resulting pruning */
int GLPKDATA::gapFindLinprog(vector<int> &usedExits) {

  usedExits.clear();

  if(objIdx.size() > 1) { printf("ERROR: gapFindLinprog only supports single-reaction objectives\n"); assert(false); }
  assert(objSense == 1);

  /* Minimum possible value of objective function (forces solutions that grow) */
  double objMin = 1.0f;

  double normal_bonus = 100.0f;
  double exit_penalty = 100.0f;

  /* Needed for re-initializing to the previous state from before */
  vector<int> origIds; 
  for(int i=0; i<objIdx.size(); i++) { origIds.push_back(rxnsUsed.rxns[objIdx[i]-1].id); }
  vector<double> origCoeff = objCoef;

  /* Get minimum and maximum possible fluxes - need this to ensure that the reaction is going the right direction */
  vector<double> minflux; vector<double> maxflux;
  FVA_SOLVE(minflux, maxflux, 0.0f);

  /* Change the LB and the UB for the original objective(s) before we change them (this is the part I'm not sure how to do if we allow multi-rxn objectives) */
  for(int i=0; i<objIdx.size(); i++) {
    lb[objIdx[i]] = objMin; ub[objIdx[i]] = 1000.0f;
  }

  /* Identify coefficients that are going the opposite way from how they should, and switch them to go the other way 
   Also, find reactions that only can go one way (either forward or reverse) and add them to the coefficient list */
  objIdx.clear(); objCoef.clear();
  for(int i=0; i<rxnsUsed.rxns.size(); i++) {
    if(minflux[i] < -1E-5 && maxflux[i] < 1E-5) { /* These get a negative coefficient because they will have a negative flux [FIXME - also need to exclude exchanges?] */
      if(rxnsUsed.rxns[i].id < _db.MINFACTORSPACING) {
	objIdx.push_back(i+1); objCoef.push_back(-normal_bonus);
      } else if(rxnsUsed.rxns[i].id >= _db.BLACKMAGICFACTOR && rxnsUsed.rxns[i].id < _db.BLACKMAGICFACTOR + _db.MINFACTORSPACING) { 
	objIdx.push_back(i+1); objCoef.push_back(exit_penalty);
      }
    } else if(minflux[i] > -1E-5 & maxflux[i] > 1E-5) { /* These get positive coefficients because they will have a positive flux */
      if(rxnsUsed.rxns[i].id < _db.MINFACTORSPACING) {
	objIdx.push_back(i+1); objCoef.push_back(normal_bonus);
      }
      else if(rxnsUsed.rxns[i].id >= _db.BLACKMAGICFACTOR && rxnsUsed.rxns[i].id < _db.BLACKMAGICFACTOR + _db.MINFACTORSPACING) { 
	objIdx.push_back(i+1); objCoef.push_back(-exit_penalty);
      }
    }
  }

  vector<double> result = this->FBA_SOLVE();
  if( result[rxnsUsed.idxFromId(_db.BIOMASS)] < _db.GROWTH_CUTOFF ) { printf("ERROR: Gapfind failed to return a solution...\n"); return -1; }

  /* Re-initialize */
  this->initialize(metsUsed, rxnsUsed, origIds, origCoeff, objSense);

  /* Add needed magic exits to the list of exits that are used */
  for(int i=0; i<result.size(); i++) {
    if( rougheq(result[i], 0.0f, _db.FLUX_CUTOFF)) { continue; }
    if( _db.DEBUGGAPFILL ) { printf("Reaction %s - gapFindLinprog result = %4.3f\n", rxnsUsed.rxns[i].name, result[i]); }
    if(rxnsUsed.rxns[i].id < _db.BLACKMAGICFACTOR || rxnsUsed.rxns[i].id > _db.BLACKMAGICFACTOR + _db.MINFACTORSPACING) { continue; }
    /* If the min flux is close to 0 ignore it - it is not needed (1E-5 is OK as long as we do the conditioning step above making all the coeffs = 1) */
    usedExits.push_back(rxnsUsed.rxns[i].id);
  }

  custom_unique(usedExits);

  return 1;
}

/************************ Private Methods ************************/

/* This code is modified from FastFVA by S. Gudmundsson and I. Theile - see  S. Gudmundsson and I. Thiele. Computationally efficient flux variability analysis, BMC Bioinformatics, 2010, 11:489
   (http://notendur.hi.is/ithiele/software/fastfva.html) 
   It is licensed under the LGPL. 

   The problem should be initialized / setup BEFORE calling this function (DO NOT call this function directly - call it through FVA_SOLVE)  */

#define FVA_SUCCESS        0
#define FVA_INIT_FAIL      1
#define FVA_MODIFIED_FAIL  2
#define TIME_RESTART_LIM   60

int GLPKDATA::FastFVA(vector<double> &minFlux, vector<double> &maxFlux, double optPercentage) {

  minFlux.clear(); minFlux.assign(numcols, 0.0f);
  maxFlux.clear(); maxFlux.assign(numcols, 0.0f);

  // Parameters for the glpk optimizer, use mostly default settings
  glp_smcp param;
  glp_init_smcp(&param);
  param.presolve = GLP_ON;
  
  // Objective (min or max)
  glp_set_obj_dir(this->problem, (this->objSense==-1) ? GLP_MIN : GLP_MAX);
  
  const double tol = _db.FLUX_CUTOFF;
  
  // Solve the initial problem
  glp_adv_basis(this->problem, 0);
  //  glp_scale_prob(this->problem, GLP_SF_EQ);

  int ret = glp_simplex(this->problem, &param);
  if (ret != 0)  {  
    /* Attempt to recover... */
    glp_std_basis(this->problem);
    int ret = glp_simplex(this->problem, &param);
    if(ret != 0) {    printGlpkError(ret); return FVA_INIT_FAIL; }
  }

  double z = glp_get_obj_val(this->problem);
  
  // Determine the value of objective function bound
  double TargetValue = 0;
  if (glp_get_obj_dir(this->problem) == GLP_MAX) {  TargetValue = floor(z/tol)*tol*optPercentage/100.0;   }
  else  { TargetValue = ceil(z/tol)*tol*optPercentage/100.0;  }

  // Add a constraint which bounds the objective, c'v >= objValue in case of max, (<= for min)
  int m = glp_add_rows(this->problem, 1);
  int n = glp_get_num_cols(this->problem);

  if (glp_get_obj_dir(this->problem) == (int)GLP_MAX) { glp_set_row_bnds(this->problem, m, GLP_LO, TargetValue, 0.0);  }
  else  { glp_set_row_bnds(this->problem, m, GLP_UP, 0, TargetValue);  }
  
  int* ind = (int*)malloc( (n+1)*sizeof(int));
  double* val = (double*)malloc( (n+1)*sizeof(double));
  for (int j = 1; j <= n; j++) {
    ind[j]=j;
    val[j]=glp_get_obj_coef(this->problem, j);
  }
  glp_set_mat_row(this->problem, m, n, ind, val);
  free(ind);
  free(val);
  
  // Zero all objective function coefficients (starting point for FVA)
  for (int j = 1; j <= n; j++)   {  glp_set_obj_coef(this->problem, j, 0.0f); }
  
  // Solve all the minimization problems first. The difference in the optimal
  // solution for two minimization problems is probably much smaller on average
  // than the difference between one min and one max solutions, leading to fewer
  // simplex iterations in each step.
  param.presolve = GLP_OFF;
  param.msg_lev = GLP_MSG_OFF;
  param.tm_lim = 1000*TIME_RESTART_LIM;
  
  for (int iRound = 0; iRound < 2; iRound++)  {
    glp_set_obj_dir(this->problem, (iRound==0) ? GLP_MIN : GLP_MAX);
    for (int k = 0; k < this->numcols; k++) {
      glp_set_obj_coef(this->problem, k+1, 1.0f);
      ret = glp_simplex(this->problem, &param);
      if (ret != 0) {
	// Numerical difficulties or timeout
	printf("Numerical problems...\n");
	printf("K: %d of %d\n", k, numcols);
	param.tm_lim = INT_MAX;
	param.presolve = GLP_ON;
	glp_adv_basis(this->problem, 0);
	ret = glp_simplex(this->problem, &param);
	if (ret != 0) {
	  return FVA_MODIFIED_FAIL;
	}
	param.presolve = GLP_OFF;
	param.tm_lim = 1000*TIME_RESTART_LIM;
      }
      /* I think this was a bug to put it above...  */
      glp_set_obj_coef(this->problem, k+1, 0.0f);      
      if (glp_get_obj_dir(this->problem) == GLP_MIN)  {  minFlux[k]= glp_get_obj_val(this->problem); }
      else{ maxFlux[k]=glp_get_obj_val(this->problem);  }
    }
  }

  /* Reset the problem, because we don't want the added row to mess things up for us... */
  setUpProblem();
  
  return FVA_SUCCESS;
}


/* Set up the glp_prob "problem" to have all the data it needs to run a simulation based on the current arrays 
 Also has the possibility of deleting the old problem and starting fresh IF you don't modify any of those private variables! 

Uses the following class variables: problem (resets before using)
 numrows, numcols
 lb, ub
 ia, ja, ar 
 totalDataSize 
 objIdx, objCoef */
void GLPKDATA:: setUpProblem() {
  glp_delete_prob(problem);
  problem = glp_create_prob();

  if(objSense == -1) { glp_set_obj_dir(problem, GLP_MIN); } 
  else { glp_set_obj_dir(problem, GLP_MAX); }

  /* Num rows and num columns */
  glp_add_rows(problem, numrows);
  glp_add_cols(problem, numcols);

  /* Assumed to be SV = 0, not SV = anything else (could change this if we need it pretty easily) */
  for(int i=1; i<numrows+1; i++) {    glp_set_row_bnds(problem, i, GLP_FX, 0.0f, 0.0f);   }
  
  /* Columns (reactions) get bounded according to the assigned LB and UB */
  for(int i=1; i<numcols+1; i++) {
    if( rougheq(lb[i] - ub[i], 0.0f, _db.FLUX_CUTOFF) == 1 ) {  glp_set_col_bnds(problem, i, GLP_FX, lb[i], lb[i]); }
    else { glp_set_col_bnds(problem, i, GLP_DB, lb[i], ub[i]); }
  }

  /* Objective function - all zeros except for the stated objectives, which have the desired coefficients */
  vector<int>::iterator it;
  for(int i=1; i<numcols+1; i++) {
    it = find(objIdx.begin(), objIdx.end(), i);
    if(it == objIdx.end()) { glp_set_obj_coef(problem, i, 0); }
    else {
      glp_set_obj_coef(problem, i, objCoef[(int)(it-objIdx.begin())]); 
    }
  }

  /* The S matrix itself */
  glp_load_matrix(problem, totalDataSize, ia, ja, ar);

  int (*func)(void*, const char *) = &suppressGLPKOutput;
  if(!_db.DEBUGFBA) { glp_term_hook(func, NULL); }

}

/******* All of these functions expect one-based indexes (one reason I made them private... )*******/
void GLPKDATA::changeLb(double newLb, int idx) {
  lb[idx] = newLb;
}
void GLPKDATA::changeUb(double newUb, int idx) {
  ub[idx] = newUb;
}
void GLPKDATA::changeObjective(vector<int> newIdx, vector<double> newCoeffs) {
  objIdx = newIdx;
  objCoef = newCoeffs;
}
void GLPKDATA::validateSense(int sense) {
  assert(sense == -1 | sense == 1);
}
int GLPKDATA::suppressGLPKOutput(void *info, const char *s) {
  /* GLPK doesn't print if this function returns 1 so that's what I do */
  return 1;
}

/* (Re)-initialize based on the given data 
 I did nto pass by reference because it causes problems with re-initializign from values already in the class (causes an easy bug to make) */
void GLPKDATA::initialize(const METSPACE &metspace, const RXNSPACE &rxnspace, vector<int> objId, vector<double> objCoeff, int sense) {

  rxnsUsed = rxnspace;
  metsUsed = metspace;

  if(objId.size() != objCoeff.size()) { printf("ERROR: In initializing GLPKDATA, provided objective IDs and objective Coefficients did not have the same size!\n"); assert(false); }
  
  validateSense(sense);
  objSense = sense;

  /****** Initialize memory ******/
  numrows = metspace.mets.size();  numcols = rxnspace.rxns.size(); 
  lb = (double*) malloc(sizeof(double) * (numcols + 1));
  ub = (double*) malloc(sizeof(double) * (numcols + 1));

  /* Total number of stoich entries - plus the one needed because GLPK starts at 1 instead of 0 */
  int total = 1;
  for(int i=0; i<numcols; i++) { total += rxnspace.rxns[i].stoich.size();  }

  ia = (int*) malloc( sizeof(int)*total);
  ja = (int*) malloc( sizeof(int)*total);
  ar = (double*) malloc( sizeof(double)*total);

  /************** Fill up data ******************/

  /* Objective info */
  objIdx.clear(); objCoef.clear();
  for(int i=0; i<objId.size(); i++) { objIdx.push_back(rxnspace.idxFromId(objId[i]) + 1);  }
  this->objCoef = objCoeff;

  /* Fill LB and UB from reaction data */
  for(int i=0; i<numcols; i++) { lb[i+1] = rxnspace.rxns[i].lb; ub[i+1] = rxnspace.rxns[i].ub; }

  /* Fill up ia (row counter), ja (column counter) and ar (stoich coeff) from the stoich data */
  int counter=0;
  for(int j=0; j<numcols; j++) {
    for(int i=0; i < rxnspace.rxns[j].stoich.size(); i++){
      /* Don't allow 0's to mess things up... */
      if(rxnspace.rxns[j].stoich[i].rxn_coeff < 1E-8 && rxnspace.rxns[j].stoich[i].rxn_coeff > -1E-8) { continue; }
      int idx = counter + 1;
      ia[idx] = metspace.idxFromId(rxnspace.rxns[j].stoich[i].met_id)+1; // Rows: Metabolites
      ja[idx] = j+1; // Columns = reactions
      ar[idx] = rxnspace.rxns[j].stoich[i].rxn_coeff;
      counter++;
    }
  }
  totalDataSize = counter;
  problem = glp_create_prob();
}

/* Print out all those lovely private variables */
void GLPKDATA::printPrivateStuff() {
  for(int i=1; i<numcols + 1; i++) {
    printf("REACTION: %s ... ", rxnsUsed.rxns[i-1].name);
    printf("LB = %4.3f; UB = %4.3f \n", lb[i], ub[i]);
  }
  for(int i=0; i< totalDataSize; i++) {
    printf("Reaction INDEX %d (NAME: %s ) and metabolite INDEX %d (NAME: %s) had coefficient %4.3f\n", 
	   ja[i+1]-1  , rxnsUsed.rxns[ja[i+1]-1].name, ia[i+1]-1, metsUsed.mets[ia[i+1] - 1].name, ar[i+1]);
  }
  printf("OBJECTIVES:\n");
  for(int i=0; i<objIdx.size(); i++){
    printf("%s (coefficient = %4.3f)\n", rxnsUsed.rxns[objIdx[i]-1].name, objCoef[i]);
  }

}

/* Print what numerical error GLPK gave us */
void GLPKDATA::printGlpkError(int errorCode) {
  printf("GLPK error code %d:\n", errorCode);
  if(errorCode == 0) {
    printf("SUCCESSFUL SOLUTION\n");
  }
  if(errorCode == GLP_EBADB) {
    printf("INITIAL BASIS INVALID\n");
  }
  if(errorCode == GLP_ESING) {
    printf("SINGULAR BASIS\n");
  }
  if(errorCode == GLP_ECOND) {
    printf("ILL-CONDITIONED BASIS\n");
  }
  if(errorCode == GLP_EBOUND) {
    printf("INVALID BOUNDS ON VARIABLES\n");
  }
  if(errorCode == GLP_EFAIL) {
    printf("SOLVER FAILED (generic) \n");
  }
  if(errorCode == GLP_EOBJLL) {
    printf("DUAL SOLUTION DECREASING WITHOUT BOUND\n");
  }
  if(errorCode == GLP_EOBJUL) {
    printf("DUAL SOLUTION INCREASING WITHOUT BOUND\n");
  }
  if(errorCode == GLP_EITLIM) {
    printf("ITERATION LIMIT EXCEEDED\n");
  }
  if(errorCode == GLP_ETMLIM) {
    printf("TIME LIMIT EXCEEDED\n");
  }
  if(errorCode == GLP_ENOPFS) {
    printf("NO PRIMAL SOLUTION [dual unbounded] (LP presolver only)\n");
  }
  if(errorCode == GLP_ENODFS) {
    printf("NO DUAL SOLUTION [primal unbounded] (LP presolver only)\n");
  }
}

