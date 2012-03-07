#include "MyConstants.h"
#include "pathUtils.h"
#include "XML_loader.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iterator>
#include <map>
#include <vector>

using std::vector;
using std::map;

bool isExternalMet(const char *name, const char *E_tag) {
  /* I assume that E_tag cannot be more than 8 characters long */
  char tempS[8] = {0};
  if(strlen(name) < strlen(E_tag)) { return false; }
  /* Check for E_tag at the end of the metabolite name */
  strncpy(tempS,name+((int)strlen(name)-strlen(E_tag)),strlen(E_tag));
  if(strcmp(tempS, E_tag) != 0 ) {  return false;  }
  return true;
}

bool anyVec(const vector<bool> &vec) {
  for(int i=0; i<vec.size(); i++) {
    if(vec[i]) {
      return true;
    }
  }
  return false;
}

/* Concatenate two vectors - don't sort or make unique (will have to do that in other functions if desired) */
/* Adds all of vector2 to the end of vector1 */
void catVectors1(vector<int> &vector1, const vector<int> &vector2) {
  int i;
  for(i=0;i<vector2.size();i++){
    vector1.push_back(vector2[i]);
  }
  return;
}

/* Intersect function for two vectors (does not assume  A and B are already unique)  */
vector<int> custom_intersect(vector<int> A, vector<int> B) {
  vector<int>::iterator iter;
  vector<int> result(A.size() + B.size());
  unsigned int i;
  custom_unique(A);
  custom_unique(B);
  iter = set_intersection(A.begin(), A.end(), B.begin(), B.end(), result.begin());
  vector<int> inter(iter - result.begin());
  for(i= 0; i< inter.size(); i++){
    inter[i] = result[i];}

  return inter;
}

vector<int> setdiff2(vector<int> A, vector<int> B) {
  if(B.size()<1){
    return A;}
  int i, j(0);	
  
  sort(A.begin(), A.end());
  sort(B.begin(), B.end());
  custom_unique(A);
  custom_unique(B);
  vector<int>::iterator iter;
  vector<int> result(A.size()+B.size());

  iter = set_difference(A.begin(), A.end(), B.begin(), B.end(), result.begin());

  vector<int> diff(iter-result.begin());
  for(i=0;i<diff.size();i++) {
    diff[i]=result[i];
  }

  return diff;
}

void setdiff(vector<int> &A, vector<int> &B, vector<int> &diff) {
  diff = setdiff2(A, B);
  return;
}

void setdiff1(vector<int> &A, const vector<int> &B) {
  vector<int> diff = setdiff2(A, B);
  A = diff;
  return;
}

int instoich(int id, const vector<STOICH> &list){
  unsigned int i;
  for(i=0;i<list.size();i++){
    if(id==list[i].met_id){ return i; }
  }
  return -1;
}

/* Convert a metabolite name to an ID. Returns -1 on failure */
int Name2Ids(const vector<METABOLITE> &metabolite, const char *met_name){
  for(int i=0;i<metabolite.size();i++){
    if(strcmp(met_name,metabolite[i].name)==0){
      return metabolite[i].id;
    }
  }
  return -1;
}

int Name2Ids(const vector<REACTION> &reaction, const char *rxn_name) {
  for(int i=0;i<reaction.size();i++){
    if(strcmp(rxn_name,reaction[i].name)==0){
      return reaction[i].id;
    }
  }
  return -1;
}

/* 0 means they are different and 1 means they are the same (based on stoich_part) - used to synonymize */
int diff2rxns(const REACTION &one, const REACTION &two){
  /* Check sizes */
  int one_stoich_size = one.stoich_part.size();
  int two_stoich_size = two.stoich_part.size();
  if(one_stoich_size==0 || two_stoich_size == 0){ return 0;}
  if(one_stoich_size!=two_stoich_size){ return 0;}
  int i,j,k[one_stoich_size];
  /* If metabolites not the same, return difference mark (0) */
  for(i=0;i<one.stoich_part.size();i++){
    k[i] = instoich(one.stoich_part[i].met_id,two.stoich_part);
    if(k[i]==-1){return 0;}
  }
  /* Make sure reactions have same ratio of reactants */
  double kk[one_stoich_size];
  for(i=0;i<one_stoich_size;i++){
    kk[i] = (double)one.stoich_part[i].rxn_coeff/(double)two.stoich_part[k[i]].rxn_coeff;
  }
  double ul = kk[0] + 0.1;
  double ll = kk[0] - 0.1;
  for(i=0;i<one_stoich_size;i++){
    if(kk[i]>ul || kk[i]<ll){ return 0;}
  }
  return 1;
}

/* Extract a smaller subset of reactions that we're interested in, given their indices and a larger list */
/* Maintains order of rxnsToSearchIds */
void pullOutRxnsbyIds(const RXNSPACE &rxnspace, const vector<int> &rxnsToSearchIds, RXNSPACE &result) {
  result.clear();
  for(int i=0;i<rxnsToSearchIds.size();i++) {
    if(!result.idIn(rxnsToSearchIds[i])){
      result.addReaction(rxnspace.rxnFromId(rxnsToSearchIds[i]));
    }
  }
  return;
}

//For synrxns/synrxnsR simultaneiously
void pullOutRxnsbyIds(const PROBLEM &ProblemSpace, const vector<int> &rxnsToSearchIds, RXNSPACE &result) {
  result.clear();
  for(int i=0;i<rxnsToSearchIds.size();i++) {
    if(rxnsToSearchIds[i] >= _db.REVFACTOR && rxnsToSearchIds[i] < _db.REVFACTOR + _db.SYNFACTOR + _db.MINFACTORSPACING){
      if(!result.idIn(rxnsToSearchIds[i])){
	result.addReaction(ProblemSpace.synrxnsR.rxnFromId(rxnsToSearchIds[i]));
      }
    }
    else{
      if(!result.idIn(rxnsToSearchIds[i])){
	result.addReaction(ProblemSpace.synrxns.rxnFromId(rxnsToSearchIds[i]));
      }
    }
  }
  return;
}

/* Pull out a METSPACE containing all of the metabolites associated with the reactions in RXNSPACE 
   (the result is a subset of "metspace"). Asserts an error if any metabolites in rxnspace are not in the metspace.*/
void pullOutMets(const RXNSPACE &rxnspace, const METSPACE &metspace, METSPACE &result) {
  result.clear();
  vector<int> metIds;
  const vector<REACTION> &rxnList = rxnspace.rxns; 
  for(int i=0;i<rxnList.size();i++) {
    for(int j=0;j<rxnList[i].stoich.size();j++) {
      metIds.push_back(rxnList[i].stoich[j].met_id);
    }
  }
  custom_unique(metIds);
  for(int i=0;i<metIds.size();i++) {
    result.addMetabolite(metspace.metFromId(metIds[i]));
  }
  return;
}

/* Based on cofactors removed from reactions in ProblemSpace.fullrxns (and placed in stoich_part)
   make a list of synonymous reactinos and dump it to ProblemSpace.synrxns (in stoich_part). */

void MakeSynList(PROBLEM &ProblemSpace){

  REACTION rtemp;
  RXNSPACE &fullspace = ProblemSpace.fullrxns;
  RXNSPACE &synspace  = ProblemSpace.synrxns;
  int *ToProcess = (int *) malloc(sizeof(int) * fullspace.rxns.size());
  for(int i=0;i<fullspace.rxns.size();i++){ToProcess[i] = 1;}
  for(int i=0;i<fullspace.rxns.size();i++){    
    if(fullspace.rxns[i].init_likelihood < -0.9f && fullspace.rxns[i].init_likelihood > -1.1f){
      ToProcess[i] = 0;}
  }
  /* Attempt to speed this up a bit... this should be more than enough memory */
  synspace.rxns.reserve(fullspace.rxns.size());

  /* Run through all pairs and match them up in SynTemp  - note the way this is set up also matches up each reaction with itself */
  for(int i=0;i<fullspace.rxns.size();i++){
    if(ToProcess[i]){
      synspace.rxns.push_back(fullspace.rxns[i]);
      for(int j=i;j<fullspace.rxns.size();j++){
	/* Check if the reactions are the same. If they are, add to synlist */
	if(ToProcess[j] && diff2rxns(fullspace.rxns[i],fullspace.rxns[j])){
	  synspace.rxns.back().syn.push_back(fullspace.rxns[j].id);
	  ToProcess[j] = 0;
	}
      }
    }
  }
  
  for(int i=0;i<synspace.rxns.size();i++){
    if(synspace.rxns[i].syn.size()>1){ 
      synspace.rxns[i].id += _db.SYNFACTOR;

      /* Merge likelihood values.
	 FIXME: Need to add likelihood priorities here. For now I assume there's never more than one - thing in a list of
	 synonyms. However, -5 should be handled correctly now... */
      vector<double> current_likelihood;
      for(int j=0;j<synspace.rxns[i].syn.size();j++){
	double curFlag = fullspace.rxnFromId(synspace.rxns[i].syn[j]).init_likelihood;
	if(curFlag < 0 && rougheq(curFlag,-5)!=1 ) {
	  synspace.rxns[i].init_likelihood = curFlag;
	  current_likelihood.clear(); current_likelihood.push_back(curFlag);
	  break;
	} else if(rougheq(curFlag, -5)==1) {
	  continue;
	} else {
	  current_likelihood.push_back(curFlag);
	}
      }

      if (current_likelihood.empty() ) {
	synspace.rxns[i].init_likelihood = -5;
      } else {
	double sum(0.0f);
	for(int j=0; j<current_likelihood.size(); j++) {
	  sum+=current_likelihood[j];
	}
	synspace.rxns[i].init_likelihood = sum/(double)current_likelihood.size();
      }

      REACTION firstSyn = fullspace.rxnFromId(synspace.rxns[i].syn[0]);
      synspace.rxns[i].stoich_part = firstSyn.stoich_part;

      /* Merge Reversibilities */
      int firstMetId = firstSyn.stoich_part[0].met_id;
      double firstCoeff = firstSyn.stoich_part[0].rxn_coeff;
      int sgn;
      /* Technically it could be 0 too but we don't care, because that means net_reversible is 0 and we'll make the rev 0 by default */
      if(firstCoeff * (double)firstSyn.net_reversible < 0.0f) { sgn = -1; } else { sgn = 1; }
      int rev = firstSyn.net_reversible; /* We want to model the reversibility on the FIRST of the synrxns */

      for(int j=0;j<synspace.rxns[i].syn.size();j++){
	REACTION curRxn = fullspace.rxnFromId(synspace.rxns[i].syn[j]);
	if(curRxn.net_reversible == 0) { rev = 0; break; }

	int curRev = curRxn.net_reversible;
	double curCoeff;
	for(int k=0; k<curRxn.stoich_part.size(); k++) {
	  if(curRxn.stoich_part[k].met_id == firstMetId) { 
	    curCoeff = curRxn.stoich_part[k].rxn_coeff;
	    break;
	  }
	}
	/* We always want the reversibility * coefficient to have the same sign 
	   or else the reaction is essentially reversible. */ 
	int curSgn;
	if(curCoeff * (double)curRev < 0.0f) { curSgn = -1; } else { curSgn = 1; }
	if(curSgn != sgn) { rev = 0; break; }
      }
      synspace.rxns[i].net_reversible = rev;
    }
  }
  /* We shouldn't really need to do this but I don't have any desire to rewrite this code, since it took forever tog et right the first time.. so here it is. */
  synspace.rxnMap();
  free(ToProcess);
  return;
}

/* Version of getProducts for a single REACTION rather than a vector of REACTIONS
 Returns product IDs (not indexes) and empty upon not finding reactant "reactantId" in the reaction */
void getProducts(const REACTION &rxn, int reactantId, vector<int> &res) {
  res.clear();
  int sgn(0);
  unsigned int i;
  bool productPresent = false;
  for(i=0;i<rxn.stoich_part.size();i++) 	{
    if(rxn.stoich_part[i].met_id == reactantId)  {
      if(rxn.stoich_part[i].rxn_coeff < 0)
	{ sgn = -1; }
      else { sgn = 1; }
      /* Test reversibility -use net_reversible to account for changes due to the algorithm */
      switch(rxn.net_reversible) {
      case 0: /* Reversible - always count the reactant as present */
	productPresent = true;
	break;
      case 1: /* Forward only - only count the reactant as present if the sign is -1 */
	if(sgn == -1) { productPresent = true; }
	break;
      case -1: /* Backward only - only count the reactant as present if the sign is +1 */
	if(sgn == 1) { productPresent = true; }
	break;
      default:
	printf("ERROR - undefined reversibility value %d in getProducts:\n",  rxn.net_reversible);
	throw;
      }
    }
  }
  
  if(!productPresent) {return;}
  
  /* Get and return objects with sign opposite of sgn (i.e. if sgn*coeff < 0) */
  for(i=0;i<rxn.stoich_part.size();i++) {
    if( (rxn.stoich_part[i].rxn_coeff * sgn) < 0 ) {
      res.push_back(rxn.stoich_part[i].met_id);
    }
  }
  return;
}

/* "Flatten" Path summary and do some validation too! */
vector<PATHSUMMARY> flattenPsum(const vector<vector<vector<PATHSUMMARY> > > &from) {
  vector<PATHSUMMARY> result;
  for(int i=0;i<from.size();i++) {
    for(int j=0;j<from[i].size();j++) {
      for(int k=0;k<from[i][j].size();k++) {
	PATHSUMMARY tmp = from[i][j][k];
	result.push_back(tmp);
      }
    }
  }
  return result;
}

// Takes flattened vector<PATHSUMMARY> and uniques paths across growth conditions
vector<PATHSUMMARY> uniquePsum(const vector<PATHSUMMARY> &from){
  vector<PATHSUMMARY> result;
  for(int i=0;i<from.size();i++){
    PATHSUMMARY temp = from[i];
    int bit = 0;
    for(int j=0;j<result.size();j++){
      if(temp==result[j]){ 
	bit = 1; 
	result[j].growthIdx.push_back(temp.growthIdx[0]);
	continue;
      }
    }
    if(bit==0){
      result.push_back(temp);
    }   
  }
  /*Good error check here would be to see if growthIdx is a unique vector*/

  return result;
}

int rougheq(const double one, const double two){
  if(one < two + 0.1 && one > two - 0.1){ return 1;}
  else{ return 0; }
}

/* Test if two variables are equal within a constant */
int rougheq(const double one, const double two, const double constant){
  if(one < two + constant && one > two - constant){ return 1;}
  else{ return 0; }
}

/* Returns the stoichiometric coefficient of a reaction given the reaction ID and met ID you want 
 Returns 0.0f if the metabolite is not in the reaction 

Looks for the metabolite in stoich, not stoich_part */
double rxnCoeff(const RXNSPACE &rxnspace, int rxnId, int metId) {
  REACTION rxn = rxnspace.rxnFromId(rxnId);
  for(int i=0; i<rxn.stoich.size(); i++) {
    if(rxn.stoich[i].met_id == metId) {
      return rxn.stoich[i].rxn_coeff;
    }
  }

  return 0.0f;
}

vector<int> getAllSecondaryIds(const PROBLEM &ProblemSpace) {
  vector<int> secondaryIds;
  for(int i=0; i<ProblemSpace.metabolites.mets.size(); i++) {
    if(ProblemSpace.metabolites.mets[i].isSecondary()) {
      secondaryIds.push_back(ProblemSpace.metabolites.mets[i].id);
    }
  }
  return secondaryIds;
}

/* whichOne = 1 --> fullrxns (stoich)
   whichOne = 2 --> synrxns (stoich_part) */
vector<int> getAllPathMets(const vector<PATHSUMMARY> &pList, const PROBLEM &problemSpace, int whichOne) {
  vector<int> metList;
  for(int i=0; i<pList.size(); i++) {
    for(int j=0; j<pList[i].rxnDirIds.size(); j++) {
      REACTION tmpRxn;
      switch(whichOne) {
      case 1:
	tmpRxn = problemSpace.fullrxns.rxnFromId(abs(pList[i].rxnDirIds[j]));
	for(int k=0; k<tmpRxn.stoich.size(); k++) {
	  metList.push_back(tmpRxn.stoich[k].met_id);
	}
	break;
      case 2:
	tmpRxn = problemSpace.synrxns.rxnFromId(abs(pList[i].rxnDirIds[j]));
	for(int k=0; k<tmpRxn.stoich_part.size(); k++) {
	  metList.push_back(tmpRxn.stoich[k].met_id);
	}
	break;
      default:
	printf("ERROR: Undefined option whichOne passed to getAllPathMets\n"); assert(false);
      }
    }
  }
  custom_unique(metList);
  return metList;
}

/* Get a vector of path IDs from a pathsummary list */
vector<int> getAllPathRxns(const vector<PATHSUMMARY> &pList) {
  vector<int> result;
  for(int i=0; i<pList.size(); i++) {
    catVectors1(result, pList[i].rxnDirIds);
  }
  custom_unique(result);
  return result;
}


void turnOffMagicExits(RXNSPACE &rxnspace, const vector<int> &exitsToKeep) {
  for(int i=0; i<rxnspace.rxns.size(); i++) {
    if(rxnspace.rxns[i].id < _db.BLACKMAGICFACTOR || rxnspace.rxns[i].id >= _db.BLACKMAGICFACTOR + _db.MINFACTORSPACING) { continue; }
    bool keep = false;
    for(int j=0; j<exitsToKeep.size(); j++) {
      if(exitsToKeep[j] == rxnspace.rxns[i].id) { keep = true; break; }
    }
    if(!keep) { rxnspace.rxns[i].lb = 0.0f; rxnspace.rxns[i].ub = 0.0f; }
  }
}

/* Modify the NGAM associated with the model PROBLEM
   NGAM is modified by changing the lower bound and upper bound of the reaction with name given by db.ATPM_name
   FIXME: Need to make sure to add the ATPM reaction to the model */
void modifyNGAM(ANSWER &model, double newAtpm) {
  int atpmId = Name2Ids(model.reactions.rxns, _db.ATPM_name);
  //printf("atpmid = %d\n",atpmId);                                                                                                                                                                           
  if(atpmId == -1) {
    printf("WARNING: ATPM reaction was never added to model!\n");
  }
  model.reactions.rxnPtrFromId(atpmId)->lb = newAtpm;
  model.reactions.rxnPtrFromId(atpmId)->ub = newAtpm;
}

/* Modify the growth-associated maintenance by changing the amount of ATP, ADP, H2O, Pi, and H                                                                                                                 
   by the amount specified in amountToChange. THe amount of H is identified based on the ATPM value.                                                                                                           
   If the speciifed amountToChange would make any of the coefficients change signs they are set to 0                                                                                                           
   instead and the maximum possible amountToChange is used (with a warning)                                                                                                                                    
                                                                                                                                                                                                               
amountToChange > 0 --> more ATPM                                                                                                                                                                               
amountToChange < 0 --> less ATPM */
void modifyGAM(ANSWER &model, double amountToChange) {
  int atpmId = Name2Ids(model.reactions.rxns, _db.ATPM_name);
  int hId = Name2Ids(model.metabolites.mets, _db.H_name);
  double numH = 0.0f;
  REACTION* biomass = model.reactions.rxnPtrFromId(_db.BIOMASS);
  REACTION ATPM = model.reactions.rxnFromId(atpmId);

  if(atpmId == -1) { printf("ERROR: ATPM reaction was never added to the model!\n"); assert(false); }
  if(hId == -1) { printf("ERROR: No H found for some reason in the list of model metabolites...\n"); assert(false); }

  for(int i=0; i<ATPM.stoich.size(); i++) {
    if(ATPM.stoich[i].met_id == hId) {
      if(ATPM.stoich[i].rxn_coeff < 0.0f) { numH = -ATPM.stoich[i].rxn_coeff; } else { numH = ATPM.stoich[i].rxn_coeff; }
    }
  }

  int atpId = Name2Ids(model.metabolites.mets, _db.ATP_name);
  int adpId = Name2Ids(model.metabolites.mets, _db.ADP_name);
  int piId = Name2Ids(model.metabolites.mets, _db.PI_name);
  int h2oId = Name2Ids(model.metabolites.mets, _db.H2O_name);

  if(atpId == -1 || adpId == -1 || piId == -1 || h2oId == -1) { printf("ERROR: ATP, ADP, PI, or H2O was not found in the model\n"); assert(false); }

  /* Test what the most that can be modified is and only go that far */
  /* Test what the most that can be modified is and only go that far */
  double mostCanModify(amountToChange);
  for(int i=0; i<biomass->stoich.size(); i++) {
    double tmpNum;
    if(biomass->stoich[i].met_id == adpId) { /* Product */
      tmpNum = biomass->stoich[i].rxn_coeff + amountToChange;
      if(tmpNum < 0.0f) {
        printf("WARNING: Attempt to modify biomass equation to have negative ATP maintenance\n");
        if(biomass->stoich[i].rxn_coeff < mostCanModify) { mostCanModify = biomass->stoich[i].rxn_coeff; }
      }
    }
    if(biomass->stoich[i].met_id == piId) { /* Product */
      tmpNum = biomass->stoich[i].rxn_coeff + amountToChange;
      if(tmpNum < 0.0f) {
        printf("WARNING: Attempt to modify biomass equation to have negative ATP maintenance\n");
        if(biomass->stoich[i].rxn_coeff < mostCanModify) { mostCanModify = biomass->stoich[i].rxn_coeff; }
      }
    }
    if(biomass->stoich[i].met_id == hId) { /* Product or non-existant (need numH) */
      tmpNum = biomass->stoich[i].rxn_coeff + numH*amountToChange;
      if(tmpNum < 0.0f) {
        printf("WARNING: Attempt to modify biomass equation to have negative ATP maintenance\n");
        if(numH == 0) { printf("ERROR: Consistency problem between numH and calcualtion of biomass coefficients\n"); assert(false); }
        if(biomass->stoich[i].rxn_coeff/numH < mostCanModify) { mostCanModify = biomass->stoich[i].rxn_coeff/numH; }
      }
    }
    if(biomass->stoich[i].met_id == atpId) { /* Reactant */
      tmpNum = biomass->stoich[i].rxn_coeff - amountToChange;
      if(tmpNum > 0.0f) {
        printf("WARNING: Attempt to modify biomass equation to have negative ATP maintenance - will not modify...\n");
        if( (-biomass->stoich[i].rxn_coeff) < mostCanModify) { mostCanModify = -biomass->stoich[i].rxn_coeff; }
      }
    }
    if(biomass->stoich[i].met_id == h2oId) { /* Reactant */
      tmpNum = biomass->stoich[i].rxn_coeff - amountToChange;
      if(tmpNum > 0.0f) {
        printf("WARNING: Attempt to modify biomass equation to have negative ATP maintenance - will not modify...\n");
        if( (-biomass->stoich[i].rxn_coeff) < mostCanModify) { mostCanModify = -biomass->stoich[i].rxn_coeff; }
      }
    }
  }

  /* Actually modify the coeffs */
  if(abs(amountToChange - mostCanModify) > 0.01){
    printf("Requested GAM modified by %4.3f...\n", amountToChange);
    printf("Actually modifying NGAM by %4.3f.\n", mostCanModify);
  }
  for(int i=0; i<biomass->stoich.size(); i++) {
    if(biomass->stoich[i].met_id == adpId) { biomass->stoich[i].rxn_coeff += mostCanModify; }
    if(biomass->stoich[i].met_id == piId)  { biomass->stoich[i].rxn_coeff += mostCanModify; }
    if(biomass->stoich[i].met_id == hId)   { biomass->stoich[i].rxn_coeff += (mostCanModify*numH); }
    if(biomass->stoich[i].met_id == atpId) { biomass->stoich[i].rxn_coeff -= mostCanModify; }
    if(biomass->stoich[i].met_id == h2oId) { biomass->stoich[i].rxn_coeff -= mostCanModify; }
  }

  return;
}
