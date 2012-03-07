#include "DataStructures.h"
#include "TableLoader.h"
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>

using std::vector;
using std::map;

/* THe reaction table is the following format:
   RXN name \t RXN ID \t REVERSIBILITY \t MET name \t MET ID \t RXN_COEFF \t KEEP_IN_SECONDARY 

   The big thing missing here is the probabilities - those should probably go into a separate table (organism specific)?
   Or maybe we can just calculate them for every organism and dump them into a database.
*/
PROBLEM readReactionTable(const char* filename) {
  PROBLEM result;
  RXNSPACE &fullrxn = result.fullrxns;
  METSPACE &fullmets = result.metabolites;

  FILE* fid = fopen(filename, "r");
  /* Do not use while while(!feof()) because it causes duplicate reading of the last line 
     which is very bad for FBA (causes duplicate entries in the matrix) and also bad in general */
  while(1) {
    char rxnName[64];
    int rxnId;
    int net_reversible;
    char metName[64];
    int metId;
    char rxnCoeffStr[64];
    double rxn_coeff;
    int secondary; /* 1 if the metabolite is a secondary in that reaction */

    /* Note to self - when using fscanf do NOT attempt to specify a precision or it will get very confused */
    int status = fscanf(fid, "%s%d%d%s%d%s%d",
	   rxnName, &rxnId, &net_reversible, metName, &metId, rxnCoeffStr, &secondary);
    if(status == EOF) { break; }

    rxn_coeff = atof(rxnCoeffStr);

    STOICH curStoich;
    sprintf(curStoich.met_name, "%s", metName);
    curStoich.met_id = metId;
    curStoich.rxn_coeff = rxn_coeff;

    if(fullrxn.idIn(rxnId)) {
      fullrxn.rxnPtrFromId(rxnId)->stoich.push_back(curStoich);
      if(!secondary) {
	fullrxn.rxnPtrFromId(rxnId)->stoich_part.push_back(curStoich);
      }
    } else {
      REACTION newrxn;
      newrxn.id = rxnId;
      newrxn.init_reversible = net_reversible;
      newrxn.net_reversible = net_reversible;
      sprintf(newrxn.name, "%s", rxnName);      
      newrxn.stoich.push_back(curStoich);
      if(!secondary) {
	newrxn.stoich_part.push_back(curStoich);
      }
      fullrxn.addReaction(newrxn);
      /* Set lb and ub appropriately according to the chosen reversibility... */
      fullrxn.changeReversibility(rxnId, net_reversible);
    }
    if(!fullmets.idIn(metId)) {
      /* Note - because we're curating now, the secondary and secondary_pair fields became useless */
      METABOLITE newmet;
      newmet.id = metId;
      sprintf(newmet.name, "%s", metName);
      fullmets.addMetabolite(newmet);
    }
    /* This is needed to ensure that anything that is treated as a secondary potentially gets a magic entrance. */
    if(secondary) {
      fullmets.metPtrFromId(metId)->secondary_lone = 1;
    }
  }

  printf("Number of fullrxns: %d\n Number of Metabolites: %d\n", (int)result.fullrxns.rxns.size(), (int)result.metabolites.mets.size());
  fclose(fid);

  /* TODO - fill all the other random crap here (is it a transporter or not? is it an exchange reaction or not? Synonym reactions. Etc... */
  return result;
}

/* Likelihood table is very simple:
   Reaction \t  Likelihood

   Any reactions without likelihood assigned will have a default likelihood of -5 (due to the reaction constructor)
   Throws a warning for any reactions that are not already in the ProblemSpace - should help weed out inconsistencies... */
void readLikelihoodTable(PROBLEM &ProblemSpace, const char* filename) {

  FILE* fid = fopen(filename, "r");
  while(1) {
    char rxnName[64];
    char likelihoodString[64];
    double likelihood;
    int status = fscanf(fid, "%s%s", rxnName, likelihoodString);
    if(status == EOF) { break; }
    likelihood = atof(likelihoodString);

    int rxnId = rxnByName(ProblemSpace.fullrxns, rxnName);
    if(rxnId == -1) { 
      printf("WARNING: Reaction %s from likelihood table not found in the reaction database\n", rxnName);
      continue; }

    ProblemSpace.fullrxns.rxnPtrFromId(rxnId)->init_likelihood = likelihood;
  }
  fclose(fid);
  return;
}

/* Find a reaction by name and return the ID (-1 if none is found) */
int rxnByName(const RXNSPACE &rxnspace, const char* name) {
  for(int i=0; i<rxnspace.rxns.size(); i++) {
    if(strcmp(rxnspace.rxns[i].name, name) == 0) {
      return rxnspace.rxns[i].id;
    }
  }
  return -1;
}
