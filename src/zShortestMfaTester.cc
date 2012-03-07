#include "DataStructures.h"
#include "Exchanges.h"
#include "ETC.h"
#include "Grow.h"
#include "MyConstants.h"
#include "pathUtils.h"
#include "Paths2Model.h"
#include "Printers.h"
#include "shortestPath.h"
#include "RunK.h"
#include "visual01.h"
#include "XML_loader.h"
#include "Annotations.h"

#include <map>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using std::vector;
using std::map;
using std::string;

void writeMilpCplex(const RXNSPACE &space, const METSPACE &metspace, int targetMetId, const char* filename);
void writeMilpSol(const RXNSPACE &space, vector<double> &fluxResult, const char* filename);
void makeMagicExits(RXNSPACE &RxnSpace, const METSPACE &exitsToMake);
void adjustMagicLikelihood(RXNSPACE &rxnspace, const GROWTH &whichCondition);
RXNSPACE splitReactions(const RXNSPACE &rxnspace);

int main(int argc, char *argv[]) {
  printf("\n\n");

  PROBLEM ProblemSpace;

  printf("Setting up problem inputs...\n");fflush(stdout);
  InputSetup(argc,argv,ProblemSpace);
  printf("...done \n");fflush(stdout);

  makeMagicExits(ProblemSpace.fullrxns, ProblemSpace.metabolites);

  for(int k=0; k<ProblemSpace.growth.size(); k++) {

    setSpecificGrowthConditions(ProblemSpace, ProblemSpace.growth[k]);

    /* Split reactions into forward and reverse parts */
    RXNSPACE revSplit = splitReactions(ProblemSpace.fullrxns);
    METSPACE metSplit = ProblemSpace.metabolites;
    calcMetRxnRelations(revSplit, metSplit);

    adjustMagicLikelihood(revSplit, ProblemSpace.growth[k]);

    /* Heavy penalty on the magic exits and entrances. */
    adjustLikelihoods(revSplit.rxns, 1.0f, -1000.0f, 1.1f, -3.0f, true);

    char rxnNameFile[64];
    sprintf(rxnNameFile, "Reaction_list_growth_%d", k);
    FILE* fid = fopen(rxnNameFile, "w");
    char rxnString[4096];
    for(int i=0; i<revSplit.rxns.size(); i++) {
      rxnString[0] = '\0';
      printRxnFormula(revSplit.rxns[i], rxnString, false);
      fprintf(fid, "%d\t%s\t%s\n", i, revSplit.rxns[i].name, rxnString);      
    }
    fclose(fid);
    
    /* Write a CPLEX file for each biomass component */  
    for(int i=0; i<ProblemSpace.growth[k].biomass.size(); i++) {
      printf("Working on biomass component number %d...\n", i);
      int bmId = ProblemSpace.growth[k].biomass[i].met_id;
      char nm[128];
      sprintf(nm, "%s_growth%d%s", metSplit.metFromId(bmId).name, k, ".cplex");
      writeMilpCplex(revSplit, metSplit, bmId, nm);
      
    /* Come up with a solution with which to seed the MILP and write that to a file */
    /*    vector<int> obj; obj.push_back(FindExchange4Metabolite(revSplit.rxns, bmId));
    vector<double> objCoeff; objCoeff.push_back(1);
    
    GLPKDATA dat(revSplit, metSplit, obj, objCoeff, 1);
    vector<double> fluxResult = dat.FBA_SOLVE();
    sprintf(nm, "%s%s", metSplit.metFromId(bmId).name, ".sol");
    writeMilpSol(revSplit, fluxResult, nm); */
    }
  }
  return 0;
}

/* Split all reactions in inputSpace into forward and reverse parts
   so that all fluxes are strictly positive. */
RXNSPACE splitReactions(const RXNSPACE &rxnspace) {
  RXNSPACE revSplit;
  for(int i=0; i<rxnspace.rxns.size(); i++) {
    REACTION tmp = rxnspace.rxns[i];
    tmp.revPair = -1;
    /* FIXME - should I attempt to put the opposite direction of all the other rxns here too? (with a penalty) */
    if(tmp.init_reversible == 1) {
      revSplit.addReaction(tmp);
    }
    /* Ensure that the reaction only goes in one direction */
    if(tmp.init_reversible == -1) {
      REACTION newRxn = tmp;
      newRxn.id += _db.REVFACTOR;
      newRxn.stoich = tmp.stoich;
      for(int j=0; j<newRxn.stoich.size(); j++) {
	newRxn.stoich[j].rxn_coeff *= -1;
      }
      newRxn.ub = -newRxn.lb;
      newRxn.lb = 0.0f;
      revSplit.addReaction(newRxn);
    }
    
    if(tmp.init_reversible == 0) {
      REACTION newRxn = tmp;
      newRxn.revPair = tmp.id;
      newRxn.id += _db.REVFACTOR;
      newRxn.stoich = tmp.stoich;
      for(int j=0; j<newRxn.stoich.size(); j++) {
	newRxn.stoich[j].rxn_coeff *= -1;
      }
      newRxn.ub = -newRxn.lb;
      newRxn.lb = 0.0f;
      revSplit.addReaction(newRxn);
      /* Also put in the forward reaction, but constrained to only go forward 
        I only assign the revPair to ONE of the reactions so that we don't get
	duplicate constraints (makes GLPK's life easier) */
      newRxn = tmp;
      newRxn.lb = 0.0f;
      revSplit.addReaction(newRxn);
    }
  }

  return revSplit;

}

void adjustMagicLikelihood(RXNSPACE &rxnspace, const GROWTH &whichCondition) {
  /* Make exchanges for things that aren't in the media / growth exits for things that aren't in biomass very unlikely (must have a penalty > all legit reactions combined to be effective) */
  for(int i=0; i<rxnspace.rxns.size(); i++) {
    if(rxnspace.rxns[i].isExchange) {
      bool ok = false;
      for(int j=0; j<whichCondition.media.size(); j++) {
	if(whichCondition.media[j].id == rxnspace.rxns[i].stoich[0].met_id) { ok = true; break; }
      }
      for(int j=0; j < whichCondition.biomass.size(); j++) {
	if(whichCondition.biomass[j].met_id == rxnspace.rxns[i].stoich[0].met_id) { ok = true; break; }
      }
      if(!ok) { rxnspace.rxns[i].init_likelihood = -3; }
    }
  }
}

/* Generate magic exits in RxnSpace for all metabolies in exitsToMake */
void makeMagicExits(RXNSPACE &RxnSpace, const METSPACE &exitsToMake) {
  /* Generate magic exits for everything... */
  for(int i=0; i<exitsToMake.mets.size(); i++) {
    int id = FindExchange4Metabolite(RxnSpace.rxns, exitsToMake.mets[i].id);
    if(id == -1) {
      /* Only allow exits - not entrances, since this formulation shouldn't require them. This drastically reduces the number
	 of integer variables and the size of the search space... */
      REACTION GrExit = GrowthExit(RxnSpace.rxns, exitsToMake.mets[i].id, 1, 1000, exitsToMake.mets[i].name);
      /* We want it to stay as -1 because the direction of exiting is positive and we need positive fluxes... 
       The below code will reverse the stoich for the reverse direction if needed. */
      //      GrExit.stoich[0].rxn_coeff = 1;
      RxnSpace.addReaction(GrExit);
    }
  }
}

/* Write a SOL file with an initial integer solution (the only thing it doesn't follow necessarily is the zi + zj constraint prohibiting loops... */
void writeMilpSol(const RXNSPACE &space, vector<double> &fluxResult, const char* filename) {
  FILE* fid = fopen(filename, "w");

  fprintf(fid, "<?xml version = \"1.0\" standalone=\"yes\"?>\n");
  fprintf(fid, "<?xml-stylesheet href=\"https://www.ilog.com/products/cplex/xmlv1.0/solution.xsl\" type=\"text/xsl\"?>\n");

  fprintf(fid, "<CPLEXSolution version=\"1.1\">\n");
  fprintf(fid, "<variables>\n");

  /* Eliminate trivial (same rxn in 2 directions) flux loops by making the fluxResult = 0
   (in the direction of less flux, if there is a net flux in a particular direction) */
  for(int i=0; i<space.rxns.size(); i++) {
    if(space.rxns[i].revPair != -1) {
      double diff = fluxResult[i] - fluxResult[space.idxFromId(space.rxns[i].revPair)];
      if(rougheq(diff, 0.0f, 0.000001)) { 
	fluxResult[i] = 0.0f;
	fluxResult[space.idxFromId(space.rxns[i].revPair)] = 0.0f;
      } else {
	if(diff > 0.0f) {
	  fluxResult[i] = diff;
	  fluxResult[space.idxFromId(space.rxns[i].revPair)] = 0.0f;
	} else {
	  fluxResult[i] = 0.0f;
	  fluxResult[space.idxFromId(space.rxns[i].revPair)] = -diff;
	}
      }
    }
  }

  printDoubleVector_rxns(space, fluxResult);

  /* No need to assign indexes. CPLEX just matches things up by name (which must match, but since the index is the same in both of the writer functions it's fine) */
  for(int i=0; i<space.rxns.size(); i++) {
    /* Zi = 0 if flux i = 0 and 1 otherwise */
    int zVal;
    double vVal;
    if(rougheq(fluxResult[i], 0.0f, 0.000001) == 1) {
      zVal = 0;
    } else { 
      zVal = 1; 
    } 
    fprintf(fid, "\t<variable name=\"z%d\" value=\"%d\" />\n", i, zVal);
    fprintf(fid, "\t<variable name=\"v%d\" value=\"%1.8f\" />\n", i, fluxResult[i]);
  }

  fprintf(fid, "</variables>\n");
  fprintf(fid, "</CPLEXSolution>");

  fclose(fid);
}

/* Write a CPLEX-formatted file containing all the necessary constraints (so that I don't have to think about how to generate a giant matrix 
   with all this stuff in it... just to test to see if this works mostly) 
*/
void writeMilpCplex(const RXNSPACE &space, const METSPACE &metspace, int targetMetId, const char* filename) {
  FILE* fid = fopen(filename, "w");
  double T = 0.00001f; /* Flux threshold (below this it is considered 0) */
  double T2 = 0.1f; /* Optimality threshold (below this it is not a solution we want to consider) */
  double M = 10000.0f; /* Maximum flux evarrrrr */

  /* Objective function: sum(zi * likelihood_i) 
   The maximum line length is 500 so if I make a new line every 10 I should be OK (no identifiers are longer than 6 digits + 2 spaces + a sign = 9 */
  fprintf(fid, "MINIMIZE\n");
  for(int i=0; i<space.rxns.size(); i++) {
    if(space.rxns[i].current_likelihood < 0.0f) { printf("ERROR: Reaction %s has current_likelhiood < 0\n", space.rxns[i].name); }
    //    fprintf(fid, "+ %1.3f z%d ", space.rxns[i].current_likelihood, i);
    fprintf(fid, "+ 1 z%d ", i);
    if( i != 0 && (i/10) * 10 == i ) { fprintf(fid, "\n\t"); }
  }
  fprintf(fid, "\n");

  /* Constraints... */
  fprintf(fid, "SUBJECT TO\n");
  /* Sv = 0 */
  for(int i=0; i<metspace.mets.size(); i++) {
    for(int j=0; j<metspace.mets[i].rxnsInvolved_nosec.size(); j++) {
      int rxnIdx = space.idxFromId(metspace.mets[i].rxnsInvolved_nosec[j]);
      for(int k=0; k<space.rxns[rxnIdx].stoich.size(); k++) {
	if(space.rxns[rxnIdx].stoich[k].met_id == metspace.mets[i].id) {
	  if(space.rxns[rxnIdx].stoich[k].rxn_coeff < 0.0f) { 
	    fprintf(fid, "- %1.8f v%d ", -space.rxns[rxnIdx].stoich[k].rxn_coeff, rxnIdx);
	  } else {
	    fprintf(fid, "+ %1.8f v%d ",  space.rxns[rxnIdx].stoich[k].rxn_coeff, rxnIdx);
	  }
	}
      }
      if( j != 0 && (j/10) * 10 == j ) { fprintf(fid, "\n\t"); }
    }
    fprintf(fid, " = 0\n");
  }
  /* Vj > 0 --> j is the exchange reaction for the target metabolite */
  int ExcId = FindExchange4Metabolite(space.rxns, targetMetId);
  int ExcIdx = space.idxFromId(ExcId);
  fprintf(fid, "v%d >= %1.3f\n", ExcIdx, T2);
  /* For the same reaction, z_forward + z_backward <= 1 */
  for(int i=0; i<space.rxns.size(); i++) {
    if(space.rxns[i].revPair == -1) { continue; }
    fprintf(fid, "z%d + z%d <= 1\n", i, space.idxFromId(space.rxns[i].revPair));
  }
  /* zj * T <= vj [if flux = 0 then zj = 0] 
     zj <= M * vj [if flux != 0 then zj = 1] 
    Note CPLEX requires the equations in standard form (constant on the right hand side) */
  for(int i=0; i<space.rxns.size(); i++) {
    fprintf(fid, "%1.8f z%d - v%d <= 0\n", T, i, i);
    fprintf(fid, "v%d - %1.8f z%d <= 0\n", i, M, i);
  }

  /* LB and UB */
  fprintf(fid, "BOUNDS\n");
  for(int i=0; i<space.rxns.size(); i++) {
    fprintf(fid, "%1.3f <= v%d <= %1.3f\n", space.rxns[i].lb, i, space.rxns[i].ub);
  }

  fprintf(fid, "BINARY\n");
  for(int i=0; i<space.rxns.size(); i++) {
    fprintf(fid, "z%d ", i);
    if( i != 0 && (i/10) * 10 == i ) { fprintf(fid, "\n\t"); }
  }
  fprintf(fid, "\n");

  /* Don't forget this! Otherwise whatever is reading it will be very angry that you don't tell it it is done... */
  fprintf(fid, "END");
  fclose(fid);
}
