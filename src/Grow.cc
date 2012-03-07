#include "DataStructures.h"
#include "Exchanges.h"
#include "Grow.h"
#include "kShortest.h"
#include "genericLinprog.h"
#include "MersenneTwister.h"
#include "MyConstants.h"
#include "pathUtils.h"
#include "Paths2Model.h"
#include "Printers.h"
#include "RunK.h"
#include "XML_loader.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

using std::map;
using std::queue;
using std::set;
using std::string;
using std::vector;

/* Wrapper script for the gapfind/gapfill algorithm that runs it on ONE chosen growth condition
   and combines the results into something usable for the genetic algorithm. 

   This will serve as the inner loop to the genetic algortihm...  */
ANSWER gapfillWrapper(const PROBLEM &problemSpace, const vector<PATHSUMMARY> &pList, const GROWTH &growth) {

  //Is this the line?
  ANSWER result;

  /* Check that the static things were filled up as required before calling this function */
  if(result.etc.empty() && result.pList.empty()) {
    printf("INTERNAL ERROR: Static variables ANSWER::etc and ANSWER::plist must be initialized before calling the inner loop (gapfillWrapper)\n");
    assert(false);
  }
  if(result.etc.empty()) { 
    printf("WARNING: No ETC chain information found - this will become an error once I figure out how to implement that within the REACTION structure...\n");
  }

  int gapfillK = _db.GAPFILL_K;

  /* Set up a network containing ALL the media conditions first (which is what this function does) 
     Magic entrances are made for all secondary lones and pairs and the network is checked for exchanges
     and transporters for all media and byproducts. Since these reactions will all be present (but possibly turned off)
     during all simulation conditions, this doesn't break anything. */
  PROBLEM baseModel = setUpGapfill(problemSpace, pList);
  
  /* Set up uptake rates and exchanges on/off for the specific growth condition passed here */
  setSpecificGrowthConditions(baseModel, growth);
  vector<GAPFILLRESULT> res = gapFindGapFill(baseModel, problemSpace, gapfillK); 
  vector<int> allEssentialExits;
  /* Identify the set of gapfill solutions that minimizes the number of essential magic exits and entrances... */
  vector<int> whichK = findSolutionsMinimizingExits(baseModel, problemSpace, res, allEssentialExits);
  custom_unique(allEssentialExits);

  printf("Final Essential exits: \n");
  for(int i=0; i<allEssentialExits.size(); i++) {
    printf("%s\t", baseModel.fullrxns.rxnFromId(allEssentialExits[i]).name);
  }
  printf("\n");
  
  /* Add the whichK to the model before running the next growth condition. Hopefully this minimizes the number of
     redundant answers... */
  for(int i=0; i<res.size(); i++) {  addGapfillResultToProblem(baseModel, problemSpace, res[i], whichK[i]);  }
  minimizeExits(baseModel);
  
  if(_db.PRINTGAPFILLRESULTS) {  PrintGapfillResult(res, problemSpace, whichK);  }

  /* Fill up ANSWER structure (note - the point of using the pointers is that all of them must point to the copy INSIDE the class */
  result.reactions = baseModel.fullrxns;
  result.metabolites = baseModel.metabolites;
  for(int i=0; i<res.size(); i++) {
    const METABOLITE* fixedPtr = result.metabolites.metPtrFromId(res[i].deadMetId);
    vector<const REACTION*> fixPtrVec;
    for(int j=0; j<res[i].deadEndSolutions[whichK[i]].size(); j++) {
      const REACTION * onePtr = result.reactions.rxnPtrFromId(res[i].deadEndSolutions[whichK[i]][j]);
      fixPtrVec.push_back(onePtr);
    }
    result.fixedGaps.push_back(fixedPtr);
    result.fixes.push_back(fixPtrVec);
  }
  for(int i=0; i<allEssentialExits.size(); i++) {
    const REACTION* essentialExit = result.reactions.rxnPtrFromId(allEssentialExits[i]);
  }

  MATLAB_out("Optimal_innerloop", result.reactions.rxns); 

  return result;
}

/* Genetic algorithm to find the best set of gapfill solutions...
 Note - this only attempts to find ONE solution for each gap. 
 However, it is possible that the algorithm gets confused because some gaps
 could be either entrance OR exit gaps. */
vector<int> findSolutionsMinimizingExits(const PROBLEM &baseModel, const PROBLEM &problemSpace, const vector<GAPFILLRESULT> &res, vector<int> &essential) {

  MTRand rng;
  vector<int> bestK;

  int sizeOfPopulation = 20; /* Number of parents to start with and size of population to maintain */
  double fracToKeep = 0.50f; /* Every time we cross over between this percent and throw out the rest...*/
  int numIterations = 5; /* For now I'll test just running it this number of times... */
  int maxNumBad = 10; /* Maximum number of times crossovers fail before doing mutations */
  int numToMutate = 8; /* Number of different gapfills to mutate per mutation round (could be this number or less) */

  /* I decided NOT to use a set here... because I need to change the scores all the time */
  vector<INNERPOPSTORE> populationSet;

  for(int i=0; i<sizeOfPopulation; i++) {
    vector<int> pop = randomK(res, rng);
    double score = -1.0f;
    INNERPOPSTORE tmp; tmp.score = score; tmp.whichK = pop;
    populationSet.push_back(tmp);
  }

  for(int iter = 0; iter < numIterations; iter++) {
    for(int i=0; i<populationSet.size(); i++) {
      if(populationSet[i].score >= -0.001f) { continue; } /* Don't re-evaluate score for the parents */
      printf("Evaluating score for population number %d in iteration %d...\n", i, iter);
      PROBLEM modified = baseModel;
      vector<int> whichK = populationSet[i].whichK;
      for(int j=0; j<res.size(); j++) { addGapfillResultToProblem(modified, problemSpace, res[j], whichK[j]); }
      populationSet[i].score = innerScore(modified, problemSpace, populationSet[i].essentialExits);
    }

    /* Keep only the top fracToKeep percent...so basically, pop_back everything beyond the cutoff 
     At most, this becomes populationSet.size(), and populationSet.begin() + populationSet.size() = populationSet.end()
    Note - lower scores are better so the default ordering is what we really want... */
    sort(populationSet.begin(), populationSet.end());
    int cutoff = fracToKeep*populationSet.size();
    while(populationSet.size() > cutoff) {  populationSet.pop_back();  }

    /* Crossover: Pick a random point in the whichK of two good solutions and cross them over like this
       1A | 1B
       2A | 2B 
       ========
       1A | 2B   and 2A | 1B
       Also keep tha parents, and remove the weakest ones (with some probability) 

       CHANGE 2-6-12 - to get more diversity I do mutation and crossover in every step. */
    while(populationSet.size() < sizeOfPopulation) {
      INNERPOPSTORE child = crossOver(populationSet, cutoff, rng);
      if(!(child.score < -1.90f)) {  
	populationSet.push_back(child);
      }

      INNERPOPSTORE child2 = mutate(populationSet, res, numToMutate, cutoff, rng);
      populationSet.push_back(child2);
    }
    
    /* Fall back on mutations if we have selected out too much diversity */
    while(populationSet.size() < sizeOfPopulation) {
      printf("Mutating...\n");
      INNERPOPSTORE child = mutate(populationSet, res, numToMutate, cutoff, rng);
      populationSet.push_back(child);
    }
    /* Make sure it's impossible to pass two identical parents to cross over */
    custom_unique(populationSet);
  }

  /* To make sure we actually get the best solution out... */
  for(int i=0; i<populationSet.size(); i++) {
    if(populationSet[i].score >= 0.0f) { continue; } /* Don't re-evaluate score for the parents */
    printf("Evaluating score for population number %d in iteration FINAL...\n", i);
    PROBLEM modified = baseModel;
    vector<int> whichK = populationSet[i].whichK;
    for(int j=0; j<res.size(); j++) { addGapfillResultToProblem(modified, problemSpace, res[j], whichK[j]); }
    populationSet[i].score = innerScore(modified, problemSpace, populationSet[i].essentialExits);
  }
  sort(populationSet.begin(), populationSet.end());

  printf("OPTIMUM FOUND - inner loop SCORE = %1.5f\n", populationSet[0].score);
  printf("WhichK:\n");
  printIntVector(populationSet[0].whichK);
  essential = populationSet[0].essentialExits;

  /* Just return the best one... */
  return populationSet[0].whichK;
}

/* Mutate a random part of a random parent and call it a new child... do this in case the crossovers fail too many times
   to keep us from getting into a pit and blowing the call stack. */
INNERPOPSTORE mutate(const vector<INNERPOPSTORE> &parents, const vector<GAPFILLRESULT> &res, int nToMutate, int mostFitCutoff, MTRand &rng) {
  int Parent = rng.randInt(mostFitCutoff - 1);
  vector<int> whichK = parents[Parent].whichK;
  for(int i=0; i<nToMutate; i++) {
    int whichToMutate = rng.randInt(whichK.size() - 1);
    int kToMake = getRandomK(res[whichToMutate], rng);
    whichK[whichToMutate] = kToMake;
  }
  INNERPOPSTORE result;
  result.whichK = whichK;
  result.score = -1.0f;
  return result;
}

INNERPOPSTORE crossOver(const vector<INNERPOPSTORE> &parents, int mostFitCutoff, MTRand &rng) {
  int firstParent(0), secondParent(0);
  while(firstParent == secondParent) {
    firstParent = rng.randInt(mostFitCutoff-1);
    secondParent = rng.randInt(mostFitCutoff-1);
  }
  vector<int> firstWhichK = parents[firstParent].whichK;
  vector<int> secondWhichK = parents[secondParent].whichK;
  assert(firstWhichK.size() == secondWhichK.size());

  int whereToCross = 0;
  /* The way the below is set up, a whereToCross of 0 would just copy one parent into the child */
  while(whereToCross == 0) {
    whereToCross = rng.randInt(firstWhichK.size() - 1);
  }
  vector<int> newWhichK;
  for(int i=0; i<whereToCross; i++) {
    newWhichK.push_back(firstWhichK[i]);
  }
  for(int i=whereToCross; i<secondWhichK.size(); i++) {
    newWhichK.push_back(secondWhichK[i]);
  }

  /* Check that the two parents aren't identical to the child */
  /* First parent */
  bool bad1 = true;
  for(int i=0; i<newWhichK.size(); i++) {
    if(newWhichK[i] != firstWhichK[i]) { bad1 = false; }
  }
  bool bad2 = true;
  for(int i=0; i<newWhichK.size(); i++) {
    if(newWhichK[i] != secondWhichK[i]) { bad2 = false; }
  }
  INNERPOPSTORE result;
  if(bad1 || bad2) {
    printf("Crossover failed to produce new children...\n");
    result.score = -2.0f;
  } else {  result.score = -1.0f;  }
  result.whichK = newWhichK;
  return result;
}


vector<int> randomK(const vector<GAPFILLRESULT> &res, MTRand &rng) {
  vector<int> rK;
  for(int i=0; i<res.size(); i++) {
    rK.push_back(getRandomK(res[i], rng));
  }
  return rK;
}

int getRandomK(const GAPFILLRESULT &res, MTRand &rng) {
  if(res.deadEndSolutions.size() < 1) { assert(false); }
  return rng.randInt(res.deadEndSolutions.size() - 1);
}

/* Calculate a score for the current model PROBLEM...  
 Lower scores are better.
 Total score = [NE + WKO + 0.0001 TC]/[NM + TKO + NR]

 NE --> Number of essential magic exits
 WKO --> Wrong knockout predictions
 TC --> Total cost

 NM --> Number of metabolites in model
 NR --> Number of reactions in the model
 TKO --> Total knockout predictions (genes not in the model are ignored)

 Nice properties: 
   - No knockout data --> just based on likelihood + parsimony
   - No magic exits needed --> just based on likelihood + knockouts
   - No magic exits or knockouts --> just based on the average likelihood
   (which reduces to the likelihood of gapfill rxns because everything else is the same in different models)

 Bad thing: 
 The contribution of likelihood has to be weighed to be less than the contribution of exits or KOs (hence the 0.0001 - this ensures that the total difference is always less than 1
 so that it essentially causes the equivalent solutions to be ordered by likelihood unless there is a DRASTIC difference, in which case
 we probably want to weigh the more-likely one higher)

 */
double innerScore(PROBLEM &model, const PROBLEM &problemSpace, vector<int> &essential) {
  essential = minimizeExits(model);
  int numKoCorrect = 0;
  int totalKo = 0;
  int numRxns = model.fullrxns.rxns.size();
  int numMets = model.metabolites.mets.size();
  printf("Number of essential exits for given model: %d\n", (int)essential.size());
  printf("Essential exits: \n");
  for(int i=0; i<essential.size(); i++) {
    printf("%s\t", model.fullrxns.rxnFromId(essential[i]).name);
  }
  printf("\n");

  double totalCost = 0.0f;
  for(int i=0; i<model.fullrxns.rxns.size(); i++) { totalCost += model.fullrxns.rxns[i].current_likelihood;  }

  /* Get knockout data and test lethality for each */
  double score = ( (double)essential.size() + (double)numKoCorrect + 0.0001f * totalCost )/
                 ( (double)numMets + (double) totalKo + (double) numRxns);
  return score;
}

/* Returns TRUE if the knockout was predicted to be lethal and FALSE otherwise */
bool knockoutLethality(PROBLEM &modified, const string &genename) {
  return false;
}


/* Given the model structure and the large problemSpace (containing all reactions in fullrxns),
   add the reactions in gapfillResult and any missing metabolites to the model.

   whichK: The k-value for the gapfill results to add */
void addGapfillResultToProblem(PROBLEM &model, const PROBLEM &problemSpace, const GAPFILLRESULT &gapfillResult, int whichK) {
  if( gapfillResult.deadEndSolutions.size() <= whichK ) {
    printf("ERROR: Asked for path %d for output metabolite %s but no such path exists!\n", 
				 whichK, problemSpace.metabolites.metFromId(gapfillResult.deadMetId).name);
    assert(false);
  }
  for(int i=0; i<gapfillResult.deadEndSolutions[whichK].size(); i++) {
    REACTION toAdd = problemSpace.fullrxns.rxnFromId(gapfillResult.deadEndSolutions[whichK][i]);
    model.fullrxns.addReaction(toAdd);
    for(int j=0; j<toAdd.stoich.size(); j++) {
      model.metabolites.addMetabolite(problemSpace.metabolites.metFromId(toAdd.stoich[j].met_id));
    }
  }
}

/* Minimizes magic exit usage in the model by running through them sequentially and seeing if they grow */
vector<int> minimizeExits(PROBLEM &model) {

  vector<int> requiredExits;
  vector<double> initialResult = FBA_SOLVE(model.fullrxns, model.metabolites);
  if(initialResult[model.fullrxns.idxFromId(_db.BIOMASS)] < _db.GROWTH_CUTOFF ) { 
    printf("ERROR: Failure to get growth after adding gapfill reactions\n");
    assert(false);
  }

  for(int i=0; i<model.fullrxns.rxns.size(); i++) {
    if(model.fullrxns.rxns[i].id >= _db.BLACKMAGICFACTOR && model.fullrxns.rxns[i].id < _db.BLACKMAGICFACTOR + _db.MINFACTORSPACING) {
      REACTION &exit = model.fullrxns.rxns[i];
      
      /* Exit is already off */
      if(rougheq(exit.lb, 0.0f, _db.FLUX_CUTOFF)==1 && rougheq(exit.ub, 0.0f, _db.FLUX_CUTOFF)==1) { continue; }
      
      /* Turn off exit and re-run FBA */
      double oldLb = exit.lb; double oldUb = exit.ub;
      exit.lb = 0; exit.ub = 0;
      vector<double> newResult = FBA_SOLVE(model.fullrxns, model.metabolites);
      if( newResult[model.fullrxns.idxFromId(_db.BIOMASS)] < _db.GROWTH_CUTOFF ) {     
	if(_db.DEBUGGAPFILL) { printf("Exit %s predicted essential\n", exit.name); }
	exit.lb = oldLb;
	exit.ub = oldUb;
	requiredExits.push_back(exit.id);
      } else {
	if(_db.DEBUGGAPFILL) { printf("Exit %s predicted nonessential and will be turned off!\n", exit.name); }
      }
    }
  }

  return requiredExits;
}

/* Add magic exits for everything (entrances for cofactors are allowed)
   Must be called based on the FULLRXN list... make sure you replace syns with fullrxns before calling this function!
   Returns a PROBLEM containing exactly the reactions found in pList, plus magic exits for all of them.
*/
PROBLEM setUpGapfill(const PROBLEM &problemSpace, const vector<PATHSUMMARY> &pList) {

  PROBLEM working;

  vector<int> rxnList = getAllPathRxns(pList);
  for(int i=0;i<rxnList.size();i++){rxnList[i] = abs(rxnList[i]);}

  /* In case there was no path found to a biomass metabolite, we still want to include it in the model! */
  rxnList.push_back(_db.BIOMASS);

  /* Set up a map from ID to direction to avoid duplicates */
  map<int, int> metId2Dir;
  for(int i=0; i<rxnList.size(); i++) {
    REACTION tmpRxn = problemSpace.fullrxns.rxnFromId(rxnList[i]);
    for(int j=0; j<tmpRxn.stoich.size(); j++) {
      METABOLITE tmpMet = problemSpace.metabolites.metFromId(tmpRxn.stoich[j].met_id);
      if(tmpMet.secondary_lone == 1 || (!tmpMet.secondary_pair.empty())) {
	metId2Dir[tmpMet.id] = 0;
      } else { 
	/* Only exits allowed for non-cofactor things */
	metId2Dir[tmpMet.id] = 1;
      }
    }
  }

  /* Set up making magic exits for every metabolite in the model */
  vector<int> exitIds;
  vector<int> dirs;
  for(map<int, int>::iterator it=metId2Dir.begin(); it!=metId2Dir.end(); it++) {
    exitIds.push_back(it->first);
    dirs.push_back(it->second);
  }

  REACTION biomass = problemSpace.fullrxns.rxnFromId(_db.BIOMASS);
  int basenum;
  makeSimulatableModel(pList, problemSpace, biomass, exitIds, dirs, working, basenum);

  /* Needed for cost calculation (scoring function) */
  adjustLikelihoods(working.fullrxns.rxns, 1.0f, -3.0f, 1.1f, -10.0f, true);

  return working;
}



/* Second method for gapfind/gapfill - crawl along each pathway as necessary to get growth... */
vector<int> gapFind(const PROBLEM &md, const vector<PATHSUMMARY> &psum) {

  PROBLEM model = md;

  /* Initialize list of entrances */
  /* Also, turn any EXITS off (keep ENTRANCES on) */
  vector<int> entranceSet;
  for(int i=0; i<model.fullrxns.rxns.size(); i++) {
    if(model.fullrxns.rxns[i].id >= _db.BLACKMAGICFACTOR && model.fullrxns.rxns[i].id < _db.BLACKMAGICFACTOR + _db.MINFACTORSPACING) {
      if(model.fullrxns.rxns[i].lb < 0.0f) { entranceSet.push_back(model.fullrxns.rxns[i].id); }
      model.fullrxns.rxns[i].ub = 0.0f;
    }
  }

  vector<int> deadEnds;
  REACTION bm = model.fullrxns.rxnFromId(_db.BIOMASS);
  for(int i=0; i<bm.stoich.size(); i++) {
    /* Test for existing production of the biomass component. */
    int grExit = FindExchange4Metabolite(model.fullrxns.rxns, bm.stoich[i].met_id);
    if(grExit == -1) { printf("ERROR: No exchange reaction found for metabolite %s which is impossible under our proposed schema...\n", bm.stoich[i].met_name); assert(false); }

    /* If we can already make the target without adding more magic exits, great. */
    vector<int> obj(1, grExit); vector<double> coeff(1, 1.0f);    
    GLPKDATA data(model.fullrxns, model.metabolites, obj, coeff, 1);
    vector<double> fbaResult = data.FBA_SOLVE();
    if(fbaResult[model.fullrxns.idxFromId(grExit)] > 0.0f) { continue; }
    
    /* Get the psum associated with the particular biomass component in question. */
    PATHSUMMARY bmPath;
    for(int j=0; j<psum.size(); j++) {
      if(psum[j].outputId == bm.stoich[i].met_id) { bmPath = psum[j]; break; }
    }

    /* Starting from the reaction with highest priority, see if the metabolites on one side of the equation can be produced (in the net) but the metabolites
       on the other cannot. This indicates a situation like A --> B + C - Turn on magic exit for C and see that it cannot be produced even though A can. That means
       B must be a flux bottleneck. */
    for(int j=0; j<bmPath.rxnPriority.size(); j++) {
      vector<STOICH> st = model.fullrxns.rxnFromId(bmPath.rxnPriority[j]).stoich;
      sort(st.begin(), st.end());

      bool reactantsMade(true);
      bool productsMade(true);

      /* Test for production of (nominal) reactants and products */
      for(int k=0; k<st.size(); k++) {
	int exitId = FindExchange4Metabolite(model.fullrxns.rxns, st[k].met_id);
	if(exitId == -1) { printf("ERROR: Missing exchange reaction...\n"); assert(false); }
	/* Turn exchange on and try to get flux through it */
        model.fullrxns.rxnPtrFromId(exitId)->ub = 1000.0f;
	GLPKDATA tmp(model.fullrxns, model.metabolites, obj, coeff, 1);
        vector<double> res = tmp.FBA_SOLVE();
	if(res[model.fullrxns.idxFromId(exitId)] < _db.FLUX_CUTOFF) { 
	  if(st[k].rxn_coeff < 0.0f) { reactantsMade = false; }
	  else { productsMade = false; }
	}
	model.fullrxns.rxnPtrFromId(exitId)->ub = 0.0f;
      }

      /* Test if reactans can be made but products cannot - and test if adding a exit of one chemical allows flux through another (on the same side of the reaction) 
       If it does, make that exit permanent and add it to our list */
      if( (productsMade && !reactantsMade) || ( reactantsMade && !productsMade) ) { 
	
      }
    } 
  }

  return deadEnds;
}

/* model is a PROBLEM structure containing all of the reactions in the model [see setUpGapfill]
   problemSpace is the complete problemSpace (containing all of the reactions); used for gapfilling
   It must already grow and have magic exits available for all metabolites. The easiest way to accomplish this is to call setUpGapfill 

    Uses the linprog technique to find gaps then Dijkstras to attempt to fill them.
    Note - this function only includes things within 2x cost of the most likely path...
    If there are no solutions for a particular output, no GAPFILLRESULT will exist for it. Otherwise, the solutions will be in
    the GAPFILLRESULT (see Datastructures.h for details) */
vector<GAPFILLRESULT> gapFindGapFill(PROBLEM &model, const PROBLEM &problemSpace, int gapfillK) {

  vector<GAPFILLRESULT> result;
  vector<int> obj(1, _db.BIOMASS); vector<double> coeff(1, 1.0f);

  if(_db.PRINTSHOULDGROW) { MATLAB_out("./outputs/Pre_gapfind", model.fullrxns.rxns); }

  GLPKDATA datainit(model.fullrxns, model.metabolites, obj, coeff, 1);

  vector<double> fluxVector1 = datainit.FBA_SOLVE();
  if( fluxVector1[model.fullrxns.idxFromId(_db.BIOMASS)] < _db.GROWTH_CUTOFF ) { printf("ERROR: No growth after gapfill setup\n"); assert(false); }

  GLPKDATA data(model.fullrxns, model.metabolites, obj, coeff, 1);

  vector<int> idVector;
  int status  = data.gapFindLinprog(idVector);
  set<int> idList; for(int i=0; i<idVector.size(); i++) { idList.insert(model.fullrxns.rxnFromId(idVector[i]).stoich[0].met_id); }

  /* Skip over the turning off of magic exits if gapfinding failed */
  if(status == -1) { printf("Gap pruning failed! Will attempt to recover using ALL possible exits...\n"); 
  } else {  
    for(int i=0; i<model.fullrxns.rxns.size(); i++) {
      /* Only work on magic exits */
      if(model.fullrxns.rxns[i].id < _db.BLACKMAGICFACTOR || model.fullrxns.rxns[i].id > _db.BLACKMAGICFACTOR + _db.MINFACTORSPACING) { continue; }    
      /* Turn off the exit if it's not on the list we obtained from gapFindLinprog() */
      set<int>::iterator it = idList.find(model.fullrxns.rxns[i].stoich[0].met_id);
      if(it == idList.end()) {
	model.fullrxns.change_Lb_and_Ub(model.fullrxns.rxns[i].id, 0.0f, 0.0f);
	if(_db.DEBUGGAPFILL) { printf("Turned off reaction %s\n", model.fullrxns.rxns[i].name); }
      } else {
	if(_db.DEBUGGAPFILL) { printf("Kept on reaction %s\n", model.fullrxns.rxns[i].name); }
      }
    }
  }

  if(_db.PRINTSHOULDGROW) { MATLAB_out("./outputs/Pre_gapfill", model.fullrxns.rxns); }

  /* DATA struct with the exits actually turned off */
  GLPKDATA data2(model.fullrxns, model.metabolites, obj, coeff, 1);
  vector<double> fluxVector = data2.FBA_SOLVE();
  if( fluxVector[model.fullrxns.idxFromId(_db.BIOMASS)] < _db.GROWTH_CUTOFF ) { printf("ERROR: Unable to get growth from specified set of exits in gapFindGapFill\n"); assert(false); }

  /* listOfLists[i][j] is the list of reactions composing the j'th possibly viable solution to gapfill for metabolite i */
  vector< vector< vector< int > > > listOfLists;
  PROBLEM tmpProblem = problemSpace;

  /* Come up with a list of gapfill solutions */
  for(set<int>::iterator it=idList.begin(); it != idList.end(); it++) {
    /* completeList[i] is the list of reactions composing the i'th potentially viable solution to a particular gapfill problem */
    vector< vector<int> > completeList;
    /* WE ONLY try to fill the gap with Dijkstras. The reasoning is that there could be a very likely solution with multiple reactions that is missed
       because there is a single-reaction (but quite unlikely) solution available */
    vector< vector< int> > dijkstrasSln = fillGapWithDijkstras(model.fullrxns, model.metabolites, tmpProblem, *it, 1, gapfillK);
    for(int j=0; j<dijkstrasSln.size(); j++) {
      /* FIXME: Why does the fillGapWithDijkstras sometimes give us empty results at the end of vectors with non-empty results? */
      if(dijkstrasSln[j].empty()) { continue; }
      if(_db.PRINTGAPFILLRESULTS) { printf("Dijkstras solution for exit of metabolite %s (magic exit): \n", 
					  model.metabolites.metFromId(*it).name);
	printRxnsFromIntVector(dijkstrasSln[j], problemSpace.fullrxns);
      }
      completeList.push_back(dijkstrasSln[j]);
    }

    /* Fill entrances for things that can have them */
    METABOLITE tmpMet = model.metabolites.metFromId(*it);
    if(tmpMet.secondary_lone == 1 || !tmpMet.secondary_pair.empty()) {
      dijkstrasSln = fillGapWithDijkstras(model.fullrxns, model.metabolites, tmpProblem, *it, -1, gapfillK);
      for(int j=0; j<dijkstrasSln.size(); j++) {
	if(_db.PRINTGAPFILLRESULTS) {
	  printf("Dijkstras solution for exit of metabolite %s (magic entrance): \n", tmpMet.name);
	  printRxnsFromIntVector(dijkstrasSln[j], problemSpace.fullrxns); }
	completeList.push_back(dijkstrasSln[j]); 
      }
    }

    /* If no solutions were found in either direction we don't want to return an empty GAPFILLRESULT */
    if(completeList.empty()) { continue; }

    /* Set up a gapfill result */
    GAPFILLRESULT res;
    res.deadMetId = *it;
    res.deadEndSolutions = completeList;
    result.push_back(res);
  }

  return result;
}

/* fillGapWithDijkstras: Fill a CNP (direction < 0) or PNC (direction > 0) gap by invoking Dijkstras algorithm

After running Dijkstras, it checks if the solution allowed a growth to occur. If it did, it is saved as a good solution. 
Returns a list of lists reaction IDs that fills the gap or 
empty vector if nothing fills the gap and lets it grow.
 */
vector<vector<int> > fillGapWithDijkstras(RXNSPACE &workingRxns, METSPACE &workingMets, PROBLEM &wholeProblem, int toFix, int direction, int K) {

  /* Cutoff - if the total cost of a gap becomes more than badCut * the shortest path, we just throw it out */
  double badCut = 2;

  /* Replace stoich_part with stoich */
  for(int i=0; i<wholeProblem.fullrxns.rxns.size(); i++) { wholeProblem.fullrxns.rxns[i].stoich_part = wholeProblem.fullrxns.rxns[i].stoich;  }

  vector<vector<int> > rxnsFillingGap;
  double minLength(10000.0f);
  METSPACE &allMets = wholeProblem.metabolites;  RXNSPACE &allRxns = wholeProblem.fullrxns;
  int origRxnSize = workingRxns.rxns.size();  int origMetSize = workingMets.mets.size();

  /* Find the magic exit for the given metabolite toFix, and turn it off */
  int meId = FindExchange4Metabolite(workingRxns.rxns, toFix);
  double oldLb(-1.0f), oldUb(-1.0f);
  if(meId == -1) {
     printf("WARNING: metabolite %d was passed to fillGapWithDijkstras but does not have an exchange reaction to turn off in workingRxns\n", toFix);
  } else {
    oldLb = workingRxns.rxnPtrFromId(meId)->lb;  oldUb = workingRxns.rxnPtrFromId(meId)->ub;
    workingRxns.rxnPtrFromId(meId)->lb = 0.0f;  workingRxns.rxnPtrFromId(meId)->ub = 0.0f;
  }

  /* Don't use reactions that are already in the network, they clearly don't work. 
   Change this in allRxns, not in workingRxns, becuase Dijktras utilizes allRxns and not workingRxns to try to fill the gap */
  for(int i=0; i<workingRxns.rxns.size(); i++) {
    // Magic exits aren't in allRxns but are found in workingRxns
    if(allRxns.idIn(workingRxns[i].id)) {
      REACTION *rxnPtr = allRxns.rxnPtrFromId(workingRxns[i].id);
      rxnPtr->old_likelihood = rxnPtr->current_likelihood;
      rxnPtr->current_likelihood = -1;
    }
  }

  /* direction > 0 means it's a PNC gap so we need to reverse available reactions to try to fill it */
  if(direction > 0) { ReverseReversible(allRxns.rxns);  }

  /* Treat everything in our working network as an "input" metabolites EXCEPT the dead end we want to fill */
  METSPACE inputs;
  for(int i=0; i<workingMets.mets.size(); i++) {
    if(workingMets.mets[i].id!=toFix) {
      inputs.addMetabolite(workingMets.mets[i]);
    }
  }

  /* Also include actual inputs[media] as inputs whether they are in our actual network or not.
     FIXME - should be dependent on the growth condition used in the PATH - this one allows it to use any of the GROWTH conditions */
  for(int i=0; i<wholeProblem.growth.size(); i++) {
    for(int j=0; j<wholeProblem.growth[i].media.size(); j++) {
      if(wholeProblem.growth[i].media[j].id!=toFix) {
	inputs.addMetabolite(allMets.metFromId(wholeProblem.growth[i].media[j].id));
      }
    }
  }

  /* The target is the output metabolite */
  METABOLITE output = allMets.metFromId(toFix);
  vector<PATH> result;
  kShortest(result, allRxns, allMets, inputs, output, K);

  /* Change the reversibilities back to the way they were (for PNC) */
  if(direction > 0) {  ReverseReversible(allRxns.rxns);  }

  double bestLikelihood(-1.0f);

  /* Check that the gapfill reactions can carry flux */
  /* TODO - need to make sure current_likelihood is filled in correctly here */
  for(int i=0; i<result.size(); i++) {  
    if(result[i].rxnIds.empty()) { break;  }    
    /* Add the gapfilling reactions and metabolites to our network (for FBA test) */
    vector<int> currentSolution = result[i].rxnIds;
    for(int j=0; j<currentSolution.size(); j++) {
      
      workingRxns.addReaction(allRxns.rxnFromId(currentSolution[j]));
      for(int k=0; k<workingRxns.rxns.back().stoich.size(); k++) {
	int metId = workingRxns.rxns.back().stoich[k].met_id;
	workingMets.addMetabolite(allMets.metFromId(metId));
      }
    }
    
    /* Check that gapfill reactions carry flux - add them if they do this and they satisfy
       the likelihood cutoff compared to the best solution */
    vector<int> obj(1, result[i].rxnIds[0]);  vector<double> coeff(1, 1.0f);
    GLPKDATA data(workingRxns, workingMets, obj, coeff, 1);
    vector<double> fbaResult = data.FBA_SOLVE();
    if( rougheq(fbaResult[workingRxns.idxFromId(result[i].rxnIds[0])], 0.0f, _db.FLUX_CUTOFF) == 0 ) {
      /* Apply cost cutoff */
      double totalCost = 0.0f;
      for(int k=0; k<currentSolution.size(); k++) { 
	totalCost += allRxns.rxnFromId(currentSolution[k]).current_likelihood;
      }
      /* First one should always be accepted */
      if(rxnsFillingGap.empty()) { 
	bestLikelihood = totalCost;
	rxnsFillingGap.push_back(currentSolution);
      } else { /* Apply cutoff */
	if(totalCost < badCut * bestLikelihood) {
	  rxnsFillingGap.push_back(currentSolution);
	}
      }
    }

    /* Remove things that were newly added */
    while(workingRxns.rxns.size() > origRxnSize) {  workingRxns.removeRxnFromBack();  }
    while(workingMets.mets.size() > origMetSize) {  workingMets.removeMetFromBack();  } 
  }

  /* Change the likelihoods back to what they were before (no longer -1) */
  for(int i=0; i<workingRxns.rxns.size(); i++) {
    if(allRxns.idIn(workingRxns[i].id)) {
      REACTION *rxnPtr = allRxns.rxnPtrFromId(workingRxns[i].id);
      rxnPtr->current_likelihood = rxnPtr->old_likelihood;
    }
  }

  /* Turn the magic exit we were filling back on */
  if(meId != -1) {
    workingRxns.rxnPtrFromId(meId)->lb = oldLb;
    workingRxns.rxnPtrFromId(meId)->ub = oldUb;
  }

  return rxnsFillingGap;
}

void setSpecificGrowthConditions(PROBLEM &model, const GROWTH &growth) {
  /* Turn off exchange reactions (aside from those that are magic exits / entrances) */
  for(int i=0; i<model.fullrxns.rxns.size(); i++) {
    if(model.fullrxns.rxns[i].id >= _db.BLACKMAGICFACTOR && model.fullrxns.rxns[i].id < _db.BLACKMAGICFACTOR + _db.MINFACTORSPACING) {  continue;    }
    if(model.fullrxns.rxns[i].isExchange) {
      model.fullrxns.rxns[i].lb = 0;
      model.fullrxns.rxns[i].net_reversible = 1;
    }
  }

  /* Set the media components to be uptaken at specified rates (negative because they are uptaken) */
  map <int,bool> meList;
  for(int i=0; i<growth.media.size(); i++) {
    int meId = FindExchange4Metabolite(model.fullrxns.rxns, growth.media[i].id);
    if(meId == 0) { 
      printf("ERROR: No exchange present for media condition %s after calling checkExchangesAndTransports \n", 
	     model.metabolites.metFromId(growth.media[i].id).name);
      assert(false);
    }
    model.fullrxns.rxnPtrFromId(meId) -> lb =  -growth.media[i].rate;
    //model.fullrxns.rxnPtrFromId(meId) -> lb = -1000.0f;
    model.fullrxns.rxnPtrFromId(meId) -> ub = 1000.0f;
    model.fullrxns.rxns[i].net_reversible = 0;
    meList[meId] = "true";
  }

  /* Set the byproducts to be output at non-specified rates */
  for(int i=0; i<growth.byproduct.size(); i++) {
    int meId = FindExchange4Metabolite(model.fullrxns.rxns, growth.byproduct[i].id);
    if(meId == 0) { 
      printf("ERROR: No exchange present for byproduct %s after calling checkExchangesAndTransports \n", 
	     model.metabolites.metFromId(growth.byproduct[i].id).name);
      assert(false);
    }
    model.fullrxns.rxnPtrFromId(meId) -> ub =  1000.0f;
    if(meList.count(meId)<1){
      model.fullrxns.rxnPtrFromId(meId) -> lb = 0.0f;
      model.fullrxns.rxns[i].net_reversible = 1;
    }
    else{
      printf("WARNING: %s both a byproduct and a media component.\n",model.metabolites.metFromId(growth.byproduct[i].id).name);
    }
  }
}

void addATPM(PROBLEM &A, ANSWER &B){
  int atpmId = Name2Ids(A.fullrxns.rxns, _db.ATPM_name);
  B.reactions.addReaction(A.fullrxns.rxnFromId(atpmId));
}
