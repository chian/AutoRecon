#include "Exchanges.h"
#include "genericLinprog.h"
#include "Grow.h"
#include "kShortest.h"
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
#include <vector>

/* Utility function identifies if we should include a reaction when unsynonymmming the PATHSUMMARY vector 
dir: <0 if we want backward directions and >0 if we want forward directions

Tested reversibilities on 3-31-11 - seems to be working as expected and yielding the correct ID's */
vector<int> findRxnsToIncludeFromSyn(const REACTION &TMPSYN, const RXNSPACE &fullrxns, int dir) {
  vector<int> result;

  /* Magic bridges are filled from the calling function in a separate step */
  if( TMPSYN.id >= _db.MAGICBRIDGEFACTOR && TMPSYN.id <= _db.MAGICBRIDGEFACTOR + _db.MINFACTORSPACING) {     
    return result;  }
  
  /* Un-Synonymize - ensure correct direction [or reversible] */
  if(TMPSYN.id >= _db.SYNFACTOR && TMPSYN.id <= _db.SYNFACTOR + _db.MINFACTORSPACING ){ 
    for(int j=0;j<TMPSYN.syn.size();j++) {

      /* Find one metabolite in common between the synrxn and its full version */
      REACTION TMPRXN = fullrxns.rxnFromId(TMPSYN.syn[j]);
      int metIdx;
      for(int k=0;k<TMPRXN.stoich_part.size();k++) {
	if(TMPRXN.stoich_part[k].met_id == TMPSYN.stoich_part[0].met_id) {
	  metIdx = k;  break;  }
      }
      
      /* Identify if the reactions go the same way */
      bool sameDir = false;
      if( ( TMPRXN.net_reversible == 1 && TMPRXN.stoich_part[metIdx].rxn_coeff*TMPSYN.stoich_part[0].rxn_coeff > 0 ) ||
	  ( TMPRXN.net_reversible == -1 && TMPRXN.stoich_part[metIdx].rxn_coeff*TMPSYN.stoich_part[0].rxn_coeff < 0) ) {
	/* Both reactions go the same way */
	sameDir = true;
      }
      /* If dir < 0 we only want things that go the OPPOSITE way (or that are reversible */
      if( TMPRXN.net_reversible == 0 || (sameDir && dir > 0 ) || (!sameDir && dir < 0)) {
	result.push_back(TMPRXN.id);
      }
    }
    return result;
  }
  
  /* Include just the reaction itself if it has no synonyms */
  return TMPSYN.syn;
}

/* Returns a list of reactions that are suggested to fill all of the magic bridges in
  the vector<PATHSUMMARY> psum, based on all of the synrxns in ProblemSpace.

  pList should be based on SYNRXNS (this is intended to be wrapped into the function that converts SYNRXNS
  into FULLRXNS)

  Returns a vector of reaction IDs (in no particular order) that are connected the correct way to
  the rest of the PATHSUMMARY, based on which way we need the bridge to go.

  Any non-secondary metabolites in the reaction must be present somewhere in the pList

  I'm HOPING that this will make our job gapfilling easier. Sigh... yeah, right. */
vector<int> fillMagicBridges(const PATHSUMMARY &currentPsum, const vector<PATHSUMMARY> &pList, 
			     const PROBLEM &ProblemSpace, int bridgeToFill, vector<int> &filledBridges) {

  // This is needed so I can keep track of the recursion and unwrap it later so I know exactly what magic bridges
  // were filled up
  filledBridges.push_back(abs(bridgeToFill));

  vector<int> suggestedRxnIds;

  /* Get a list of metabolites from reactions in the pList (needed later...) */
  vector<int> allRxnIds = getAllPathRxns(pList);
  for(int i=0; i<allRxnIds.size(); i++) { allRxnIds[i] = abs(allRxnIds[i]); }

  RXNSPACE modelSynRxns;  pullOutRxnsbyIds(ProblemSpace, allRxnIds, modelSynRxns);
  METSPACE modelSynMets;  pullOutMets(modelSynRxns, ProblemSpace.metabolites, modelSynMets);
  /* Temporary list of all metabolites - needed because we need to keep track of links between metaboiltes
     and reactions... */
  METSPACE allMets = ProblemSpace.metabolites;
  calcMetRxnRelations(ProblemSpace.fullrxns, allMets);

  REACTION bridge = ProblemSpace.synrxns.rxnFromId(abs(bridgeToFill));
  
  /* We want the filling reaction to contain [metToConsume --> metToProduce] */
  int metToProduce; int metToConsume;
  for(int i=0; i<bridge.stoich_part.size(); i++) {
    if(bridgeToFill < 0) { 
      if(bridge.stoich_part[i].rxn_coeff < 0) { metToProduce = bridge.stoich_part[i].met_id;   }
      else                                    { metToConsume = bridge.stoich_part[i].met_id;   }
    } else {
      if(bridge.stoich_part[i].rxn_coeff < 0) { metToConsume = bridge.stoich_part[i].met_id;   }
      else                                    { metToProduce = bridge.stoich_part[i].met_id;   }
    }
  }

  vector<int> rxnsWithProdMet = allMets.metFromId(metToConsume).rxnsInvolved_nosec;

  vector<int> tmpSuggested;
  vector<bool> allCof;
  for(int i=0; i<rxnsWithProdMet.size(); i++) {
    REACTION tmp = ProblemSpace.fullrxns.rxnFromId(rxnsWithProdMet[i]);

    /* don't try to fill with the biomass equation. That's bad. */
    if(tmp.id == _db.BIOMASS) { continue; }

    /* First test if the metToProduce and metToConsume are present at all */
    bool producePresent(false); bool consumePresent(false);
    double produceCoeff(0.0f); double consumeCoeff(0.0f);
    for(int j=0; j<tmp.stoich.size(); j++) {
      if(tmp.stoich[j].met_id == metToProduce) {
	producePresent = true;
	produceCoeff = tmp.stoich[j].rxn_coeff;
      }
      if(tmp.stoich[j].met_id == metToConsume) {
	consumePresent = true;
	consumeCoeff = tmp.stoich[j].rxn_coeff;
      }
    }
    if( !(producePresent && consumePresent) ) {
      continue;
    }
    /* Make sure the two metabolites are on OPPOSITE sides of the equation */
    if(produceCoeff * consumeCoeff > 0.0f) { continue; }

    /* Now test if they're the correct direction */
    bool dirGood = true;
    for(int j=0; j < tmp.stoich.size(); j++) {
      if(tmp.stoich[j].met_id == metToProduce) { 
	if(tmp.net_reversible != 0) {
	  if(tmp.net_reversible * tmp.stoich[j].rxn_coeff < 0 ) { dirGood = false; break; }
	}
      }
      if(tmp.stoich[j].met_id == metToConsume) {
	if(tmp.net_reversible != 0) {
	  if(tmp.net_reversible * tmp.stoich[j].rxn_coeff > 0) { dirGood = false; break; }
	}
      }
    }

    if(!dirGood) { continue; }

    /* Look for all of the metabolites in the new reaction and make sure that any non-secondary metabolites
       are present in the path list.
       [FIXME - I dont check that the direction is correct, will I have to?] */

    bool metsIn = true;
    for(int j=0; j<tmp.stoich.size(); j++) {
      METABOLITE tmpMet = ProblemSpace.metabolites.metFromId(tmp.stoich[j].met_id);
      if(tmpMet.secondary_lone == 1 | (!tmpMet.secondary_pair.empty())) { continue; }
      if(!modelSynMets.idIn(tmpMet.id)) { metsIn = false; break; }
    }

    if(!metsIn) { continue; }

    /* Done with all the tests, so just return the ID */
   tmpSuggested.push_back(tmp.id);

    /* ALmost always if the magic bridge is fixable by only one reaction, 
       the correct answer(s) is / are reactions that only
       contain cofactors. I looked at all the ones that we run into from DIjkstras
       and the only exception was accoa/coa, which was filled by PDH. */
    bool oneAllCof = true;
    for(int j=0; j<tmp.stoich.size(); j++) {
      METABOLITE tmpMet = ProblemSpace.metabolites.metFromId(tmp.stoich[j].met_id);
      if(tmpMet.secondary_lone == 1 | (!tmpMet.secondary_pair.empty())) { continue; }
      oneAllCof = false;
      break;
    }
    allCof.push_back(oneAllCof);
  } 

  /* If any of these reactions are all secondaries include them at the expense of the others */
  bool anyCof = anyVec(allCof);
  if(anyCof) {
    for(int i=0; i<allCof.size(); i++) {
      if(allCof[i]) { suggestedRxnIds.push_back(tmpSuggested[i]); }
    }
  } else {  suggestedRxnIds = tmpSuggested;  }

  if(suggestedRxnIds.empty()) {
    vector<int> inputs = Load_Inputs_From_Growth(ProblemSpace.growth[currentPsum.growthIdx[0]]);
    /* This is the difference from the paths we had before... 
     We are now adding the metabolite we're trying to consume to our list */
    inputs.push_back(metToConsume);
    RXNSPACE allrxns = ProblemSpace.synrxns;
    METSPACE allmets = ProblemSpace.metabolites;

    METSPACE inputSpace(allmets, inputs);
    METABOLITE output = allmets.metFromId(metToProduce);
    int K = 5;
    vector<PATH> result;

    kShortest(result, allrxns, allmets, inputSpace, output, K);

    /* Figure out what is connected to the rest of the network based on the DIjkstras spitouts */
    vector<int> gapfillRxnDirIds;
    for(int i=0; i<result.size(); i++) {
      /* Throw out the trivial solution that contains only the magic exit we're trying to fill
       I kept this solution around so I Can check my sanity that I actually get such a trivial solution
       after calling Dijkstras, and to avoid having to removeRxn */
      if(result[i].rxnIds.size() == 1 && result[i].rxnIds[0] >= _db.MAGICBRIDGEFACTOR &&
	 result[i].rxnIds[0] < _db.MAGICBRIDGEFACTOR + _db.MINFACTORSPACING) { continue; }

      /* Test inputs - all must be in the network (may not be a necessary check) */
      bool inputsIn = true;
      for(int j=0; j<result[i].inputIds.size(); j++) { 
	if(!modelSynMets.idIn(result[i].inputIds[j])) { inputsIn = false; break; }
      }
      if(!inputsIn) { continue; }

      /* Test dead ends (do the extra outputs lead anywhere?) - more important check */
      bool deadIn = true;
      for(int j=0; j<result[i].deadEndIds.size(); j++) {
	if(!modelSynMets.idIn(result[i].deadEndIds[j])) { deadIn = false; break; }
      }
      if(!deadIn) { continue; }

      /* All consistency checks passed so tag the rxnDirIDs for further analysis */
      vector<int> rxnDirIds;
      for(int j=0; j<result[i].rxnIds.size(); j++) {
	rxnDirIds.push_back(result[i].rxnIds[j] * result[i].rxnDirection[j]);
      }
      catVectors1(gapfillRxnDirIds, rxnDirIds);
      
      /* FIXME: Take a look at the other solutions and decide if we want to do this or not.
	 For now.. only take the first (most likely) solution as the right one. */
      break;

    }

    /* This path came from a list of synrxns, and therefore it is very likely that there will be
       synrxns (8000000-series) and/or magic bridges (4000000-series) to decipher here.
       Magic bridges need to be deciphered with a recursive call to this function.
       Synonym reactions can be dealt with using the findRxnsToIncludeFromSyn function. */
    for(int i=0; i<gapfillRxnDirIds.size(); i++) {
      REACTION TMPSYN = ProblemSpace.synrxns.rxnFromId(abs(gapfillRxnDirIds[i]));
      /* This should never happen, but because I'm paranoid... 
       ACtually I should be testing if there's EVER a repeat in the recursion and saving all the
      magic bridges that are ever fixed. If I get into trouble this is probably why. */
      if(abs(gapfillRxnDirIds[i]) == bridge.id) { continue; }

      /* Normal reaction */
      if(TMPSYN.id < _db.MINFACTORSPACING) { suggestedRxnIds.push_back(TMPSYN.id); continue; }
      /* Synonym */
      if(TMPSYN.id >= _db.SYNFACTOR && TMPSYN.id < _db.SYNFACTOR + _db.MINFACTORSPACING) {
	vector<int> newIds = findRxnsToIncludeFromSyn(TMPSYN,ProblemSpace.fullrxns,gapfillRxnDirIds[i]);
	catVectors1(suggestedRxnIds, newIds);
	continue;
      }
      /* Recurse */
      if(TMPSYN.id >= _db.MAGICBRIDGEFACTOR && TMPSYN.id < _db.MAGICBRIDGEFACTOR + _db.MINFACTORSPACING) {
	vector<int> newIds = fillMagicBridges(currentPsum, pList, ProblemSpace, gapfillRxnDirIds[i], filledBridges);
	catVectors1(suggestedRxnIds, newIds);
      }
    }
  }

  if(_db.DEBUGBRIDGES) {
    printf("Suggested reaction(s) to fill magic exit %s:\n", bridge.name);
    for(int i=0; i<suggestedRxnIds.size(); i++) {
      printf("%s(%1.3f)\t", ProblemSpace.fullrxns.rxnFromId(suggestedRxnIds[i]).name, ProblemSpace.fullrxns.rxnFromId(suggestedRxnIds[i]).init_likelihood);
    }
    printf("\n");
  }

  return suggestedRxnIds;
}

/* Take the synrxn ID's in a given PATHSUMMARY and replace them with the names from fullreaction 
(according to synrxns.syn) throughout the PATHSUMMARY psum and return a new PATHSUMMARY
with the ID's fixed. Only take things that go the correct direction.

allPsum: The vector of ALL the synrxns pathsummaries (needed to replace magic bridges reasonably sensibly)

CommonRxnIds tested on 3-11-11 and seems to be replaced correctly at least for the cast I tried */

PATHSUMMARY replaceWithRealRxnIds(const PATHSUMMARY &psum, const vector<PATHSUMMARY> &allPsum, 
				  PROBLEM &ProblemSpace) {
  PATHSUMMARY result = psum;
  vector<int> rxnIds = psum.rxnDirIds;
  vector<bool> toErase(ProblemSpace.synrxns.rxns.size(), false);

  int origSize = rxnIds.size();

  /* The different reactions chosen to replace a synrxn are all given adjacent priorities (it really doesn't
     matter which one is picked first because we're traversing the whole list anyway) 

     Since we're traversing rxnIds and not rxnPriority, the resulting list will be backwards... */
  vector<int> backwardsPriority;

  for(int i=0;i<origSize;i++) {
    REACTION TMPSYN = ProblemSpace.synrxns.rxnFromId(abs(rxnIds[i]));
    /* Don't replace if theres no synonyms - but DO replace rxnIds[i] with its 
       absolute value so we actually pull the right reaction out */
    if(TMPSYN.id < _db.MINFACTORSPACING) { 
      rxnIds[i] = abs(rxnIds[i]);
      backwardsPriority.push_back(rxnIds[i]);
      continue; 
    }

    if(TMPSYN.id >= _db.MAGICBRIDGEFACTOR && TMPSYN.id < _db.MAGICBRIDGEFACTOR + _db.MINFACTORSPACING) {
      /* Pass direction and ID to the sub-function so I can call 
	 findRxnsToIncludeFromSyn from there and match things up correctly */
      vector<int> filledBridges;
      vector<int> toInclude = fillMagicBridges(psum, allPsum, ProblemSpace, rxnIds[i], filledBridges);

      // This adds toInclude to the END of the list and messes with the order.
      catVectors1(rxnIds, toInclude);
      for(int j=0; j<toInclude.size(); j++) {
	toErase.push_back(false);
      }
      toErase[i] = true;

      // In our list of priorities we preserve order.
      for(int j=0; j<toInclude.size(); j++) { backwardsPriority.push_back(toInclude[j]); }
    }

    if((TMPSYN.id >= _db.SYNFACTOR && TMPSYN.id < _db.SYNFACTOR + _db.MINFACTORSPACING) ||
       (TMPSYN.id >= _db.SYNFACTOR + _db.REVFACTOR && TMPSYN.id < _db.SYNFACTOR + _db.REVFACTOR + _db.MINFACTORSPACING) ) {
      /* Get the relevant synonyms and tack them on to the end. Set up the old name 
	 (synonym name) to be deleted and set the real IDs not to be deleted */
      vector<int> toInclude = findRxnsToIncludeFromSyn(TMPSYN,ProblemSpace.fullrxns,rxnIds[i]);
      catVectors1(rxnIds, toInclude);
      
      toErase[i] = true;
      for(int j=0;j<toInclude.size();j++) {
	toErase.push_back(false);
      }

      for(int j=0; j<toInclude.size(); j++) { backwardsPriority.push_back(toInclude[j]); }
    }
  }

  /* Take the ID's that aren't set to be deleted and add them to the results */
  result.rxnDirIds.clear();

  for(int i=0;i<rxnIds.size();i++) {
    if(!toErase[i]) {
      result.rxnDirIds.push_back(rxnIds[i]);
    }
  }

  reverse(backwardsPriority.begin(), backwardsPriority.end());
  result.rxnPriority = backwardsPriority; 
  custom_unique(result.rxnDirIds);

  return result;
}

/* Add appropriately needed rxns from synrxnsR to the ProblemSpace                                                                                                                                             
   FIXME: THese should be moved into the paths2model.h / paths2model.cc files */
void addR(PATHSUMMARY &psum, PROBLEM &ProblemSpace) {
  vector<int> rxnIds = psum.rxnDirIds;
  int origSize = rxnIds.size();
  for(int i=0;i<origSize;i++){
    if(ProblemSpace.synrxns.idIn(abs(rxnIds[i]))){continue;}
    if(abs(rxnIds[i]) >= _db.REVFACTOR && abs(rxnIds[i]) < _db.REVFACTOR + _db.MINFACTORSPACING){
      REACTION TEMP = ProblemSpace.synrxnsR.rxnFromId(abs(rxnIds[i]));
      ProblemSpace.synrxns.addReaction(TEMP);
      ProblemSpace.fullrxns.addReaction(TEMP);
    }
    if(abs(rxnIds[i]) >= _db.REVFACTOR + _db.SYNFACTOR
       && abs(rxnIds[i]) < _db.REVFACTOR + _db.SYNFACTOR + _db.MINFACTORSPACING){
      REACTION TEMP = ProblemSpace.synrxnsR.rxnFromId(abs(rxnIds[i]));
      ProblemSpace.synrxns.addReaction(TEMP);
      for(int j=0;j<TEMP.syn.size();j++){
        REACTION TEMP2 = ProblemSpace.fullrxns.rxnFromId(TEMP.syn[j]);
        TEMP2.id += _db.REVFACTOR;
        TEMP2.net_reversible *= -1;
        ProblemSpace.fullrxns.addReaction(TEMP2);
      }
    }
  }
}

void addR(vector<vector<vector<PATHSUMMARY> > > &psum, PROBLEM &ProblemSpace) {
  for(int i=0;i<psum.size();i++){
    for(int j=0;j<psum[i].size();j++){
      for(int k=0;k<psum[i][j].size();k++){
        addR(psum[i][j][k],ProblemSpace);
      }
    }
  }
}

/* Perform setup steps that I always seem to get wrong...

   1: Add biomass equation to the rxnspace
   2: Check for exhanges and transporters for anything in any of the growth conditions specified by pList
   3: Add magic exits (leave blank if none) identified by metabolite ID with the specified direction (-1 or +1)
   4: Synchronize metabolite and reaction lists 

   baseNum is the number of reactions present in the RXNSPACE before adding magic (which is added LAST, 
   so we can just iterate until the end) */

void makeSimulatableModel(const vector<PATHSUMMARY> &pList, const PROBLEM &ProblemSpace, const REACTION &biomass, 
			  const vector<int> &magicIds, const vector<int> &magicDirs, PROBLEM &working, int &baseNum) {

  working.clear();
  const METSPACE &BaseMets = ProblemSpace.metabolites;
  const RXNSPACE &BaseRxns = ProblemSpace.fullrxns;
  const vector<GROWTH> &growth = ProblemSpace.growth;

  METSPACE &workingMets = working.metabolites;
  RXNSPACE &workingRxns = working.fullrxns;

  vector<int> rxnIdList = getAllPathRxns(pList);
  for(int i=0;i<rxnIdList.size();i++){rxnIdList[i] = abs(rxnIdList[i]);}
  pullOutRxnsbyIds(BaseRxns, rxnIdList, workingRxns);
  /* Add the biomass equation to BaseRxns - we never added it previously, because each path is typically
     tied to one particular output */
  workingRxns.addReaction(biomass);

  /* Checking exchanges and transporters must be done BEFORE calculating baseNum, otherwise
     we make fluxes = 0 for things that are not magic exits... */
  vector<int> growthIdx;
  for(int i=0; i<pList.size(); i++) {
    growthIdx.push_back(pList[i].growthIdx[0]);
  }
  custom_unique(growthIdx);

  /* I use a map here because we are interested in making sure that we don't duplicate metabolite ID's with different directions.
     When in doubt, make it 0 */
  map<int, int> id2Dir;

  for(int i=0; i<growthIdx.size(); i++) {
    int currentIdx = growthIdx[i];
    for(int j=0;j<growth[currentIdx].byproduct.size();j++) {
      id2Dir[growth[currentIdx].byproduct[j].id] = 1;
    }
    for(int j=0;j<growth[currentIdx].media.size();j++) {
      /* Note - because we do this second, if that ID is already in the map its reversibility 
	 gets overwritten with 0. This is what we want... if something is a media AND a byproduct
	 we be conservative and allow it to be either... */
      id2Dir[growth[currentIdx].media[j].id] = 0;
    }
  }

  vector<int> metIds,  byDirs;
  for(map<int, int>::iterator it=id2Dir.begin(); it!=id2Dir.end(); it++) {
    metIds.push_back(it->first);
    byDirs.push_back(it->second);
  }

  checkExchangesAndTransporters(working, ProblemSpace, metIds, byDirs);

  /* Deal with duplicates in the magic IDs/ magic Dirs */
  id2Dir.clear();
  for(int i=0; i<magicIds.size(); i++) {
    if(id2Dir.find(magicIds[i]) == id2Dir.end() ) {
      /* magicDirs can be direction = -1 (magic entrance), or 1 (magic exit), or 0 (reversible) */
      if(magicDirs[i] < 0) {
	id2Dir[magicIds[i]] = -1;
      } else if(magicDirs[i] == 0) {
	id2Dir[magicIds[i]] = 0;
      } else  {
	id2Dir[magicIds[i]] = 1;
      }
    } else {
      /* Make reversible if it goes in both directions */
      if(id2Dir[magicIds[i]]*magicDirs[i] < 0) {
	id2Dir[magicIds[i]] = 0;
      }
    }
  }

  vector<int> newMagic = magicIds;
  vector<int> newDirs;
  /* Decide the direction of the magic stuff */
  custom_unique(newMagic);
  for(int i=0; i<newMagic.size(); i++) {
    if(id2Dir.find(newMagic[i])==id2Dir.end()) {
      printf("Badness - for some reason magicIds[i] was not found in id2Dir\n");
      throw;
    }
    newDirs.push_back(id2Dir[newMagic[i]]);
  }

  baseNum = workingRxns.rxns.size();

  /* Add magic exchanges */
  RXNSPACE excRxnList;
  GetExchangeReactions(newMagic, newDirs, ProblemSpace, excRxnList);
  workingRxns.addReactionVector(excRxnList.rxns);

  pullOutMets(workingRxns, BaseMets, workingMets);

}

/* Utility function that:
1: Checks to see if all of the metabolites with ID's in metsToCheck exist in workingMetSet
2: If NO: Add them from allMets (it is an error if not in allMets either - should never happen since
the media is checked against the possible metabolites on the front end)
3: Checks to see if the metabolites in metsToCheck have exchange reactions in workingRxnSet
4: If not: Looks for exchange rxn in allRxns, if present, adds it to working reactions, and otherwise, 
creates one for you 
5: Does 3 and 4 for transport reactions 

This should elminiarte the need for the warning in the transport reaction adder but I may change it into an error later

If you call this with "fullrxns" or "synrxns" as both arguments it will just add the magic stuff to the reactions
that is necessary - useful for updating the fullreactions. However, I'm not sure that would work with a pass by reference... 
you may need to make an extra copy in the calling function. */
void checkExchangesAndTransporters(PROBLEM &working, const PROBLEM &ProblemSpace, const vector<int> &metsToCheck, 
				   const vector<int> &dirsToCheck) {

  RXNSPACE &workingRxns = working.fullrxns;
  METSPACE &workingMets = working.metabolites;

  const METSPACE &allMets = ProblemSpace.metabolites;

  RXNSPACE excRxns;
  RXNSPACE transporters;

  // DO NOT UNIQUE METSTOCHECK BECAUSE IT SORTS THE VALUES
  //custom_unique(metsToCheck);

  // Existence check for metabolites in metsToCheck - add things that are missing to the metabolite list 
  for(int i=0; i<metsToCheck.size();i++) {
    if(!workingMets.idIn(metsToCheck[i])) {
      if(!allMets.idIn(metsToCheck[i])) { printf("ERROR: Metaboilte %d not found in allMets although it was in metsToCheck in checkReactionsAndTransporters.\n", metsToCheck[i]);	throw;   }
      workingMets.addMetabolite(allMets.metFromId(metsToCheck[i]));
    }

    // Also add the pair (in case for example a metabolite that wasn't needed was included in the media 
    // Could cause issues if someone passes both a metabolite and its pair as mets to check so we check for that 
    // condition here  
    int pairId = inOutPair(metsToCheck[i], allMets);
    if(pairId == -1) { printf("ERROR: Metabolite %s (%d) has no internal metabolite\n", allMets.metFromId(metsToCheck[i]).name, metsToCheck[i]); throw; }
    if(!workingMets.idIn(pairId)) {
      workingMets.addMetabolite(allMets.metFromId(pairId));
    }
  }

  // Look for exchange reactions in the WHOLE LIST and add any to workingRxnSet that aren't already there.
  // If we need to, we'll add magic too. 
  GetExchangeReactions(metsToCheck, dirsToCheck, ProblemSpace, excRxns);
  for(int i=0;i<excRxns.rxns.size();i++) {
    int excId = excRxns.rxns[i].id;
    if(!workingRxns.idIn(excId)) {
      workingRxns.addReaction(excRxns.rxns[i]);
    } else {
      // Keep the reaction the same but check the reversibility and modify as needed 
      workingRxns.rxnPtrFromId(excId)->net_reversible = excRxns.rxns[i].net_reversible;
    }
  }

  // Look for transport reactions in the WHOLE LIST and add any to workingRxnSet that aren't already there. 
  // If we need to we'll add magic too EXCEPT for H. Also, add cofactors to transporters to metabolite list if necessary... 

  GetTransportReactions(metsToCheck, dirsToCheck, ProblemSpace, transporters.rxns);
  for(int i=0; i<transporters.rxns.size();i++) {
    int transId = transporters[i].id;
    if(!workingRxns.idIn(transId)) {
      workingRxns.addReaction(transporters[i]);
    } else {
      // Keep the reaction the same but check the reversibility and modify as needed
      workingRxns.rxnPtrFromId(transId)->net_reversible = transporters[i].net_reversible;
    }
    // Check for new metabolites added because of transporters and add them to the metabolite list
    for(int j=0;j<transporters[i].stoich.size();j++) {
      if(!workingMets.idIn(transporters[i].stoich[j].met_id)) {
	workingMets.addMetabolite(allMets.metFromId(transporters[i].stoich[j].met_id));
      }
    }
  }
  return;
}
