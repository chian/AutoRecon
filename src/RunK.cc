#include "DataStructures.h"
#include "kShortest.h"
#include "genericLinprog.h"
#include "MyConstants.h"
#include "Paths2Model.h"
#include "pathUtils.h"
#include "Printers.h"
#include "RunK.h"
#include "shortestPath.h"
#include "visual01.h"
#include "XML_loader.h"

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <map>
#include <omp.h>
#include <set>
#include <vector>

using std::vector;
using std::map;
using std::set;

/* Does all input parsing and hum-drum vector filling */
void InputSetup(int argc, char *argv[], PROBLEM &ProblemSpace) {

  if(_db.DEBUGSYN){printf("entering InputSetup\n");}
  if(argc!=4){
    printf("Usage: ./a.out num_threads likelihoods.xml input.xml\n");
    exit(0);}

  omp_set_num_threads(atoi(argv[1]));
  /* likelihoods.xml */
  char* docName = argv[2];
  /* input.xml */
  char* docName2 = argv[3];

  /* Parse XML files and do some of the initial setup steps (which should be moved here) */
  parseALL(docName, docName2, ProblemSpace);
  ProblemSpace.fullrxns.rxnMap();

  /* Ensure that exchangers and transprots exist for everything in the GROWTH conditions (media and byproducts) */
  vector<GROWTH> &growth = ProblemSpace.growth;
  vector<int> dirs;
  vector<int> metsToCheck;
  for(int i=0;i<growth.size();i++) {
    for(int j=0;j<growth[i].media.size();j++) {
      dirs.push_back(0); /* Note we want to allow media to be imported OR exported... (especially important for H+ because transporters depend on it) */
      metsToCheck.push_back(growth[i].media[j].id);
    }
    for(int j=0;j<growth[i].byproduct.size();j++) {
      dirs.push_back(1); /* ... but anything NOT in the media must be only exported */
      metsToCheck.push_back(growth[i].byproduct[j].id);
    }
  }

  PROBLEM TMPproblem(ProblemSpace.fullrxns, ProblemSpace);
  checkExchangesAndTransporters(ProblemSpace, TMPproblem, metsToCheck, dirs);
  TMPproblem.clear();

  MakeSynList(ProblemSpace);
  makeMagicBridges(ProblemSpace);
  ProblemSpace.synrxns.rxnMap();

  if(_db.DEBUGSYN){
    for(int i=0; i<ProblemSpace.synrxns.rxns.size();i++) {
      if(ProblemSpace.synrxns.rxns[i].syn.size() > 1) {
	printf("SYN %d: ", ProblemSpace.synrxns.rxns[i].id);
	for(int j=0;j<ProblemSpace.synrxns.rxns[i].syn.size();j++) {
	  printf("%s(%d) ;  ", 
		 ProblemSpace.fullrxns.rxnFromId(ProblemSpace.synrxns.rxns[i].syn[j]).name, 
		 ProblemSpace.synrxns.rxns[i].syn[j]);
	}
	printf("\n");
      }
    }
  }

  if(_db.DEBUGSYN) {
    printf("METRXNRELATIONS:\n");
    for(int i=0;i<ProblemSpace.metabolites.mets.size();i++) {
      printf("%s:  ", ProblemSpace.metabolites.mets[i].name);
      printIntVector(ProblemSpace.metabolites.mets[i].rxnsInvolved_nosec);
    }
  }

  /* Load synrxnsR (syn reactions in opposite direction) */
  GoodReversible(ProblemSpace);

  for(int i=0;i<ProblemSpace.metabolites.mets.size();i++){
    if(ProblemSpace.metabolites.mets[i].noncentral==1){
      ProblemSpace.cofactors.addMetabolite(ProblemSpace.metabolites.mets[i]);
    }
  }

  /* Load cofactor list */
  for(int i=0; i<ProblemSpace.metabolites.mets.size(); i++) {
    if(ProblemSpace.metabolites.mets[i].secondary_lone == 1) {
      ProblemSpace.secondaryLones.addMetabolite(ProblemSpace.metabolites.mets[i]);
    }
  }

  /* Needed for ETC chain - find out what reactions each noncentral cofactor (ProblemSpace.cofactors) comes from */
  calcMetRxnRelations(ProblemSpace.fullrxns,ProblemSpace.cofactors);
  calcMetRxnRelations(ProblemSpace.fullrxns,ProblemSpace.metabolites);

  /* Add the full biomass equation */
  REACTION fullBiomass;
  /* FIXME: Only looks at 0'th growth for now... */
  fullBiomass = MakeBiomassFromGrowth(ProblemSpace.growth[0]);
  ProblemSpace.fullrxns.addReaction(fullBiomass);


  /* Ensure that GAM and NGAM can be properly tuned (applies to fullrxns only) */
  setUpMaintenanceReactions(ProblemSpace);

  /* Adjust likelihoods into Dijkstras-compatible distances (including accounting for special flags for spontaneous,
     etc. */
  if(_db.PARSIMONY) {
    for(int i=0; i<ProblemSpace.synrxns.rxns.size(); i++) { ProblemSpace.synrxns.rxns[i].current_likelihood = 1;}
    for(int i=0; i<ProblemSpace.synrxnsR.rxns.size(); i++) { ProblemSpace.synrxnsR.rxns[i].current_likelihood = 1; }
    for(int i=0; i<ProblemSpace.fullrxns.rxns.size(); i++) { ProblemSpace.fullrxns.rxns[i].current_likelihood = 1; }
  } else {
    /* void adjustLikelihoods(vector<REACTION> &rxnList, double spont_likely, double black_magic_likely,
       double hard_include_likely, double no_likely, bool adjustNonspecial) */
    adjustLikelihoods(ProblemSpace.synrxns.rxns, 1.0f, -3.0f, 1.1f, -10.0f, true);
    adjustLikelihoods(ProblemSpace.synrxnsR.rxns,  1.0f, -3.0f, 1.1f, -10.0f, true);
    adjustLikelihoods(ProblemSpace.fullrxns.rxns, 1.0f, -3.0f, 1.1f, -10.0f, true);
  }

  return;
}

/* Function for finding shortest paths to an output */
void Run_K(PROBLEM &ProblemSpace, vector<MEDIA> &media, RXNSPACE &rxnspace, int outputId, int K, vector<PATHSUMMARY> &result, int direction, int growthIdx){

  vector<int> inputIds;
  vector<PATH> kpaths;
  int Kq(K);

  for(int i=0;i<media.size();i++){
    inputIds.push_back(media[i].id);
  }

  if(_db.DEBUGRUNK) {
    printf("Number of reactions passed to Run_K: %d\n", (int)rxnspace.rxns.size());
    printf("Output ID: %d\n", outputId);
    printf("Input IDs:"); printIntVector(inputIds);
  }

  METSPACE inputs(ProblemSpace.metabolites, inputIds);  
  kShortest(kpaths,rxnspace,ProblemSpace.metabolites, inputs,
	    ProblemSpace.metabolites.metFromId(outputId), Kq);


  /* Optional Intermediate Print Statement */
  if(_db.DEBUGPATHS){
    printPathResults(kpaths,ProblemSpace,rxnspace);
  }

  if(_db.VISUALIZEPATHS) {
    #pragma omp parallel for
    for(int j=0;j<kpaths.size();j++){
      char s[128];
      sprintf(s, "./outputs/path_%s_k%d.dot", ProblemSpace.metabolites.metFromId(outputId).name, j);
      VisualizePath2File(s, s, kpaths[j], ProblemSpace, 1);
    }
  }

  for(int i=0;i<kpaths.size();i++){
    PATHSUMMARY temp;
    if(direction==1){
      temp = kpaths[i];
      temp.growthIdx.push_back(growthIdx);
      temp.k_number = i;
      /* outputId can be from 0 to 10000 or so, growthIdx from 0 to 100 or so and K from 0 to 10 or so. So this should avoid any possible conflicts */
      temp.id = 10000000*temp.k_number + 100000*growthIdx + temp.outputId;
      result.push_back(temp);
    }
    if(direction==-1){
      temp / kpaths[i];
      temp.growthIdx.push_back(growthIdx);
      temp.k_number = i;
      temp.id = 10000000*temp.k_number + 100000*growthIdx + temp.outputId;
      result.push_back(temp);
    }
  }
  return;
}

/* Function for finding shortest paths to an output using kShortest2
   Again, bad practice to copy a function over like this, but it gives
   me a completely separate space to work in. */
void Run_K2(PROBLEM &ProblemSpace, vector<MEDIA> &media, int outputId, int K, int startingK,
	    vector<PATHSUMMARY> &result, int direction, int growthIdx){
  //printf("enter\n");fflush(stdout);
  vector<int> inputIds;
  vector<PATH> kpaths;
  int Kq(K);
  RXNSPACE &rxnspace = ProblemSpace.synrxnsR;

  for(int i=0;i<media.size();i++){
    inputIds.push_back(media[i].id);
  }

  if(_db.DEBUGRUNK) {
    printf("Number of reactions passed to Run_K2: %d\n", (int)rxnspace.rxns.size());
    printf("Output ID: %d\n", outputId);
    printf("Input IDs:"); printIntVector(inputIds);
  }

  METSPACE inputs(ProblemSpace.metabolites, inputIds);  
  kShortest2(kpaths,rxnspace,ProblemSpace.metabolites, inputs,
	     ProblemSpace.metabolites.metFromId(outputId), Kq,
	     ProblemSpace.synrxnsR);

  /* Optional Intermediate Print Statement */
  if(_db.DEBUGPATHS){  printPathResults(kpaths,ProblemSpace,rxnspace);  }

  for(int i=0;i<kpaths.size();i++){
    PATHSUMMARY temp;
    if(direction==1){
      temp = kpaths[i];
      temp.growthIdx.push_back(growthIdx);
      temp.k_number = i + startingK;
      /* outputId can be from 0 to 10000 or so, growthIdx from 0 to 100 or so and K from 0 to 10 or so. So this should avoid any possible conflicts */
      temp.id = 1000000*temp.outputId + 1000*temp.growthIdx[0] + temp.k_number;
      result.push_back(temp);
    }
    if(direction==-1){
      temp / kpaths[i];
      temp.growthIdx.push_back(growthIdx);
      temp.k_number = i + startingK;
      temp.id = 1000000*temp.outputId + 1000*temp.growthIdx[0] + temp.k_number;
      result.push_back(temp);
    }
  }
  return;
}

/* Run K-shortest in teh forward direction */
void FirstKPass(PROBLEM &ProblemSpace, int K, vector<vector<vector<PATHSUMMARY> > > &psum){
   /* K-Shortest ROUND 1*/
  for(int i=0;i<ProblemSpace.growth.size();i++){
    vector<int> outputIds = Load_Outputs_From_Growth(ProblemSpace, i);
    vector<vector<PATHSUMMARY> > tempP2;
    for(int j=0;j<outputIds.size();j++){
      if(_db.DEBUGPATHS) {
	printf("FirstKPass: growth %d of %d   output %d(%s) of %d\n",i+1,(int)ProblemSpace.growth.size(),j+1,
	       ProblemSpace.metabolites.metFromId(outputIds[j]).name,
	       (int)outputIds.size());
      }
      vector<PATHSUMMARY> tempP1;
      Run_K(ProblemSpace,ProblemSpace.growth[i].media,ProblemSpace.synrxns,outputIds[j],K,tempP1,1, i);
      if(_db.DEBUGPATHS) { printf("FirstKPass: %d paths found\n",(int)tempP1.size()); }
      tempP2.push_back(tempP1);
    }
    psum.push_back(tempP2);
  }
  return;
}

/* Run K-shortest in the reverse direction (everything except transporters is allowed to go the opposite way, but with a penalty)  */
void SecondKPass(PROBLEM &ProblemSpace, int K, vector<vector<vector<PATHSUMMARY> > > &psum){

  /* Output List */
  vector<vector<int> > outputIds;
  vector<GROWTH> &growth         = ProblemSpace.growth;

  /* Note - the likelihoods have been modified here to penalize those reactions with modified reversibilities.
     The likelihoods are changed to BLACK MAGIC */
  RXNSPACE &synrxns             = ProblemSpace.synrxns;
  RXNSPACE &synrxnsR            = ProblemSpace.synrxnsR;
  RXNSPACE &fullRxns            = ProblemSpace.fullrxns;

  vector<int> dum;
  
  for(int i=0;i<growth.size();i++){
    vector<int> outputId = Load_Outputs_From_Growth(ProblemSpace, i);
    outputIds.push_back(outputId);
  }

  vector<int> idsToMakeReversible;

  /* Rerun KShortest on the outputs that could not be found */
  for(int i=0;i<psum.size();i++){
    for(int j=0;j<psum[i].size();j++){
      if(_db.DEBUGPATHS) {
	printf("SecondKPass: growth %d of %d   output %s (%d of %d)\n",i+1,(int)psum.size(),
	       ProblemSpace.metabolites.metFromId(outputIds[i][j]).name, 
	       j+1,(int)psum[i].size());
      }

      /* Less paths were originally found than requested. */
      if(psum[i][j].size() < K){
	vector<PATHSUMMARY> temp_psum;
	Run_K2(ProblemSpace,growth[i].media,outputIds[i][j],K-psum[i][j].size(), psum[i][j].size(),
	       temp_psum,1,i);

	for(int l=0;l<temp_psum.size();l++){  psum[i][j].push_back(temp_psum[l]);  }
      } /* if psum[i][k].size() < K */
    } 
  } /* for i=1:psum.size() */

  return;
}

REACTION MagicExit(const vector<REACTION> &reaction, int met_id, const char* name){
  REACTION newRxn;
  /* Irreversible by default */
  newRxn = MagicExit(reaction,met_id,name,1);
  return newRxn;
}

/* Could probably flush this out a bit more. R: Reversibility 
Changed to return a REACTION to make it more portable / save on memory in various places */
REACTION MagicExit(const vector<REACTION> &reaction, int met_id, const char* name, int R){
  REACTION rxn_add;
  char *temp = (char *) malloc(sizeof(char) * 64);
  STOICH stoich_add;

  rxn_add.init_likelihood = -3;
  rxn_add.isExchange = 1;
  rxn_add.net_reversible = R;
  rxn_add.id = _db.BLACKMAGICFACTOR + met_id;

  sprintf(temp,"%s%s","MagicExit_",name);
  strcpy(rxn_add.name,temp);


  stoich_add.met_id = met_id;
  stoich_add.rxn_coeff = -1;
  sprintf(stoich_add.met_name,"%s",name);
  rxn_add.stoich.push_back(stoich_add);
  rxn_add.stoich_part.push_back(stoich_add);

  return rxn_add;
}

/* Could probably flush this out a bit more - should never use this except to add exchanges
for secondary metabolites. The difference from MagicExit is that this one includes a flux
bound AND adds a transporter. MagicExit just uses the default. Also the init_likelihood is different.  */
REACTION GrowthExit(const vector<REACTION> &reaction, int met_id, int reversible, double fluxBound, const char* name){
  REACTION rxn_add;
  STOICH stoich_add;
  char* temp = (char*) malloc( sizeof(char)*64 ); 

  rxn_add.init_likelihood = -2;
  rxn_add.isExchange = 1;

  /* could be an issue in case of long met names */
  sprintf(temp, "%s%s", "GrowthExit_", name); 
  strcpy(rxn_add.name,temp);
 //printf("Made GrowthExit: %s\n",rxn_add.name);
  rxn_add.init_reversible = reversible;
  rxn_add.net_reversible = reversible;
  rxn_add.id = _db.BLACKMAGICFACTOR + met_id;
  stoich_add.met_id = met_id;
  stoich_add.rxn_coeff = -1;
  /* No need to cut off here, name same from a size-64 string and being printed into the same */
  sprintf(stoich_add.met_name,"%s",name);
  rxn_add.stoich.push_back(stoich_add);
  rxn_add.stoich_part.push_back(stoich_add);
  if(reversible > 0) { rxn_add.lb = 0.0f; } else { rxn_add.lb = -fluxBound; }
  if(reversible < 0) { rxn_add.ub = 0.0f; } else { rxn_add.ub = fluxBound;  }

  free(temp);
  return rxn_add;
}

/* More general function than AddBiomass - add a arbitrary combination of metabolites into a "biomass equation" with coefficients and IDs
as specified in the vector<stoich> */
REACTION MakeObjRxn(const vector<STOICH> &stoichVec) {
  REACTION obj;
  obj.id = _db.BIOMASS;
  strcpy(obj.name, "BIOMASS_CUST");
  obj.stoich = stoichVec;
  obj.stoich_part = stoichVec;
  obj.net_reversible = 1;
  obj.init_reversible = 1;
  /* Set these to make it consistent with the reversibility */
  obj.lb = 0.0f;  obj.ub = 1000.0f;
  obj.init_likelihood = -2;
  obj.current_likelihood = 1.0f;

  return obj;
}

REACTION MakeBiomassFromGrowth(const GROWTH &growth){
  return MakeObjRxn(growth.biomass);
}

/* Given an metabolite with id met_id, finds the corresponding 
   metabolite in the opposite compartment. Returns the matching met ID or -1 if it fails to find a match */
int inOutPair(int met_id, const METSPACE &metspace){ 
  unsigned int i;
  int flag;
  bool isExternal = isExternalMet(metspace.metFromId(met_id).name, _db.E_tag);
  char tempS[64] = {0};
  if(isExternal){
    strncpy(tempS,metspace.metFromId(met_id).name,(int)strlen(metspace.metFromId(met_id).name)-3);
  } else {
    strcat(tempS,metspace.metFromId(met_id).name);
    strcat(tempS,_db.E_tag);
  }

  for(i=0;i<metspace.mets.size();i++){
    if(strcmp(tempS,metspace.mets[i].name)==0){
      return metspace.mets[i].id;
    } 
  }

  return -1;
}

vector<int> Load_Outputs_From_Growth(const PROBLEM &parentSpace, int growthIndex) {
  /* Load outputs from a GROWTH
     START = -1 means pull out everything 
     Put this into a function and call it to avoid possible consistency problems between Dijkstra runs */
  const GROWTH &growth = parentSpace.growth[growthIndex];
  vector<int> outputIds;
  for(int j=0;j<growth.biomass.size();j++) {
    outputIds.push_back(growth.biomass[j].met_id);}
  for(int j=0;j<growth.byproduct.size();j++) {
    outputIds.push_back(growth.byproduct[j].id);}
  for(int j=0; j<parentSpace.secondaryLones.mets.size(); j++) {
    outputIds.push_back(parentSpace.secondaryLones.mets[j].id);
  }
  custom_unique(outputIds);
  return outputIds;
}

vector<int> Load_Inputs_From_Growth(const GROWTH &growth) {
  vector<int> inputIds;
  for(int j=0; j<growth.media.size(); j++) {
    inputIds.push_back(growth.media[j].id);
  }
  return inputIds;
}

/* Add opposite-direction version of every reactino to synrxnsR 
 Give it a BLACK MAGIC flag to make it less likely to occur than the "bad direction" version */
void GoodReversible(PROBLEM &ProblemSpace){
  ProblemSpace.synrxnsR.rxns.reserve(ProblemSpace.synrxns.rxns.size()*2);
  ProblemSpace.synrxnsR.addReactionVector(ProblemSpace.synrxns.rxns);
  for(int i=0;i<ProblemSpace.synrxns.rxns.size();i++){
    if(ProblemSpace.synrxnsR.rxns[i].isExchange==0){
      if(ProblemSpace.synrxnsR.rxns[i].init_reversible != 0) {
	/* BLACK MAGIC [this causes secondKpass to 
	 give a significant likelihood penalty to reactions that have
	changed reversibility - hopefully this prevents too many of those
	from appearing!] */
	REACTION temp;
	temp = ProblemSpace.synrxnsR.rxns[i];
	temp.id += _db.REVFACTOR;
	temp.init_likelihood = -3; /* BLACK MAGIC */
	strcat(temp.name,"_REV");
	ProblemSpace.synrxnsR.addReaction(temp);
	ProblemSpace.synrxnsR.changeReversibility(temp.id, temp.init_reversible * -1);
      }
    }
  }
  return;
}

/* Reverse reversibility of all reactions (for gapfilling) */
void ReverseReversible(vector<REACTION> &reaction){
  for(int i=0;i<reaction.size();i++){
    reaction[i].net_reversible *= -1;
  }
  return;
}

int in(int id, const vector<int> &list){
  for(int i=0;i<list.size();i++){
    if(list[i]==id){ return 1;}
    if(list[i]==-id){ return -1;}
  }
  return 0;
}

/* We need to convert the switches for special reactions into likelihoods for use with Dijkstras. Dijkstras DOES NOT WORK with negative edge (except things I hard-coded excpetions for)
Modified 11-27-10: Input now already is in terms of probabilities.
HARD NO INCLUDE (-1):  I keep at -1 and deal with it directly in the shortest path algorithm...IMPORTANT for K-shortest paths.
HARD INCLUDE (-2): 1 or higher (maybe 2?)
BLACK MAGIC (-3): black_magic_likely [suggested 0 or somethign negative (-1 or -2)]
SPONTANEOUS (-4): spont_likely [suggested: close to 1]
NO LIKELIHOODS (-5): no_likely [suggested: some negative number, probably -1 or -2]

Everything else: just leave it the way it is

Suggested sample: adjustLikelihoods(rxnList, 1.0f, -2.0f, 1.1f, -2.0f, true) */

void adjustLikelihoods(vector<REACTION> &rxnList, double spont_likely, double black_magic_likely,
                       double hard_include_likely, double no_likely, bool adjustNonspecial) {
  unsigned int i;
  for(i=0;i<rxnList.size();i++) {
    if(rxnList[i].init_likelihood < -4.1) {
      rxnList[i].current_likelihood = hard_include_likely - no_likely;
    }
    else if(rxnList[i].init_likelihood < -3.1) {
      rxnList[i].current_likelihood = hard_include_likely - spont_likely;
    }
    else if(rxnList[i].init_likelihood < -2.1) {
      rxnList[i].current_likelihood = hard_include_likely - black_magic_likely;
    }
    else if(rxnList[i].init_likelihood < -1.1) {
      rxnList[i].current_likelihood = hard_include_likely - hard_include_likely;
    }
    else if (rxnList[i].init_likelihood < -0.1) {
      /* hard DO NOT INCLUDE - DO NOT DO ANYTHING (dealt with directly in Dijkstras algorithm) */
    }
    else {
      // FIX 2-21-11: Subtract init_likelihood instead of current_likelihood?
      if(adjustNonspecial) {
        rxnList[i].current_likelihood = hard_include_likely - rxnList[i].init_likelihood;
      }
    }
  }
}

/* Same as below, but with stoich instead of stoich_part */
void calcMetRxnRelations(const RXNSPACE &rxnspace, METSPACE &metspace) {

  unsigned int rxnSize = rxnspace.rxns.size();
  unsigned int metSize = metspace.mets.size();
  DEBUGFLAGS db;

  /* Clear out any existing rxnsInvolved */
  for(int i=0;i<metSize;i++) {
    metspace.mets[i].rxnsInvolved_nosec.clear();
  }

  for(int i=0;i<rxnSize;i++) {
    const REACTION &currentRxn = rxnspace.rxns[i];
    for(int j=0;j<currentRxn.stoich.size();j++) {
      //      if(currentRxn.id < _db.MAGICBRIDGEFACTOR || currentRxn.id > _db.MAGICBRIDGEFACTOR + _db.MINFACTORSPACING){
      if(metspace.idIn(currentRxn.stoich[j].met_id)) {
	metspace.mets[metspace.idxFromId(currentRxn.stoich[j].met_id)].rxnsInvolved_nosec.push_back(currentRxn.id);
      }
      // }
    }
  }

}

/* Same as above, but using stoich_part instead of stoich [so that reactions in which a metabolite is
   a secondary_lone or secondary_pair are excluded] */
void calcMetRxnRelations_nosec(const RXNSPACE &rxnspace, METSPACE &metspace) {

  unsigned int rxnSize = rxnspace.rxns.size();
  unsigned int metSize = metspace.mets.size();

  /* Clear out any existing rxnsInvolved */
  //printf("Note - calcMetRxnRelations_nosec clears out any existing relationships!\n");
  for(int i=0;i<metSize;i++) {
    metspace.mets[i].rxnsInvolved_nosec.clear();
  }

  for(int i=0;i<rxnSize;i++) {
    const REACTION& currentRxn = rxnspace.rxns[i];
    for(int j=0;j<currentRxn.stoich_part.size();j++) {
      int metIdx = metspace.idxFromId(currentRxn.stoich_part[j].met_id);
      metspace.mets[metIdx].rxnsInvolved_nosec.push_back(currentRxn.id);
    }
  }
}
