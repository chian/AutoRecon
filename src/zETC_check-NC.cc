#include "DataStructures.h"
#include "RunK.h"
#include "shortestPath.h"
#include "genericLinprog.h"
#include "pathUtils.h"
#include "visual01.h"
#include "Grow.h"
#include "Paths2Model.h"
#include "XML_loader.h"
#include "ETC.h"
#include "Exchanges.h"
#include "RandomMagic.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <glpk.h>
#include <map>
#include <vector>

using std::vector;
using std::map;
using std::string;

int main(int argc, char *argv[]) {

  /* Debug Statements */
  int DEBUG = 0; // For error checking, -1 mean none

  /* XML Load Variables */
  PROBLEM ProblemSpace;

  /* Will be TRUE if the ETC was found to be connected to the rest of the network and FALSE otherwise... */
  vector<bool> isConnected;

  /* K shortest Variables */
  int K = 4;
  int Kd = 1;
  vector<vector<vector<PATHSUMMARY> > > psum;

  /* Turn off GLPK printing */
  int (*func)(void*, const char *) = &suppressGLPKOutput;
  glp_term_hook(func, NULL);

  /* Setup done... */
  printf("Entering InputSetup...\n");fflush(stdout);
  InputSetup(argc,argv,ProblemSpace);
  printf("...exiting InputSetup\nEntering FirstKPass...\n");fflush(stdout);

  vector<GROWTH> &growth = ProblemSpace.growth;
  vector<int> dirs;
  vector<int> metsToCheck;
  for(int i=0;i<growth.size();i++) {
    for(int j=0;j<growth[i].media.size();j++) {
      dirs.push_back(0); /* Note we want to allow media to be imported OR exported... (especially important for H+ because transporters depend on it) */
      metsToCheck.push_back(growth[i].media[j].id);
    }
    for(int j=0;j<growth[i].byproduct.size();j++) {
      dirs.push_back(1); /* ... but anythign NOT in the media must be only exported */
      metsToCheck.push_back(growth[i].byproduct[j].id);
    }
  }

  PROBLEM TMPproblem(ProblemSpace.fullrxns, ProblemSpace);
  checkExchangesAndTransporters(ProblemSpace, TMPproblem, metsToCheck, dirs);
  TMPproblem.clear();

  /* Add the full biomass equation (note - this is needed to test if something is produced or destroyed in the
     biomass equation, but be careful of ID conflicts with the "mini" BIOMASS reactions in individual PLIST's! */
  REACTION fullBiomass;
  /* FIXME: Only looks at 0'th growth for now... */
  fullBiomass = MakeBiomassFromGrowth(ProblemSpace.growth[0]);
  ProblemSpace.fullrxns.addReaction(fullBiomass);

  /* Re-adjust likelihoods in case we added something that is tagged as BLACK MAGIC */
  adjustLikelihoods(ProblemSpace.fullrxns.rxns, 1.0f, -2.0f, 1.1f, -2.0f, true);
  adjustLikelihoods(ProblemSpace.synrxns.rxns, 1.0f, -2.0f, 1.1f, -2.0f, true);
  adjustLikelihoods(ProblemSpace.synrxnsR.rxns,  1.0f, -2.0f, 1.1f, -2.0f, true);

  FirstKPass(ProblemSpace,K,psum);

  printf("...exiting FirstKPass\nEntering SecondKPass...\n");fflush(stdout);
  SecondKPass(ProblemSpace,K,psum);
  if(DEBUG == 0) {
    PrintPathSummary2(psum,ProblemSpace); }

  /* ETC */
  vector<NETREACTION> ETCout;
  double cutoff = 0.0;

  printf("Entering ETC...");fflush(stdout);
  ETC(ProblemSpace,ETCout,cutoff,psum);
  printf("...exiting ETC\n");

  // BEGIN - Print ETC results 
  RXNSPACE &rxnspace = ProblemSpace.fullrxns;

  for(int i=0;i<ETCout.size();i++){
    printf("chain %d of %d: ",i,(int)ETCout.size());
    for(int j=0;j<ETCout[i].rxnDirIds.size();j++){
      printf("%s(%d)  ",rxnspace.rxnFromId(abs(ETCout[i].rxnDirIds[j])).name,ETCout[i].rxnDirIds[j]);
    }
    printf("\n");
  }

  vector<NETREACTION> ETCout2 = ETC_dir_check(ProblemSpace.fullrxns,ETCout);

  printf("\nFilter ETCs for Direction:\n");
  for(int i=0;i<ETCout2.size();i++){
    printf("chain %d of %d (%d): ",i,(int)ETCout2.size(),ETCout2[i].rxn.net_reversible);
    for(int j=0;j<ETCout2[i].rxnDirIds.size();j++){
      printf("%s(%1.3f)  ",rxnspace.rxnFromId(abs(ETCout2[i].rxnDirIds[j])).name,
	     rxnspace.rxnFromId(abs(ETCout2[i].rxnDirIds[j])).init_likelihood);
    }
    printf("\n");fflush(stdout);
  }
  // END - Print ETC results 
  printf("Seeing which ETCs connect...\n");fflush(stdout);

  /* From now on we're dealing only with the fullrxns, not the synrxns. Therefore, we should
     adjust the rxnsInvolved to use fullrxns, not the synrxns */
  calcMetRxnRelations(ProblemSpace.fullrxns, ProblemSpace.metabolites);

  /* Convert to the "Flat" convention. Flat pathsummary contains pointers to the original PATHSUMMARY in psum */
  vector<PATHSUMMARY> flat;
  flat = flattenPsum(psum);

  vector<PATHSUMMARY> unSynPsum;
  for(int i=0;i<flat.size();i++) {
    unSynPsum.push_back(replaceWithRealRxnIds(flat[i], ProblemSpace));
  }

  ETC_connect(ProblemSpace, ETCout2, unSynPsum);

  printf("Seeing which ETCs connect...done\n");fflush(stdout);
  

  printf("No seg faults!\n");
  return 0;
}

