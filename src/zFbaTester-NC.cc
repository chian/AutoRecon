#include "DataStructures.h"
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
#include "score.h"
#include "MersenneTwister.h"

#include <map>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using std::vector;
using std::map;
using std::string;

int main(int argc, char *argv[]) {
  printf("\n\n");

  /* Main problem structure containing all the metabolites and reactions
     (in fullrxns), their synonym versions (synrxns), and a reversed synonym version
     useful for finding paths despite reversibility bottlenecks (synrxnsR) 
  
      All three sets of reactions share the same metabolites (metabolites) */
  PROBLEM ProblemSpace;

  /* Number of K-shortest paths (forward if possible, backward otherwise) */
  int K = _db.INITIAL_K;

  /* psum[i][j][k]: Growth condition i, output j, k'th shortest path */
  vector<vector<vector<PATHSUMMARY> > > psum;

  printf("Setting up problem inputs...\n");fflush(stdout);
  InputSetup(argc,argv,ProblemSpace);
  printf("...done \n");fflush(stdout);
  
  /* I tried to move this to inputSetup as well but it gave me a compile error
     so here it is*/
  printf("Annotating genes...\n");
  map<string, vector<VALUESTORE> > annoteToRxns = GeneAnnotations(ProblemSpace.fullrxns, _db.ANNOTE_CUTOFF_1, _db.ANNOTE_CUTOFF_2);
  printf("...done\n");

  printf("Finding paths in forward direction...\n");
  FirstKPass(ProblemSpace,K,psum);
  if(_db.DEBUGPATHS) { PrintPathSummary(psum); }
  printf("...done\n");

  printf("Finding paths in the reverse direction (if necessary)...\n");fflush(stdout);
  SecondKPass(ProblemSpace,K,psum);
  if(_db.DEBUGPATHS) { PrintPathSummary2(psum, ProblemSpace.synrxnsR); }
  printf("...done\n");

  /* Treat any metabolites that have no paths found as secondary_lones
     This could make magic entrances under all conditions but we should be able to prune out if any are only
     needed conditionally... */
  for(int i=0; i<psum.size(); i++) {
    vector<int> outputIds = Load_Outputs_From_Growth(ProblemSpace, i);
    for(int j=0; j<psum[i].size(); j++) { 
      if(psum[i][j].size() == 0) {
	printf("WARNING: No paths found to output %s under growth condition %d so making it a secondary_lone \n", ProblemSpace.metabolites.metFromId(outputIds[j]).name, i);
	ProblemSpace.metabolites.metPtrFromId(outputIds[j])->secondary_lone = 1;
      }
    }
  }

  /* Find ETC chains here but connect them within the gapfillWrapper */
  vector<NETREACTION> bigout;
  ETC(ProblemSpace, bigout, -1000.0f);

  printf("Unsynonymizing path summaries...\n");
  /* Add reversible versions of the reactions that had to be reversed */
  addR(psum,ProblemSpace);

  /* From now on we're dealing only with the fullrxns, not the synrxns. Therefore, we should
     adjust the rxnsInvolved to use fullrxns, not the synrxns */
  calcMetRxnRelations(ProblemSpace.fullrxns, ProblemSpace.metabolites);

  /* Convert to the "Flat" convention. Flat pathsummary contains pointers to the original PATHSUMMARY in psum */
  //  vector<PATHSUMMARY> flat, flat_temp;
  //  flat_temp = flattenPsum(psum);
  //  flat = uniquePsum(flat_temp);

  /* This is just temporary...I use it to help with filling magic bridges but beyond that it isn't really
     needed. Elsewhere, we should just use the unsynPsum below */
  vector<PATHSUMMARY> flat = flattenPsum(psum);

  vector<vector<vector<PATHSUMMARY> > > unSynPsum = psum;
  //  for(int i=0;i<flat.size();i++) {  unSynPsum.push_back(replaceWithRealRxnIds(flat[i], flat, ProblemSpace));  }
  for(int i=0; i<psum.size(); i++) {
    for(int j=0; j<psum[i].size(); j++) {
      for(int k=0; k<psum[i][j].size(); k++) {
	unSynPsum[i][j][k] = replaceWithRealRxnIds(psum[i][j][k], flat, ProblemSpace);
      }
    }
  }

  vector<PATHSUMMARY> flat_unsyn = flattenPsum(unSynPsum);

  printf("...done\n");

  printf("Making pre-gapfill results files...\n");
  /* Visualize if requested */
  if(_db.VISUALIZEPATHS) {  visualizePathSummary2File("./outputs/UnSynPaths", "UnSynPaths", flat_unsyn, ProblemSpace, false);  }
  if(_db.OUTPUTPATHRESULTS) {
    PATHS_rxns_out("./outputs/PATH_rxns_report", flat_unsyn, ProblemSpace);
    PATHS_mets_out("./outputs/PATH_mets_report", flat_unsyn, ProblemSpace);
  }
  printf("...done\n");

  printf("Performing gapfilling...\n");
  /**********************************
   NICK: This is the function you will call when you are running the genetic algorithm
   I think ETC-connect will actually go in here
   (ETC finding should be done elsewhere though, like in the main function)
  *********************************/

  initializeAnswer(bigout, unSynPsum, annoteToRxns);
  ANSWER ans1 = gapfillWrapper(ProblemSpace, flat_unsyn, ProblemSpace.growth[0]);
  addATPM(ProblemSpace,ans1);

  vector<ANSWER> ans;
  ans.push_back(ans1);
  SCORE1 score1;
  printf("measureScore...\n");
  //double sc = measureScore(score1, ans, ProblemSpace.growth);
  //printf("%f\n", score1.score);

  /* TEST */
  vector<ANSWER> *test = &ans;

  optimizeAMs(ans, ProblemSpace, 50);


  /*
  for(int i=0;i<100;i++){
    for(int j=0;j<100;j++){
      SCORE1 score12;
      score12.ngam = 1 * i;
      score12.gam  = 1 * j;
      //printf("measureScore for NGAM = %1.3f, GAM = %1.3f ...\n",score12.ngam,score12.gam);
      double sc = measureScore(score12, ans, ProblemSpace.growth);
      for(int k=0;k<ProblemSpace.growth.size();k++){
	printf("%f\t%f\t%d\t%f\t%f\t%f\n", score12.ngam, score12.gam, k, ProblemSpace.growth[k].growth_rate, score12.growthRate[k], score12.growthScore[k]);
      }
    }
  }
  */

    //  gapFindGapFill(allMagic, ProblemSpace);
  printf("...done.\n");

  //vector<double> fbaResult = FBA_SOLVE(ans,flat,ProblemSpace); //make it take ANSWER and go from there...

  printf("No seg faults!\n");
  return 0;

}

