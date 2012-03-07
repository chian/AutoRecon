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

  /* Generate a list of fullrxn and synrxn information to curate */
  FILE* fid = fopen("FULL_SYNRXNS", "w");
  FILE* fid2 = fopen("RXN_LIKELIHOODS", "w");
  RXNSPACE &tmpfull = ProblemSpace.fullrxns;
  /*  char fullString[4086];
  char partString[4086];
  for(int i=0; i<tmpfull.rxns.size(); i++) {
    fullString[0] = '\0';
    partString[0] = '\0';
    printRxnFormula(tmpfull.rxns[i], fullString, true);
    printRxnFormula(tmpfull.rxns[i], partString, false);
    fprintf(fid, "%s\t%d\t%1.3f\t%s\t%s\n", 
	    tmpfull.rxns[i].name, tmpfull.rxns[i].net_reversible, tmpfull.rxns[i].current_likelihood,
	    fullString, partString);
	    } */

  /* Make a computer-friendly format so I can easily curate, add/remove things */
  for(int i=0; i<tmpfull.rxns.size(); i++) {
    for(int j=0; j<tmpfull.rxns[i].stoich.size(); j++) {
      bool keep = false;
      for(int k=0; k<tmpfull.rxns[i].stoich_part.size(); k++) {
	if(tmpfull.rxns[i].stoich[j].met_id == tmpfull.rxns[i].stoich_part[k].met_id) { keep = true; break; }
      }
      /* Reaction name - Reaction ID - Reversibility - Metabolite name - Metabolite ID - reaction coefficient - Secondary or not [1 = secondary] */
      fprintf(fid, "%s\t%d\t%d\t%s\t%d\t%1.8f\t%d\n",
	      tmpfull.rxns[i].name, tmpfull.rxns[i].id, tmpfull.rxns[i].net_reversible, 
	      tmpfull.rxns[i].stoich[j].met_name, tmpfull.rxns[i].stoich[j].met_id, tmpfull.rxns[i].stoich[j].rxn_coeff, keep?0:1);
      fprintf(fid2, "%s\t%1.4f\n", tmpfull.rxns[i].name, tmpfull.rxns[i].init_likelihood);
    }
  }

  fclose(fid); fclose(fid2);

  printSynRxns(ProblemSpace.synrxns, ProblemSpace.fullrxns);

  /* I tried to move this to inputSetup as well but it gave me a compile error
     so here it is*/
  printf("Annotating genes...\n");
  map<string, vector<VALUESTORE> > annoteToRxns = GeneAnnotations(ProblemSpace.fullrxns, _db.ANNOTE_CUTOFF_1, _db.ANNOTE_CUTOFF_2);
  printf("...done\n");

  //  printREACTIONvector(ProblemSpace.synrxns.rxns, 1);

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

  PrintPathSummary2(unSynPsum, ProblemSpace.fullrxns);

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
  ANSWER ans = gapfillWrapper(ProblemSpace, flat_unsyn, ProblemSpace.growth[0]);

    //  gapFindGapFill(allMagic, ProblemSpace);
  printf("...done.\n");

  //vector<double> fbaResult = FBA_SOLVE(ans,flat,ProblemSpace); //make it take ANSWER and go from there...

  printf("No seg faults!\n");
  return 0;

}

