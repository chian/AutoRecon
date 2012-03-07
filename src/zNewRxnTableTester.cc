#include "Annotations.h"
#include "DataStructures.h"
#include "Grow.h"
#include "MyConstants.h"
#include "Paths2Model.h"
#include "Printers.h"
#include "RunK.h"
#include "TableLoader.h"

#include <cstdio>
#include <vector>

using std::vector;

int main(int argc, char *argv[]) {
  printf("\n\n");

  /* Until we can get growth stuff sorted out I still read that from the XML file */
  PROBLEM tmp;
  InputSetup(argc,argv,tmp);

  PROBLEM ProblemSpace = readReactionTable("FULL_SYNRXNS_MOD");
  readLikelihoodTable(ProblemSpace, "RXN_LIKELIHOODS");
  /* FIXME - will read the list of ETC cofactor pairs from a file... */

  MakeSynList(ProblemSpace);
  ProblemSpace.growth = tmp.growth;
  ProblemSpace.cofactors = tmp.cofactors;

  /* Needed for ETC chain - find out what reactions each noncentral cofactor (ProblemSpace.cofactors) comes from */
  calcMetRxnRelations(ProblemSpace.fullrxns,ProblemSpace.cofactors);
  calcMetRxnRelations(ProblemSpace.fullrxns,ProblemSpace.metabolites);

  tmp.clear();

  adjustLikelihoods(ProblemSpace.synrxns.rxns, 1.0f, -3.0f, 1.1f, -10.0f, true);
  //  adjustLikelihoods(ProblemSpace.synrxnsR.rxns,  1.0f, -3.0f, 1.1f, -10.0f, true);
  adjustLikelihoods(ProblemSpace.fullrxns.rxns, 1.0f, -3.0f, 1.1f, -10.0f, true);  

  //  printREACTIONvector(ProblemSpace.synrxns.rxns, 1);
  printSynRxns(ProblemSpace.synrxns, ProblemSpace.fullrxns);

  map<string, vector<VALUESTORE> > annoteToRxns = GeneAnnotations(ProblemSpace.fullrxns, _db.ANNOTE_CUTOFF_1, _db.ANNOTE_CUTOFF_2);

  int K=1;
  vector<vector<vector<PATHSUMMARY> > > psum;
  FirstKPass(ProblemSpace,K,psum);  
  PrintPathSummary(psum);

  /* From now on we're dealing only with the fullrxns, not the synrxns. Therefore, we should
     adjust the rxnsInvolved to use fullrxns, not the synrxns */
  calcMetRxnRelations(ProblemSpace.fullrxns, ProblemSpace.metabolites);

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

  printf("Performing gapfilling...\n");
  vector<NETREACTION> bigout;
  initializeAnswer(bigout, unSynPsum, annoteToRxns);
  ANSWER ans = gapfillWrapper(ProblemSpace, flat_unsyn, ProblemSpace.growth[0]);

  return 0;
}
