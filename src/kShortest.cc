#include "DataStructures.h"
#include "kShortest.h"
#include "shortestPath.h"
#include "MyConstants.h"
#include "pathUtils.h"
#include "Printers.h"
#include "RunK.h"

#include <cstdio>
#include <omp.h>
#include <queue>
#include <map>
#include <set>
#include <vector>

using std::set;
using std::queue;
using std::vector;
using std::map;
using std::priority_queue;

/* This version of kShortest is made for the reversibility check in SecondKPass when we are trying to close
   gaps that cannot be closed without altering a reaction. NOTE: It is almost never a good idea to copy large
   chunks of code like this just to change a few lines, but this helps me avoid breaking anything while I
   figure this out. --NICK */
void kShortest2(vector<PATH> &result, RXNSPACE &rxnspace, METSPACE &metspace, METSPACE &inputs, METABOLITE output, int K,
		RXNSPACE &truedir) {

  priority_queue<GRAPHSTORE> L;
  GRAPHSTORE tmp;
  calcMetRxnRelations_nosec(rxnspace, metspace);

  set<BADIDSTORE> badIds;
  PATH onePath = findShortestPath(rxnspace, metspace, inputs, output, badIds);
  vector<PATH> tmpPath;
  tmpPath.push_back(onePath);

  int currentK = 0;
  PATH badPath;

  /* Since this is the first one we leave the excluded reactions empty */
  tmp.path = onePath;
  L.push(tmp);

  vector<int> previousGraphRxns;

  while(true) {

    if(L.size() == 0) { 
      //      printf("WARNING: The total number of K-shortest solutions available was less than K! Output name: %s \n\n", output.name);
      return;
    }

    GRAPHSTORE currentGraph = L.top(); /* Note - automatically takes out the shortest one */
    L.pop();

    /* Test if the current graph is the same as the previous best. If it is, just skip it because we already
       found all the subgraphs for that, and we don't want repeat solutions...
       We only need to test the previous best because ones before that will be different by induction 

       Note: repeat solutions occur if two rxns are part of the same pathway, we remove either one and it finds the same
       optimal away to get around them in either case. */
    /* TWO NOTES: (by NICK)
       1. This looks to me as if it only prevents tandem repeats, i.e., only if it is the same as the last one, not
          if it's ever been done before (which would be more efficient, I would think)
       2. How do you know, in the next iteration, if it's okay to exclude the reaction from the first case? You could
          also have differences depending on which pair of reactions you exclude even though both exclusions gave
	  you the same result in this run.
       NICK TO FIX - remove repeat stopper and get rid of repeats after all the calculations
     */
    
    bool toSkip;
    if(previousGraphRxns.size() == currentGraph.path.rxnIds.size()) {
      toSkip = true;
      for(int i=0; i<previousGraphRxns.size(); i++) {
	if(previousGraphRxns[i]!=currentGraph.path.rxnIds[i]) {
	  toSkip = false;
	  break;
	}
      }
    } else { 
      toSkip = false; 
    }
    if(toSkip) { goto dontsave; }
    

    /* Now that it passed our sanity check, lets consider it our next optimum and then find all the subgraphs again like before */

    result.push_back(currentGraph.path);
    currentK++;

  dontsave:
    /* Already found K shortest */
    if(currentK == K) { return; }

    /* No more graphs to check (meaning the total number of shortest paths is less than K) */
    if(_db.DEBUGPATHS) { printf("Working on the %dth shortest...\n", currentK + 1); }

    vector<int> &currentRxnList = currentGraph.path.rxnIds;
    vector<int> &excludedRxnIds = currentGraph.excludedRxnIds;

    /* Set all of the specified reactions to not be included */
    for(int i=0; i<excludedRxnIds.size();i++) {
      rxnspace.rxnPtrFromId(excludedRxnIds[i])->old_likelihood = rxnspace.rxnPtrFromId(excludedRxnIds[i])->current_likelihood;
      rxnspace.rxnPtrFromId(excludedRxnIds[i])->current_likelihood = -1;
    }

    /* Compute shortest paths for next iteration */
    tmp = currentGraph;
    /* Not allowed to save a reference as PRIVATE so I have to do this: */
    RXNSPACE tmpRxn = rxnspace;

    vector<GRAPHSTORE> temporaryList(currentRxnList.size(), currentGraph);

#pragma omp parallel for shared(L, currentRxnList, temporaryList) firstprivate(tmp, onePath, tmpRxn)
    for(int i=0; i<currentRxnList.size();i++) {
      set<BADIDSTORE> dum;
      int dir = truedir.rxnFromId(currentRxnList[i]).net_reversible;
      if(dir!=tmp.path.rxnDirection[i] || 1){
	tmp.excludedRxnIds.push_back(currentRxnList[i]);
	tmpRxn.rxnPtrFromId(currentRxnList[i])->old_likelihood = tmpRxn.rxnPtrFromId(currentRxnList[i])->current_likelihood;
	tmpRxn.rxnPtrFromId(currentRxnList[i])->current_likelihood = -1;
	tmp.path = findShortestPath(tmpRxn, metspace, inputs, output, dum);
      }
      else{
	tmp.path = badPath;
      }	
      temporaryList[i] = tmp;
      tmp.excludedRxnIds.pop_back();
      tmpRxn.rxnPtrFromId(currentRxnList[i])->current_likelihood = tmpRxn.rxnPtrFromId(currentRxnList[i])->old_likelihood;
    }

    for(int i=0; i<temporaryList.size(); i++) {
      L.push(temporaryList[i]);
    }

    /* Reset excluded reactions for the next run */
    for(int i=0; i < excludedRxnIds.size();i++) {
      rxnspace.rxnPtrFromId(excludedRxnIds[i])->current_likelihood = rxnspace.rxnPtrFromId(excludedRxnIds[i])->old_likelihood;
    }

    previousGraphRxns = currentGraph.path.rxnIds;
    
  } /* Until the end...... */

}

/* K-shortest on multiple outputs */
void kShortest2(vector<vector<PATH> > &result, RXNSPACE &rxnspace, METSPACE &metspace, METSPACE &inputs, METSPACE &outputs, int K,
	       RXNSPACE &truedir) {
  result.clear();
  for(int i=0; i<outputs.mets.size(); i++) {
    vector<PATH> tmpRes;
    kShortest2(tmpRes, rxnspace, metspace, inputs, outputs.mets[i], K, truedir);
    result.push_back(tmpRes);
  }
}

/* K-shortest on multiple outputs */
void kShortest(vector<vector<PATH> > &result, RXNSPACE &rxnspace, METSPACE &metspace, METSPACE &inputs, METSPACE &outputs, int K) {
  result.clear();
  for(int i=0; i<outputs.mets.size(); i++) {
    vector<PATH> tmpRes;
    kShortest(tmpRes, rxnspace, metspace, inputs, outputs.mets[i], K);
    result.push_back(tmpRes);
  }
}

void kShortest(vector<PATH> &result, RXNSPACE &rxnspace, METSPACE &metspace, METSPACE &inputs, METABOLITE output, int K) {

  set<BADIDSTORE> badIds;
  priority_queue<GRAPHSTORE> L;
  GRAPHSTORE tmp;
  calcMetRxnRelations_nosec(rxnspace, metspace);

  /* The badIds will be the same all the time if they are needed [to troubleshoot no path
     found instances] so we only save the results of the first one (K=1) */
  PATH onePath = findShortestPath(rxnspace, metspace, inputs, output, badIds);
  int currentK = 0;

  /* Since this is the first one we leave the excluded reactions empty */
  tmp.path = onePath;
  L.push(tmp);

  vector<int> previousGraphRxns;

  while(true) {

    if(L.size() == 0) { 
      //      printf("WARNING: The total number of K-shortest solutions available was less than K! Output name: %s \n\n", output.name);
      return;
    }

    GRAPHSTORE currentGraph = L.top(); /* Note - automatically takes out the shortest one */
    L.pop();

    /* Test if the current graph is the same as the previous best. If it is, just skip it because we already
       found all the subgraphs for that, and we don't want repeat solutions...
       We only need to test the previous best because ones before that will be different by induction 

       Note: repeat solutions occur if two rxns are part of the same pathway, we remove either one and it finds the same
       optimal away to get around them in either case. */
    
    bool toSkip;
    if(previousGraphRxns.size() == currentGraph.path.rxnIds.size()) {
      toSkip = true;
      for(int i=0; i<previousGraphRxns.size(); i++) {
	if(previousGraphRxns[i]!=currentGraph.path.rxnIds[i]) {
	  toSkip = false;
	  break;
	}
      }
    } else { 
      toSkip = false; 
    }
    if(toSkip) { goto dontsave; }
    
    /* Now that it passed our sanity check, lets consider it our next optimum and then find all the subgraphs again like before */
    result.push_back(currentGraph.path);
    currentK++;

  dontsave:
    /* Already found K shortest */
    if(currentK == K) { return; }

    /* No more graphs to check (meaning the total number of shortest paths is less than K) */
    //printf("Working on the %dth shortest...\n", currentK + 1);

    vector<int> &currentRxnList = currentGraph.path.rxnIds;
    vector<int> &excludedRxnIds = currentGraph.excludedRxnIds;

    /* Set all of the specified reactions to not be included */
    for(int i=0; i<excludedRxnIds.size();i++) {
      rxnspace.rxnPtrFromId(excludedRxnIds[i])->old_likelihood = rxnspace.rxnPtrFromId(excludedRxnIds[i])->current_likelihood;
      rxnspace.rxnPtrFromId(excludedRxnIds[i])->current_likelihood = -1;
    }

    /* Compute shortest paths for next iteration */
    tmp = currentGraph;
    /* Not allowed to save a reference as PRIVATE so I have to do this: */
    RXNSPACE tmpRxn = rxnspace;

    vector<GRAPHSTORE> temporaryList(currentRxnList.size(), currentGraph);

    #pragma omp parallel for shared(L, currentRxnList, temporaryList) firstprivate(tmp, onePath, tmpRxn)
    for(int i=0; i<currentRxnList.size();i++) {
      set<BADIDSTORE> dum;
      tmp.excludedRxnIds.push_back(currentRxnList[i]);
      tmpRxn.rxnPtrFromId(currentRxnList[i])->old_likelihood = tmpRxn.rxnPtrFromId(currentRxnList[i])->current_likelihood;
      tmpRxn.rxnPtrFromId(currentRxnList[i])->current_likelihood = -1;
      tmp.path = findShortestPath(tmpRxn, metspace, inputs, output, dum);
      temporaryList[i] = tmp;
      tmp.excludedRxnIds.pop_back();
      tmpRxn.rxnPtrFromId(currentRxnList[i])->current_likelihood = tmpRxn.rxnPtrFromId(currentRxnList[i])->old_likelihood;
    }

    for(int i=0; i<temporaryList.size(); i++) {
      if(!(temporaryList[i].path.outputId==-1)) {	L.push(temporaryList[i]);
      }
    }

    /* Reset excluded reactions for the next run */
    for(int i=0; i < excludedRxnIds.size();i++) {
      rxnspace.rxnPtrFromId(excludedRxnIds[i])->current_likelihood = rxnspace.rxnPtrFromId(excludedRxnIds[i])->old_likelihood;
    }

    previousGraphRxns = currentGraph.path.rxnIds;
    
  } /* Until the end...... */

}


