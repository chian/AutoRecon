#include <cassert>
#include <cstdio>
#include <set>
#include <vector>
#include <map>
#include <queue>

#include "DataStructures.h"
#include "shortestPath.h"
#include "pathUtils.h"
#include "Printers.h"
#include "XML_loader.h"

using std::set;
using std::map;
using std::queue;
using std::vector;
using std::priority_queue;

/* Graph is encoded in metabolite.rxnsInvolved_nosec */
PATH findShortestPath(const RXNSPACE &rxnspace, const METSPACE &metspace, const METSPACE &inputs, const METABOLITE &output, set<BADIDSTORE> &badIds) {

  badIds.clear();
  //printf("entering findShortestPath\n");
  /*Idx is basically the index of a metabolite still in Q 
    (the set of non-optimal points), and value is the
    length of the shortest path yet found */
  priority_queue<VALUESTORE> nodeList;

  /* Bookkeeping of various sorts */
  vector<bool> isDoneAlready(metspace.mets.size(), false);
  vector<double> values(metspace.mets.size(), 0.0f);
  map<int, bool> metsExplored;

  /* For tracing back - maps met ID to reaction ID */
  map<int, int> precursorRxnIds;

  /* Initialize */
  for(int i=0; i<metspace.mets.size();i++) {
    VALUESTORE tmp;
    tmp.id = metspace.mets[i].id;
    if(inputs.idIn(tmp.id)) {
      /* Input */
      tmp.value = 0.0f;
      precursorRxnIds[tmp.id] = -1; /* -1 signifies the end of a pathway when we're backtracing*/
      nodeList.push(tmp);
    } else {
      tmp.value = 1000000.0f;
      precursorRxnIds[tmp.id] = -2; /* -2 signifies that the metabolite has not been reached by Dijkstras (this gets updated 
				       to the precursor reaction ID later in the algorithm - running into this in the backtracing step is an error) */
    }
    values[i] = tmp.value;
  }

  if(nodeList.size() == 0) {
    printf("WARNING: No inputs found! Will return an empty path\n");
    PATH dum;
    return dum;
  }

  int outputIdx = metspace.idxFromId(output.id);

  while(nodeList.size() > 0) {

    vector<int> reactionList;
    VALUESTORE tmp = nodeList.top();
    int tmpValIdx = metspace.idxFromId(tmp.id);
    nodeList.pop(); /* Actually remove the highest value (since top() doens't remove it) */

    if(isDoneAlready[tmpValIdx]) {
      /* We happened to find another path with a higher likelihood than a path we'd found before.
	 In such a case we normally want to continue on, however we must check and make sure that there are still
	 things left to check, otherwise there is no solution */
      if(nodeList.size() == 0) {
	//printf("No path found to output [found from duplicate] \n");
	//	printIntMap_mets(metspace, precursorRxnIds);
	PATH tmpPath;
	return tmpPath;
      }      
      continue;
    }
    isDoneAlready[tmpValIdx] = true;

    /* We have found optimal path to the specified output already */
    if(isDoneAlready[outputIdx]) {
      // printf("Optimal path found to output\n");
      break;
    }

    /* Find all products of reactions starting with the given metabolite. Update V(P) */
    reactionList = metspace.mets[tmpValIdx].rxnsInvolved_nosec;

    for(int i=0; i<reactionList.size();i++) {
      vector<int> modifiedMetIdx;
      vector<int> products;
      double reactantValue(0.0f);
      const REACTION* currentRxn = rxnspace.rxnPtrFromId(reactionList[i]);

      /* DO NOT INCLUDE flag */
      if(currentRxn->current_likelihood> -1.1 && currentRxn->current_likelihood < -0.9) { continue; }

      int tmpIdx;
      for(int j=0; j<currentRxn->stoich_part.size();j++) {
	if(currentRxn->stoich_part[j].met_id == tmp.id) {
	  tmpIdx = j; 
	  break; 
	}
      }

      /* Note - it is NOT sufficient to just let the queue do its thing, we MUST explicitly identify all of
	 the reactants as already having been reached optimally. Otherwise the code will incorrectly allow
	 just one reactant to be present before labeling the products. 

	 Try to keep track of reactions that we find blocked here so that we can go and try to un-block them later
      */
      bool notAllInputsPresent = false;
      BADIDSTORE tmpBad;
      for(int j=0; j<currentRxn->stoich_part.size(); j++) {
	if(currentRxn->stoich_part[j].rxn_coeff * currentRxn->stoich_part[tmpIdx].rxn_coeff > 0.0f) {
	  if(!isDoneAlready[metspace.idxFromId(currentRxn->stoich_part[j].met_id)]) { 
	    notAllInputsPresent = true; 
	    tmpBad.badRxnId = currentRxn->id;
	    tmpBad.badMetIds.push_back(currentRxn->stoich_part[j].met_id);
	  }
	}
      }
      if(notAllInputsPresent) {
	badIds.insert(tmpBad);
	continue;
      } else { badIds.erase(tmpBad); }


      for(int j=0; j<currentRxn->stoich_part.size();j++) {
	if(currentRxn->stoich_part[j].rxn_coeff * currentRxn->stoich_part[tmpIdx].rxn_coeff > 0) { 
	  reactantValue += values[metspace.idxFromId(currentRxn->stoich_part[j].met_id)];
	}
      }     

      /* Check for other negative likelihoods. If any exist Dijkstras will die a horrible, horrible death */
      if(currentRxn->current_likelihood < 0) {
	printf("ERROR: in Dijkstras algorithm - NEGATIVE LIKELIHOOD %1.2f for REACTION %d\n", currentRxn->current_likelihood, currentRxn->id);
	assert(currentRxn->current_likelihood >= 0);
      }

      getProducts(*currentRxn, tmp.id, products);

      for(int j=0; j<products.size();j++) {
	double newValue = reactantValue + currentRxn->current_likelihood;
	int prodIdx = metspace.idxFromId(products[j]);
        if(values[prodIdx] > newValue ) {
	  values[prodIdx] = newValue;
	  precursorRxnIds[products[j]] = reactionList[i];
	  modifiedMetIdx.push_back(prodIdx);
	}
      }

      /* Add updated values to the heap */
      for(int j=0; j<modifiedMetIdx.size();j++) {
	VALUESTORE mod;
	mod.id = metspace.mets[modifiedMetIdx[j]].id;
	mod.value = values[modifiedMetIdx[j]];
	nodeList.push(mod);
      }
    } /* For each reaction in reactionList */

    if(nodeList.size() == 0) {
      //printf("No path found to output\n");
      PATH tmpPath;
      //      printRxnsFromIntSet(badIds, rxnspace);
      return tmpPath;
    }

  } /* While nodeList.size() > 0 */

  /* Trace the path back */
  queue<int> dfsList; dfsList.push(output.id);
  vector<int> rxnIds;
  vector<int> inputIdList;
  vector<int> rxnDirections; /* -1 if used in negative direction and +1 in positive direction */
  /* Note - at the end of this, "explored" should be TRUE for anything in the path and FALSE for everything else */

  tracePath(metspace, rxnspace, precursorRxnIds, metsExplored, rxnIds, inputIdList, rxnDirections, dfsList);

  PATH tracedPath;
  /* Explored contains all of the metabolites (including inputs and outputs!) that were crossed on the way 
   Now we translate that back into a list of indicies */
  tracedPath.rxnIds = rxnIds;
  tracedPath.rxnDirection = rxnDirections;
  tracedPath.inputIds = inputIdList;
  tracedPath.outputId = output.id;
  tracedPath.totalLikelihood = values[metspace.idxFromId(output.id)];
  for(int i=0; i<metsExplored.size();i++) {
    /* FIXME: We want to include inputs but NOT outputs here? */
    if(metsExplored[i]) {
      tracedPath.metsConsumedIds.push_back(metspace.mets[i].id);
    }
  }
  /* Assign priorities - the first element in the resulting list has the highest priority, the second element has the second highest, etc... 
   Since the synthesis reaction is first in the rxnIds list and the inputs are (should be) last... we need to reverse it as seen here */
  for(int i=rxnIds.size()-1; i>=0; i--) {
    tracedPath.rxnPriority.push_back(rxnIds[i]);
  }

  /* Identify dead ends */
  vector<int> allMets;
  for(int i=0; i<tracedPath.rxnIds.size(); i++) {
    const REACTION* rxn = rxnspace.rxnPtrFromId(tracedPath.rxnIds[i]);
    for(int j=0; j < rxn->stoich_part.size(); j++) {
      allMets.push_back(rxn->stoich_part[j].met_id);
    }
  }

  vector<int> toRemove = tracedPath.metsConsumedIds;
  /* Also remove the output ID (inputs should already be in metsConsumed) */
  toRemove.push_back(output.id);
  tracedPath.deadEndIds = setdiff2(allMets, toRemove);

  return tracedPath;
}

/*    precursorRxnIds: map between metabolite IDs and reaction IDs that they came from (-1 for inputs)
      metsExplored: needed internally. true for anything that has already been traced by the DFS
      rxnIds: ID of any reactions that have been traced
      inputIds: ID of any inputs (-1) found while tracing
      rxnDirections: Direction of any reaction traced relative to its reversibility (indexed the same way as rxnIds)
      nodeList: Queue used for DFS

 Note that if we get an indexing out of bounds here that indicates I messed up the code... */

void tracePath(const METSPACE &metspace, const RXNSPACE &rxnspace, map<int, int>  &precursorRxnIds, map<int, bool> &metsExplored, vector<int> &rxnIds, vector<int> &inputIds, 
	       vector<int> &rxnDirections, queue<int> &nodeList) {

  /* Termination condition - nothing left in the queue! */
  if(nodeList.size() == 0) { 
    return; 
  }

  int currentNodeId = nodeList.front();
  nodeList.pop();
  int currentEdgeId = precursorRxnIds[currentNodeId];

  metsExplored[currentNodeId] = true;

  if(currentEdgeId == -1) { 
    inputIds.push_back(currentNodeId);
    /* This is needed because otherwise after the return it ignores the rest of the nodeList! */
    tracePath(metspace, rxnspace, precursorRxnIds, metsExplored, rxnIds, inputIds, rxnDirections, nodeList);
    return;
  }
  assert(currentEdgeId != -2);

  int sgn(0);
  const REACTION* rxn = rxnspace.rxnPtrFromId(currentEdgeId);
  for(int i=0; i<rxn->stoich_part.size();i++) {
    if(rxn->stoich_part[i].met_id == currentNodeId) {
      if(rxn->stoich_part[i].rxn_coeff < 0) { sgn = -1; } else { sgn = 1; };
      break;
    }
  }

  /* Note - sgn should now always be -1 or +1 */
  vector<int> oppositeIds = opposite(rxn, currentNodeId, sgn);

  rxnDirections.push_back(sgn);
  rxnIds.push_back(currentEdgeId);

  for(int i=0; i<oppositeIds.size(); i++) {
    /* This should prevent loops */
    if(metsExplored[oppositeIds[i]]) { 
      continue; 
    }
    nodeList.push(oppositeIds[i]);
  }

  tracePath(metspace, rxnspace, precursorRxnIds, metsExplored, rxnIds, inputIds, rxnDirections, nodeList);
}

/* Find metabolites opposite of metId in a reaction (based on stoich_part) */
/* 6-14-11: CONFIRMED working by printout */
vector<int> opposite(const REACTION* rxn, int metId, int sgn) {

    vector<int> finalList;
    for(int i=0; i<rxn->stoich_part.size();i++) {
      if(sgn*rxn->stoich_part[i].rxn_coeff < 0) {
	finalList.push_back(rxn->stoich_part[i].met_id);
      }
    }
    return finalList;
}
