#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

/* Personal libaries */
#include "Grow.h"
#include "MyConstants.h"
#include "shortestPath.h"
#include "visual01.h"
#include "XML_loader.h"

using std::vector;
using std::map;
using std::string;

/* Generates a .dot file to be viewed after a 
   "dot -Tpng dotfile.dot > pngfile.png" command */
void Paths2Dot(FILE *dotput, const METSPACE &metspace,
	       const RXNSPACE &rxnspace, const char *label){

  const vector<METABOLITE> &incl_mets = metspace.mets;
  const vector<REACTION> &incl_rxns = rxnspace.rxns;
  int numMETs = incl_mets.size();
  int numRXNs = incl_rxns.size();

  /* Header for DOT file */
  fprintf(dotput,"digraph pathway_plotter { \n");
  fprintf(dotput,"sep = 1 ;\n");
  fprintf(dotput,"size = \"8,11!\" ;\n");
  fprintf(dotput,"ranksep = 1 ;\n");
  fprintf(dotput,"rankdir = LR ;\n");
  fprintf(dotput,"mincross = 2.0 ;\n");
  fprintf(dotput,"label = \"%s\" ; \n",label);

  /* MET nodes - creates all MET nodes involved in reactions */
  for(int i=0;i<numMETs;i++){
    if(!incl_mets[i].rxnsInvolved_nosec.empty()){
      if(incl_mets[i].input!=1 && incl_mets[i].output!=1){

	if(incl_mets[i].secondary_pair.empty() && incl_mets[i].secondary_lone == 0) {
	  fprintf(dotput,"\"%s\" [shape=ellipse,margin=0,regular=1,style=filled,fillcolor=orange1,fontsize=12,height=0,width=0] ; \n",incl_mets[i].name);
	}

	if(incl_mets[i].secondary_lone != 0) { continue; }

	/* Only make nodes for secondary pair metabolites if they are part of magic bridge fixes...*/
	bool printName = false;
	/*	for(int j=0; j<incl_rxns.size(); j++) {
	  	  for(int k=0; k<incl_rxns[j].magicBridgeFix.size(); k++) {
	    if(incl_rxns[j].magicBridgeFix[k]==incl_mets[i].id) { 
	      printName = true; break; }
	  }
	  if(printName) { break; }
	}
	if(printName) { */
	  fprintf(dotput,"\"%s\" [shape=ellipse,margin=0,regular=1,style=filled,fillcolor=orange1,fontsize=12,height=0,width=0] ; \n",incl_mets[i].name);
	  /*	} */
      }
    }
  }

  /* Double circles for INPUTS & OUTPUTS */
  for(int i=0;i<numMETs;i++){
    if(incl_mets[i].input==1 || incl_mets[i].output==1){
      fprintf(dotput,"node [shape=doublecircle,label=\"%s\",margin=0,regular=1,style=filled,fillcolor=orange1,fontsize=12,height=0,width=0] \"%s\" ; \n",incl_mets[i].name,incl_mets[i].name);}}

  /* Double circles around BIOMASS */
  int bm = -9;
  for(int i=0;i<numRXNs;i++){
    if(incl_rxns[i].id==_db.BIOMASS){
      fprintf(dotput,"node [shape=doublecircle,label=\"%s\",margin=0,regular=1,style=filled,fillcolor=orange1,fontsize=12,height=0,width=0] \"%s\" ; \n",incl_rxns[i].name,incl_rxns[i].name);
      bm = i; break;}}

  /* BLANK metabolites - used to set up the dangling secondaries that make things MUCH easier to read. */
  /* FIXME: Need to specially deal here with fixes for magic bridges because we don't want to
     separate cofactor pairs from those */
  map<int, int> metId2Count;
  for(int i=0; i<numRXNs; i++) {
    bool pairsDone = false;
    for(int j=0; j<incl_rxns[i].stoich_part.size(); j++) {
      METABOLITE MET = metspace.metFromId(incl_rxns[i].stoich_part[j].met_id);
      /* Keep track of how many blanks at the end of the secondary edges we need to keep track of */
      if(MET.secondary_lone == 1 | hasValidSecondaryPair(incl_rxns[i], metspace, MET.id)) {
	// Secondary lones always will be "dangling" off by definition
	if(metId2Count.find(MET.id)==metId2Count.end()) { metId2Count[MET.id] = 0; }
	metId2Count[MET.id] += 1;
	int count = metId2Count[MET.id];
	fprintf(dotput,"\"%s_BLANK%d\" [shape=ellipse,label=\"\",margin=0,regular=1,style=filled,fillcolor=white,fontsize=6,height=0,width=0] ; \n",MET.name,count);
      }
    }
  }

  int good = 0;
  /* RXN boxes */
  for(int i=0;i<numRXNs;i++){
    int k=0; int l=0;
    for(int j=0;j<incl_rxns[i].stoich_part.size();j++){
      if(incl_rxns[i].stoich_part[j].rxn_coeff>0){k++;} 
      if(incl_rxns[i].stoich_part[j].rxn_coeff<0){l++;}
    }
    if(l>0 && k==0){
      fprintf(dotput,"subgraph cluster_%s { \n",incl_rxns[i].name);
      fprintf(dotput,"label = \"%s\" ; \n",incl_rxns[i].name);
      fprintf(dotput,"margin = 0 ; \n");
      fprintf(dotput,"fontsize = 12 ; \n");
      //fprintf(dotput,"color = blue ; \n");
      fprintf(dotput,"node [label=\"\",margin=0,regular=1,style=filled,fillcolor=white,fontsize=6,height=0,width=0] \"%s%s\" ; \n",incl_rxns[i].name,"_in");
      fprintf(dotput,"node [label=\"\",margin=0,regular=1,style=filled,fillcolor=white,fontsize=6,height=0,width=0] \"%s%s\" ; \n",incl_rxns[i].name,"_out");
      fprintf(dotput,"} \n");
      good = 1;
    }
    if(k>0 && l==0){
      fprintf(dotput,"subgraph cluster_%s { \n",incl_rxns[i].name);
      fprintf(dotput,"label = \"%s\" ; \n",incl_rxns[i].name);
      fprintf(dotput,"margin = 0 ; \n");
      fprintf(dotput,"fontsize = 12 ; \n");
      //fprintf(dotput,"color = blue ; \n");
      fprintf(dotput,"node [label=\"\",margin=0,regular=1,style=filled,fillcolor=white,fontsize=6,height=0,width=0] \"%s%s\" ; \n",incl_rxns[i].name,"_in");
      fprintf(dotput,"node [label=\"\",margin=0,regular=1,style=filled,fillcolor=white,fontsize=6,height=0,width=0] \"%s%s\" ; \n",incl_rxns[i].name,"_out");
      fprintf(dotput,"} \n");
      //fprintf(dotput,"\"%s%s\" -> \"%s%s\" [dir=none,label=\"\",margin=0,regular=1,fontsize=8] ; \n",incl_rxns[i].name,"_in",incl_rxns[i].name,"_out");
      good = 1;
    }
    if(l>0 && k>0){
      fprintf(dotput,"subgraph cluster_%s { \n",incl_rxns[i].name);
      fprintf(dotput,"label = \"%s\" ; \n",incl_rxns[i].name);
      fprintf(dotput,"margin = 0 ; \n");
      fprintf(dotput,"fontsize = 12 ; \n");
      //fprintf(dotput,"color = blue ; \n");
      fprintf(dotput,"node [label=\"\",margin=0,regular=1,style=filled,fillcolor=white,fontsize=6,height=0,width=0] \"%s%s\" ; \n",incl_rxns[i].name,"_in");
      fprintf(dotput,"node [label=\"\",margin=0,regular=1,style=filled,fillcolor=white,fontsize=6,height=0,width=0] \"%s%s\" ; \n",incl_rxns[i].name,"_out");
      fprintf(dotput,"} \n");
      //fprintf(dotput,"\"%s%s\" -> \"%s%s\" [dir=none,label=\"\",margin=0,regular=1,fontsize=8] ; \n",incl_rxns[i].name,"_in",incl_rxns[i].name,"_out");
      good = 1;
    }
    if(good == 1) {
      if(incl_rxns[i].net_reversible==0) {
	      fprintf(dotput,"\"%s%s\" -> \"%s%s\" [dir=both,label=\"\",margin=0,regular=1,fontsize=8] ; \n",incl_rxns[i].name,"_in",incl_rxns[i].name,"_out"); }
      else if(incl_rxns[i].net_reversible==1) {
	      fprintf(dotput,"\"%s%s\" -> \"%s%s\" [dir=forward,label=\"\",margin=0,regular=1,fontsize=8] ; \n",incl_rxns[i].name,"_in",incl_rxns[i].name,"_out");
      }
      else if(incl_rxns[i].net_reversible==-1) {
	      fprintf(dotput,"\"%s%s\" -> \"%s%s\" [dir=back,label=\"\",margin=0,regular=1,fontsize=8] ; \n",incl_rxns[i].name,"_in",incl_rxns[i].name,"_out");
      }
    }
  }

  /* Links between METs and RXNs */

  /* Biomass special */
  if(bm >= 0) {
    for(int i=0;i<incl_rxns[bm].stoich.size();i++){
      fprintf(dotput,"\"%s\" -> \"%s\" [weight=1] ; \n",
	      incl_rxns[bm].stoich[i].met_name,incl_rxns[bm].name);}
  }

  /* Normal Reactions */
  for(int i=0;i<numRXNs;i++){
    int k=0; int l=0;
    for(int j=0;j<incl_rxns[i].stoich_part.size();j++){
      if(incl_rxns[i].stoich_part[j].rxn_coeff>0){k++;} 
      if(incl_rxns[i].stoich_part[j].rxn_coeff<0){l++;}
    }
    if(k>0 && l>0){
      for(int j=0;j<incl_rxns[i].stoich_part.size();j++){
	if(incl_rxns[i].net_reversible==0 
	   && incl_rxns[i].stoich_part[j].rxn_coeff<0){ /* reversible arrows */
	  METABOLITE MET = metspace.metFromId(incl_rxns[i].stoich_part[j].met_id);
	  if(MET.secondary_lone == 1) {
	    int count = metId2Count[MET.id];
	    metId2Count[MET.id] -= 1;
	    fprintf(dotput,"\"%s_BLANK%d\" -> \"%s%s\" [dir=none,weight=1,label=\"%s\",fontsize=14] ; \n",
		    incl_rxns[i].stoich_part[j].met_name, count,
		    incl_rxns[i].name,"_in",
		    incl_rxns[i].stoich_part[j].met_name);
	  } else if(hasValidSecondaryPair(incl_rxns[i], metspace, MET.id)) {
	    /* We DON'T need to print the secondaries here since we're looping through the whole thing anyway
	       and we'll hit the other half of the pair on the way around (and the pairs are labeled both ways) */
	    int count = metId2Count[MET.id];
	    metId2Count[MET.id] -= 1;
	    fprintf(dotput,"\"%s_BLANK%d\" -> \"%s%s\" [dir=none,weight=1,label=\"%s\",fontsize=14] ; \n",
		    incl_rxns[i].stoich_part[j].met_name, count,
		    incl_rxns[i].name,"_in",
		    incl_rxns[i].stoich_part[j].met_name);
	  } else {
	    fprintf(dotput,"\"%s\" -> \"%s%s\" [dir=none,weight=1] ; \n",
		    incl_rxns[i].stoich_part[j].met_name,
		    incl_rxns[i].name,"_in");
	  }
	}
	if(incl_rxns[i].net_reversible==0 
	   && incl_rxns[i].stoich_part[j].rxn_coeff>0){ /* reversible arrows */
	  METABOLITE MET = metspace.metFromId(incl_rxns[i].stoich_part[j].met_id);
	  if(MET.secondary_lone == 1) {
	    int count = metId2Count[MET.id];
	    metId2Count[MET.id] -= 1;
	    fprintf(dotput,"\"%s%s\" -> \"%s_BLANK%d\" [dir=none,weight=1,label=\"%s\",fontsize=14] ; \n",
		    incl_rxns[i].name, "_out", 
		    incl_rxns[i].stoich_part[j].met_name, count,
		    incl_rxns[i].stoich_part[j].met_name);
	  } else if(hasValidSecondaryPair(incl_rxns[i], metspace, MET.id)) {
	    int count = metId2Count[MET.id];
	    metId2Count[MET.id] -= 1;
	    fprintf(dotput,"\"%s%s\" -> \"%s_BLANK%d\" [dir=none,weight=1,label=\"%s\",fontsize=14] ; \n",
		    incl_rxns[i].name, "_out", 
		    incl_rxns[i].stoich_part[j].met_name, count,
		    incl_rxns[i].stoich_part[j].met_name); 
	  } else {
	    fprintf(dotput,"\"%s%s\" -> \"%s\" [dir=none,weight=1] ; \n",
		    incl_rxns[i].name, "_out", incl_rxns[i].stoich_part[j].met_name);
	  }
	}
	if((incl_rxns[i].stoich_part[j].rxn_coeff<0 
	    && incl_rxns[i].net_reversible==1)
	   ||(incl_rxns[i].stoich_part[j].rxn_coeff>0 
	      && incl_rxns[i].net_reversible==-1)){ /* into rxn_in */

	  METABOLITE MET = metspace.metFromId(incl_rxns[i].stoich_part[j].met_id);
	  if(MET.secondary_lone == 1) {
	    int count = metId2Count[MET.id];
	    metId2Count[MET.id] -= 1;
	    fprintf(dotput,"\"%s_BLANK%d\" -> \"%s%s\" [dir=none,weight=1,label=\"%s\",fontsize=14] ; \n",
		    incl_rxns[i].stoich_part[j].met_name, count,
		    incl_rxns[i].name,"_in",
		    incl_rxns[i].stoich_part[j].met_name);
	  } else if(hasValidSecondaryPair(incl_rxns[i], metspace, MET.id)) {
	    /* We DON'T need to print the secondaries here since we're looping through the whole thing anyway
	       and we'll hit the other half of the pair on the way around (and the pairs are labeled both ways) */
	    int count = metId2Count[MET.id];
	    metId2Count[MET.id] -= 1;
	    fprintf(dotput,"\"%s_BLANK%d\" -> \"%s%s\" [dir=none,weight=1,label=\"%s\",fontsize=14] ; \n",
		    incl_rxns[i].stoich_part[j].met_name, count,
		    incl_rxns[i].name,"_in",
		    incl_rxns[i].stoich_part[j].met_name);
	  } else {
	    // Not a secondary metabolite
	    fprintf(dotput,"\"%s\" -> \"%s%s\" [dir=none,weight=1] ; \n",
		    incl_rxns[i].stoich_part[j].met_name,
		    incl_rxns[i].name,"_in");
	  }
	}
	
	if((incl_rxns[i].stoich_part[j].rxn_coeff<0 
	    && incl_rxns[i].net_reversible==-1)
	   ||(incl_rxns[i].stoich_part[j].rxn_coeff>0 
	      && incl_rxns[i].net_reversible==1)){ /* out of rxn_out */

	  METABOLITE MET = metspace.metFromId(incl_rxns[i].stoich_part[j].met_id);
	  if(MET.secondary_lone == 1) {
	    int count = metId2Count[MET.id];
	    metId2Count[MET.id] -= 1;
	    fprintf(dotput,"\"%s%s\" -> \"%s_BLANK%d\" [dir=none,weight=1,label=\"%s\",fontsize=14] ; \n",
		    incl_rxns[i].name, "_out", 
		    incl_rxns[i].stoich_part[j].met_name, count,
		    incl_rxns[i].stoich_part[j].met_name); 
	  } else if(hasValidSecondaryPair(incl_rxns[i], metspace, MET.id)) {
	    int count = metId2Count[MET.id];
	    metId2Count[MET.id] -= 1;
	    fprintf(dotput,"\"%s%s\" -> \"%s_BLANK%d\" [dir=none,weight=1,label=\"%s\",fontsize=14] ; \n",
		    incl_rxns[i].name, "_out", 
		    incl_rxns[i].stoich_part[j].met_name, count,
		    incl_rxns[i].stoich_part[j].met_name); 
	  } else {
	    // Not a secondary metabolite
	    fprintf(dotput,"\"%s%s\" -> \"%s\" [dir=none,weight=1] ; \n",
		    incl_rxns[i].name, "_out", 
		    incl_rxns[i].stoich_part[j].met_name);
	  }
	}
      }
    }
  }

  fprintf(dotput,"}\n");

  return;
}

/* Test if reaction RXN has metablite metId and its opposite, it is not a magic bridge fix, 
   AND it is not a magic bridge. This indicates ones that we want to REMOVE */
bool hasValidSecondaryPair(const REACTION &rxn, const METSPACE &workingMets, int met_id) {

  if(rxn.id >= _db.MAGICBRIDGEFACTOR && rxn.id < _db.MAGICBRIDGEFACTOR + _db.MINFACTORSPACING) { return false; }

  /* Check if met_id is present in the magicBridgeFix */
  /*  for(int i=0; i<rxn.magicBridgeFix.size(); i++) {
    if(rxn.magicBridgeFix[i] == met_id) { return false; }
    } */

  METABOLITE met = workingMets.metFromId(met_id);

  /* Is met_id in rxn? */
  bool good = false;
  for(int i=0; i<rxn.stoich_part.size(); i++) {
    if(met.id == rxn.stoich_part[i].met_id) { good = true; break; }
  }
  if(!good) { return false; }

  /* Is at least ONE of the secondary_pairs of met_id in the rxn? */
  good = false;
  for(int i=0; i<met.secondary_pair.size(); i++) {
    for(int j=0; j<rxn.stoich_part.size(); j++) {
      if(rxn.stoich_part[j].met_id == met.secondary_pair[i]) { good = true; break; }
    }
    if(good) { break; }
  }

  if(good) { return true; } else { return false; }

}

/* Take a PATH structure and convert it to inputs needed for the dot converter */
/* useSyns = true --> use synRxns
   useSyns = false --> use fullRxns */
void convertPathForVisuals(const PATH &path, const PROBLEM &ProblemSpace,  RXNSPACE &workingRxns, METSPACE &workingMets, int useSyns) {
  switch(useSyns){
  case 0:
    pullOutRxnsbyIds(ProblemSpace.fullrxns, path.rxnIds, workingRxns);
    pullOutMets(workingRxns, ProblemSpace.metabolites, workingMets);
    /* For fullrxns it's possible someone forgot to copy stoich over to stoich_part.
       Do that for them here. 

       This basically assumes that if you call the function with FULLRXNS, then you intend to graph with the cofactors present, 
    since it doesn't really make much sense otherwise (if you use FULLRXNS without cofactors, 
    you'll have a lot of rxns that look identical taking up space) */
    for(int i=0; i<workingRxns.rxns.size(); i++) {
      workingRxns.rxns[i].stoich_part = workingRxns.rxns[i].stoich;
    }
    break;
  case 1: 
    pullOutRxnsbyIds(ProblemSpace.synrxns, path.rxnIds, workingRxns);
    pullOutMets(workingRxns, ProblemSpace.metabolites, workingMets);
    break;
  case 2:
    pullOutRxnsbyIds(ProblemSpace.synrxnsR, path.rxnIds, workingRxns);
    pullOutMets(workingRxns, ProblemSpace.metabolites, workingMets);
    break;
  }
    
  return;
}

/* Pass the rxnspace so that we can distinguish synrxns from fullrxns 
 Note that either way we are using STOICH_PART and not STOICH. 

useSynRxns = true --> pull out reactions from ProblemSpace.synrxns
useSynRxns = false --> pull out reactions from ProblemSpace.fullrxns 

Which you use depends on where the PATH came from. */
void VisualizePath2File(const char *fileName, const char *label, const PATH &path, 
			const PROBLEM &ProblemSpace, int useSynRxns){
  //use SynRxns, 0 = fullrxns, 1 = synrxns, 2 = synrxnsR
  RXNSPACE workingRxns;
  METSPACE workingMets;

  FILE *dotput = fopen(fileName, "w");
  if(dotput==NULL) { printf("Unable to open DOT file in VisualizePath2File \n");fflush(stdout);
    return;}

  convertPathForVisuals(path, ProblemSpace, workingRxns, workingMets, useSynRxns);
  Paths2Dot(dotput, workingMets, workingRxns, label);
  fclose(dotput);

  printf("Exiting VisualizePath2File\n");
  return;
}

/* Visualize the fullrxns

   useSyn: true - use synRxns (must already have stoich_part filled)
          false - use fullrxns (assumes that you want to use stoich, and copies that over to stoich_part before calling the paths2Dot file)

*/
void visualizePathSummary2File(const char* fileBase, const char* label, const vector<PATHSUMMARY> &psum, const PROBLEM &ProblemSpace, int useSyn) {
  
  for(int i=0; i<psum.size(); i++) {
    char fileName[128]; sprintf(fileName, "%s%s.dot", fileBase, ProblemSpace.metabolites.metFromId(psum[i].outputId).name);
    FILE* dotput = fopen(fileName, "w");

    vector<PATHSUMMARY> tmpP; tmpP.push_back(psum[i]);

    /* Get a list of metabolites from reactions in the pList (needed later...) */
    vector<int> allRxnIds = getAllPathRxns(tmpP);
    for(int i=0; i<allRxnIds.size(); i++) { allRxnIds[i] = abs(allRxnIds[i]); }
    
    RXNSPACE modelRxns; 
    switch(useSyn) {
    case 0:
      pullOutRxnsbyIds(ProblemSpace.fullrxns, allRxnIds, modelRxns);
      // Copy stoich to stoich_part per assumption that you actually want cofactors if you're using fullrxns
      for(int i=0; i<modelRxns.rxns.size(); i++) {
	modelRxns.rxns[i].stoich_part = modelRxns.rxns[i].stoich;}
      break;
    case 1:
      pullOutRxnsbyIds(ProblemSpace.synrxns, allRxnIds, modelRxns);
      // stoich_part should already be in synrxns so no need to copy it over
      break;
    case 2: 
      pullOutRxnsbyIds(ProblemSpace.synrxnsR, allRxnIds, modelRxns);
      // stoich_part should already be in synrxns so no need to copy it over
      break;
    }

    METSPACE modelMets;  pullOutMets(modelRxns, ProblemSpace.metabolites, modelMets);

    Paths2Dot(dotput, modelMets, modelRxns, label);
    fclose(dotput);
  }
}

void VisualizeProblem(PROBLEM &small_model, const int i){
  char s[128];
  sprintf(s, "flux_o%d.dot", i);
  FILE *dotput = fopen(s,"w");
  Paths2Dot(dotput,small_model.metabolites.mets,small_model.fullrxns.rxns,s);
}

void VisualizeAllSolutions(const PROBLEM &ProblemSpace, const vector< vector<PATH> > &paths, int useSyns){
  for(int i=0;i<paths.size();i++){
    for(int j=0;j<paths[i].size();j++){
      char s[128];
      sprintf(s, "path_o%d_k%d.dot", i, j);
      VisualizePath2File(s,s,paths[i][j],ProblemSpace, useSyns);
    }
  }
  return;
}
