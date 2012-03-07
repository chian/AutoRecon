#include "Annotations.h"
#include"DataStructures.h"
#include "MyConstants.h"
#include"XML_loader.h"
#include"pathUtils.h"
#include "Printers.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <string>

using std::string;

/* value for VALUESTORE will be the probability and ID is the reaction ID */
map<string, vector<VALUESTORE> > GeneAnnotations(RXNSPACE &rxnspace, double cut, double cut2){
  map<string, vector<VALUESTORE> > gene_map; /* Map from a gene name to a list of reaction IDs */

  vector<REACTION> &small_list = rxnspace.rxns;

  for(int i=0;i<small_list.size();i++){
    sort(small_list[i].annote.begin(), small_list[i].annote.end());
    /* Remove things that are less than cut% of the top hit */
    if(small_list[i].annote.size()>0){
      double tophit = small_list[i].annote[0].probability;
      for(int j=small_list[i].annote.size()-1;j>0;j--){
	if(small_list[i].annote[j].probability < cut*tophit){
	  small_list[i].annote.pop_back();
	}
      }
    }
    /* Place annotations into gene map |  key = gene name, stoich = vector rxns ids */
    for(int j=0;j<small_list[i].annote.size();j++){
      /* Note - you can't put arrays [c strings] into STL containers... so I just made the gene names all strings to limit the number of conversions I have to do... */
      string tempS = small_list[i].annote[j].genename;
      double prob = small_list[i].annote[j].probability;
      int RxnId = small_list[i].id;
      VALUESTORE tmpVal; tmpVal.id = RxnId; tmpVal.value = prob;

      vector<VALUESTORE> tempV = gene_map[tempS]; //super slow
      tempV.push_back(tmpVal);
      gene_map[tempS] = tempV; //hopefully that now stores it in the map
    }
  }

  /* Now that the map is constructed... identify annotations within cut2% of the maximum for each gene - throw out the rest */
  map<string, vector<VALUESTORE> >::iterator it;
  for(it=gene_map.begin(); it!=gene_map.end(); it++) {
    string tmpS = it->first;
    vector<VALUESTORE> tempV = it->second;
    if(tempV.empty()) {  printf("ERROR: Somehow an element of the map got created that shouldn't have been in GeneAnnotations...\n");  assert(false); }
    /* Sorts by probability */
    sort(tempV.begin(), tempV.end());
    double maxProb = tempV[0].value;
    for(int i=0; i<tempV.size(); i++) {
      if(tempV[i].value < maxProb * cut2) {
	removeAnnotation(rxnspace, tmpS, tempV[i]);
      }
    }
  }

  if(_db.PRINTANNOTATIONS) {
    for(int i=0;i<small_list.size();i++){
      for(int j=0;j<small_list[i].annote.size();j++){
	printf("%s %s %f\n",small_list[i].name,small_list[i].annote[j].genename.c_str(),small_list[i].annote[j].probability);
      }
    }
  }
  if(_db.SAVEANNOTATIONS) {
    ANNOTATIONS_out("./outputs/reaction_annotations", small_list);
  }

  return gene_map;
}

void removeAnnotation(RXNSPACE &rxnspace, const string &genename, const VALUESTORE &val) {

  /* NOTE- for now, I don't attempt to change the likelihood scores themselves. This is intentional - 
     if I change the likelihood scores then the penalty is probably too high for these reactions, and
     we'd like to go back later and add these annotations back in if needed. */
  vector<ANNOTATION> tmpAnnotation;
  int idx = rxnspace.idxFromId(val.id);
  for(int i=0; i<rxnspace.rxns[idx].annote.size(); i++) {
    /* Note - compare to 0 means they are the same... */
    if(rxnspace.rxns[idx].annote[i].genename.compare(genename) == 0) {
      continue;
    }
    tmpAnnotation.push_back(rxnspace.rxns[idx].annote[i]);
  }
  rxnspace.rxns[idx].annote = tmpAnnotation;
  return;
}
