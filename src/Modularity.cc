#include"DataStructures.h"
#include"XML_loader.h"
#include"pathUtils.h"

/*Functions*/

void AllHardIncludes(const vector<REACTION> &biglist, vector<REACTION> &small_list){
  small_list.clear();
  for(int i=0;i<biglist.size();i++){
    if(rougheq(biglist[i].init_likelihood,-2.0) || rougheq(biglist[i].init_likelihood,-4.0)){
      small_list.push_back(biglist[i]);
    }
  }
  return;
}

void SynIncludes(PROBLEM &Model){
  MakeSynList(Model);
  makeMagicBridges(Model);
  return;
}


