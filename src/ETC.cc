#include "Exchanges.h"
#include "ETC.h"
#include "Grow.h"
#include "MyConstants.h"
#include "pathUtils.h"
#include "Paths2Model.h"
#include "Printers.h"
#include "RunK.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <map>
#include <set>
#include <vector>

using std::vector;
using std::map;
using std::set;

/* Electron Transport Chain Functions */
/* Ids cannot include zero */
void decompose(const vector<int> &rxnDirIds, vector<int> &rxnIds,  vector<int> &rxnDirs){
  int i;
  rxnIds.clear();
  rxnDirs.clear();
  for(i=0;i<rxnDirIds.size();i++){
    rxnIds.push_back(abs(rxnDirIds[i]));
    if(rxnDirIds[i]>0){rxnDirs.push_back(1);}
    else{rxnDirs.push_back(-1);}
  }
  return;
}

void flip(vector<int> &rxnDirs){
  for(int i=0;i<rxnDirs.size();i++){
    rxnDirs[i] *= -1;}
  return;
}

NETREACTION flip(const NETREACTION &one){
  int i;
  NETREACTION tempN;
  vector<int> tempI;
  tempI = one.rxnDirIds;
  tempN.rxnDirIds = tempI;
  tempN.rxn = one.rxn;
  for(i=0;i<one.rxn.stoich_part.size();i++){
    tempN.rxn.stoich_part[i].rxn_coeff *= -1;
  }
  return tempN;
}

/* Identify Electrion Transport Chain Reactions*/
void identifyETCrxns(const PROBLEM &ProblemSpace, vector<int> &ETCrxnIds){
  ETCrxnIds.clear();
  int k;
  int HE_ids,Na_ids,NaE_ids;
  int HE_coe,Na_coe,NaE_coe;
  const char *MetNAME;
  char tempS[4] = {0};
  int MetID;

  const vector<METABOLITE> &metabolite = ProblemSpace.metabolites.mets;
  const vector<REACTION> &reaction = ProblemSpace.fullrxns.rxns;

  HE_ids  = Name2Ids(metabolite,_db.H_plus_E);
  Na_ids  = Name2Ids(metabolite,_db.Na_name);
  NaE_ids = Name2Ids(metabolite,_db.Na_plus_E);

  for(int i=0;i<reaction.size();i++){
    if(reaction[i].transporter==1){
      HE_coe = 0;
      Na_coe = 0;
      NaE_coe = 0;
      k = 0;
      for(int j=0;j<reaction[i].stoich.size();j++){
	int MetID   = reaction[i].stoich[j].met_id;
	const char* MetNAME = reaction[i].stoich[j].met_name;
	if((int)strlen(MetNAME)>=3){
	  strncpy(tempS,MetNAME+((int)strlen(MetNAME)-3),3);	
	  if(strcmp(tempS,_db.E_tag)==0){ k++;}
	}
	if(HE_ids==MetID){
	  HE_coe = 1;
	  k--;}
	if(NaE_ids==MetID){
	  NaE_coe = 1;}
	if(Na_ids==MetID){
	  Na_coe = 1;}
      }
      if(HE_coe==1 && k==0){
	ETCrxnIds.push_back(reaction[i].id);
      }
      else if(Na_coe==1 && NaE_coe==1 && k==0){
	ETCrxnIds.push_back(reaction[i].id);
      }
    }
  }

  if(_db.DEBUGETC){printf("identifyETCrxns result: ");printRxnsFromIntVector(ETCrxnIds,ProblemSpace.fullrxns);}
  return;
}

void identifyETCrxns(const PROBLEM &ProblemSpace, vector<int> &ETCrxnIds, double cutoff){
  vector<int> ETCrxnIdsTemp;

  const RXNSPACE &rxnspace = ProblemSpace.fullrxns;

  identifyETCrxns(ProblemSpace, ETCrxnIdsTemp);
  for(int i=0;i<ETCrxnIdsTemp.size();i++){
    double currentLikelihood = rxnspace.rxnFromId(ETCrxnIdsTemp[i]).init_likelihood;
    if(currentLikelihood > cutoff && currentLikelihood < 1.0f) {
      ETCrxnIds.push_back(ETCrxnIdsTemp[i]);
    }
  }
  return;
}

void convert(NETREACTION &net, REACTION add, int rxnDirId){
  net.rxn = add;
  net.rxnDirIds.clear();
  net.rxnDirIds.push_back(rxnDirId);
}

/* Lists cofactor pairs present in fromnet. PairOfIds[i] is the i'th cofactor pair present there */
vector<vector<int> > PossiblePairs(const vector<STOICH> &fromnet,  const METSPACE &cofactors){
  vector<vector<int> > PairOfIds;
  for(int i=0;i<fromnet.size();i++){
    if(cofactors.idIn(fromnet[i].met_id)){
      for(int j=i+1;j<fromnet.size();j++){
    	if(cofactors.metFromId(fromnet[i].met_id).secondary_pair[0]==fromnet[j].met_id){
	  if(rougheq(fromnet[i].rxn_coeff,-fromnet[j].rxn_coeff)){
	    vector<int> tempI;
	    tempI.push_back(fromnet[i].met_id);
	    tempI.push_back(fromnet[j].met_id);
	    PairOfIds.push_back(tempI);
	  }
	}
      }
    }
  }
  return PairOfIds;
}
  
/* Adds a single reaction with ID rxnDirId to the NETREACTION if there is an overlapping cofactor pair with
   compatible production / consumption ratios
   Return values: 0 if not pairable, 1 is successful pair, 2 if chain terminate, and -1 if pair has unequal ratios */
int netrxn(NETREACTION &net, const PROBLEM &ProblemSpace, int rxnDirId){

  const RXNSPACE &BaseRxns = ProblemSpace.fullrxns;
  const METSPACE &BaseMets = ProblemSpace.metabolites;

  REACTION add = BaseRxns.rxnFromId(abs(rxnDirId));

  /* Eliminate External Non-players */
  vector<STOICH> fromnet = net.rxn.stoich;
  vector<STOICH> fromadd = add.stoich;

  /* There can't be a cofactor pair if there are less than two metabolites in the reactions so say so */
  if(fromadd.size()<2){return 0;}
  if(fromnet.size()<2){return 2;}

  /* Lists of cofactor pairs in the NETREACTION and the query reaction */
  vector<vector<int> > PairsNet = PossiblePairs(fromnet, ProblemSpace.cofactors);
  vector<vector<int> > PairsAdd = PossiblePairs(fromadd, ProblemSpace.cofactors);

  /* Find overlapping pairs of cofactors between the NETRACTION and the query reaction */
  vector<vector<int> > PairsOverlap;
  for(int i=0;i<PairsNet.size();i++){
    for(int j=0;j<PairsAdd.size();j++){
      for(int k=0;k<2;k++){
	for(int l=0;l<2;l++){
	  if(PairsNet[i][k] == PairsAdd[j][l]){
	    vector<int> tempI;
	    tempI.push_back(i);
	    PairsOverlap.push_back(tempI);
  }}}}}

  /* Not allowed to have more than one overlapping pair */
  if(PairsOverlap.size()>2){
    printf("Exception Caught in netrxn - multiple overlapping pairs found !\n");
    return -1;
  }
  if(PairsOverlap.size()==0){
    return 0;
  }

  /* Identify location of the two cofactors in the NETREACTION (m1 / n1) and
     the query reaction (m2 / n2) */
  int metId1 = PairsNet[PairsOverlap[0][0]][0];
  int metId2 = PairsNet[PairsOverlap[0][0]][1];
  int m1,m2,n1,n2;

  for(int i=0;i<fromnet.size();i++){
    if(fromnet[i].met_id==metId1){m1=i;}
    if(fromnet[i].met_id==metId2){n1=i;}
  }
  for(int i=0;i<fromadd.size();i++){
    if(fromadd[i].met_id==metId1){m2=i;}
    if(fromadd[i].met_id==metId2){n2=i;}
  }

  /* Check for compatible ratios */
  double ratio1 = -fromnet[m1].rxn_coeff / fromadd[m2].rxn_coeff;
  double ratio2 = -fromnet[n1].rxn_coeff / fromadd[n2].rxn_coeff;

  if(!rougheq(ratio1,ratio2)){
    printf("ERROR: In netrxn: reactions found with unequal ratios: ratios are %f and %f\n",ratio1,ratio2); 
    return -1;}

  for(int i=0;i<fromadd.size();i++){ fromadd[i].rxn_coeff *= ratio1; }

  for(int i=0;i<fromadd.size();i++){
    int k=1;
    for(int j=0;j<fromnet.size();j++){
      if(fromadd[i].met_id==fromnet[j].met_id){
	fromnet[j].rxn_coeff += fromadd[i].rxn_coeff;
	k=0;
	break;
      }	
    }
    if(k){fromnet.push_back(fromadd[i]);}
  }

  net.rxn.stoich.clear();

  for(int j=0;j<fromnet.size();j++){
    if(!rougheq(fromnet[j].rxn_coeff,0.0f)){
      net.rxn.stoich.push_back(fromnet[j]); 
    }
  }

  if(ratio1>0) {  net.rxnDirIds.push_back(rxnDirId);
  }  else if(ratio1<0){  net.rxnDirIds.push_back(-rxnDirId); }

  return 1;
}

/* Recursive function that will add to net reaction 
DO NOT make net a pass by reference */
NETREACTION ETC_add(const PROBLEM &ProblemSpace, NETREACTION net, vector<NETREACTION> &bigout){

  NETREACTION net_out;

  const RXNSPACE &rxnspace = ProblemSpace.fullrxns;
  const METSPACE &metspace = ProblemSpace.metabolites;
  const METSPACE &cofspace = ProblemSpace.cofactors;

  /* Find Candiates */
  if(_db.DEBUGETC){
    printf("ETC_add: In Rxns: ");
    printRxnsFromIntVector(net.rxnDirIds,rxnspace);fflush(stdout);}

  vector<STOICH> fromnet = net.rxn.stoich;
  vector<vector<int> > PairsNet = PossiblePairs(fromnet,cofspace);

  /* Terminate if we have reached the point where no more reactions are connected */
  if(PairsNet.size() < 1){
    if(bigout.size() > 0){
      if(bigout.back().rxnDirIds.size()!=net.rxnDirIds.size()){
	bigout.push_back(net);}
      else{
	int j = 0;
	for(int i=0;i<net.rxnDirIds.size();i++){
	  if(abs(net.rxnDirIds[i])==abs(bigout.back().rxnDirIds[i])){j++;}
	}
	if(j!=net.rxnDirIds.size()){ bigout.push_back(net);}
      }
    }
    else{  bigout.push_back(net);  }
    return net;
  }

  if(_db.DEBUGETC){printf("ETC_add: pass check... continue\n");}

  /* For each pair of cofactors found in the current NETREACTION... */
  for(int i=0;i<PairsNet.size();i++){
    /* Identify all reactions connected to the cofactors in the current cofactor pair 
     and then those that have both (i.e. those that have the whole pair) */
    vector<int> rxnSet1 = cofspace.metFromId(PairsNet[i][0]).rxnsInvolved_nosec;
    vector<int> rxnSet2 = cofspace.metFromId(PairsNet[i][1]).rxnsInvolved_nosec;
    vector<int> rxnSet=custom_intersect(rxnSet1,rxnSet2);
    
    vector<int> tempI;
    for(int j=0;j<net.rxnDirIds.size();j++){  tempI.push_back(abs(net.rxnDirIds[j])); }
    setdiff1(rxnSet,tempI);

    for(int j=0;j<rxnSet.size();j++){
      net_out = net;
      /* Try to add rxnSet[j] to the net reaction */
      int status = netrxn(net_out,ProblemSpace,rxnSet[j]);
      if(status==1){
	net_out = ETC_add(ProblemSpace,net_out,bigout);
      }
    }
  }

  /* Should be net_out? */
  return net;

}

/* Finds all possible ETCs (without regard to the actual PATHs) - we will try to connect them
   to the results of K-shortest paths later */
void ETC(const PROBLEM &ProblemSpace, vector<NETREACTION> &bigout, double cutoff){
  NETREACTION net;
  vector<int> ETC_start_ids;

  const METSPACE &BaseMets = ProblemSpace.metabolites;
  const RXNSPACE &BaseRxns = ProblemSpace.fullrxns;

  /* Populate list of starting points */
  identifyETCrxns(ProblemSpace,ETC_start_ids,cutoff);
  /* Start searching down chains */
  for(int i=0;i<ETC_start_ids.size();i++){
    convert(net,BaseRxns.rxnFromId(ETC_start_ids[i]),ETC_start_ids[i]);
    net = ETC_add(ProblemSpace,net,bigout);
  }

  /* Unique Chains - problem: should we do this with stoich or stoich_part? */
  for(int i=0;i<bigout.size();i++){
    sort(bigout[i].rxnDirIds.begin(),bigout[i].rxnDirIds.end());
    sort(bigout[i].rxn.stoich_part.begin(),bigout[i].rxn.stoich_part.end());
  }

  for(int i=0;i<bigout.size();i++){
    NETREACTION element = bigout[i];
    for(int j=i+1;j<bigout.size();j++){
      if(element==bigout[j]){
	bigout.erase(bigout.begin()+j);
      }
    }
  }

  vector<NETREACTION> bigout2 = ETC_dir_check(ProblemSpace.fullrxns, bigout);

  if(_db.PRINTETCRESULTS) { printNetReactionVector(bigout2, ProblemSpace); }

  bigout = bigout2;

  return;
}

vector<NETREACTION> ETC_dir_check(const RXNSPACE &rxnspace, const vector<NETREACTION> &bigout){
  vector<NETREACTION> bigout2;

  for(int i=0;i<bigout.size();i++){
    /* check if it will work forwards */
    vector<int> rxnIds,rxnDirs;
    decompose(bigout[i].rxnDirIds,rxnIds,rxnDirs);

    int numForwardGood = 0;
    for(int j=0;j<bigout[i].rxnDirIds.size();j++){
      //printf("ETC_dir_check: %d %d %d\n",i,k,bigout[i].rxn.net_reversible);
      if(rxnspace.rxnFromId(rxnIds[j]).net_reversible==0){numForwardGood++;}
      if(rxnspace.rxnFromId(rxnIds[j]).net_reversible==rxnDirs[j]){numForwardGood++;}
    }
    if(numForwardGood==bigout[i].rxnDirIds.size()){
      bigout2.push_back(bigout[i]);
      bigout2.back().rxn.net_reversible = 1;
      bigout2.back().rxn.lb = 0.0f;
      bigout2.back().rxn.ub = 1000.0f;
    }
    /* check if it will work backwards */
    flip(rxnDirs);
    int numRevGood = 0;
    for(int j=0;j<bigout[i].rxnDirIds.size();j++){
      if(rxnspace.rxnFromId(rxnIds[j]).net_reversible==0){numRevGood++;}
      if(rxnspace.rxnFromId(rxnIds[j]).net_reversible==rxnDirs[j]){numRevGood++;}
    }
    if(numRevGood==bigout[i].rxnDirIds.size()){
      bigout2.push_back(flip(bigout[i]));
      bigout2.back().rxn.net_reversible = -1;
      bigout2.back().rxn.lb = -1000.0f;
      bigout2.back().rxn.ub = 0.0f;
    }
  }
  return bigout2;
}
