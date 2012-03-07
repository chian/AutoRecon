#include "Exchanges.h"
#include "MyConstants.h"
#include "Printers.h"
#include "RunK.h"

#include<cassert>
#include<cstdio>
#include<cstdlib>
#include<vector>

/* Functions */

// Loads all exchange reactions into the RXNSPACE and into REACTION list for matrix inclusion
void LoadAllExchanges(PROBLEM &ProblemSpace){
  vector<GROWTH> &growth = ProblemSpace.growth;
  vector<REACTION> &fullrxns = ProblemSpace.fullrxns.rxns;
  vector<REACTION> &exchanges = ProblemSpace.exchanges.rxns; //holds new reactions
  vector<int> foodlist;
  for(int i=0;i<growth.size();i++){
    for(int j=0;j<growth[i].media.size();j++){
      foodlist.push_back(growth[i].media[j].id);
    }
    for(int j=0;j<growth[i].byproduct.size();j++){
      foodlist.push_back(growth[i].byproduct[j].id);
    }
  }
  custom_unique(foodlist);
  vector<int> rev(foodlist.size(),1); 
  GetExchangeReactions(foodlist,rev,ProblemSpace,ProblemSpace.exchanges);
  return;
} 

// Zeros out the exchanges (reset)
void ResetFood(vector<REACTION> &reactions){
  for(int i=0;i<reactions.size();i++){
    if(reactions[i].isExchange==1){
      if(reactions[i].stoich[i].rxn_coeff<0){
	reactions[i].net_reversible = 1;
      }
      else{
	reactions[i].net_reversible = -1;
      }
    }
  }
}

// Adjusts Flux limits on Exchange Reactions for each Media condition
void FeedTheBeast(RXNSPACE &inModel, GROWTH &growth){
  ResetFood(inModel.rxns);
  for(int i=0;i<growth.media.size();i++){
    int j = FindExchange4Metabolite(inModel.rxns,growth.media[i].id);
    assert(j!=-1);

    printf("%d %d\n",growth.media[i].id,j);
    inModel.change_Lb_and_Ub(j, -growth.media[i].rate, 1000.0f);
  }
  return;
}

void FeedEnergy(PROBLEM &Model,double flux_bound){
  for(int i=0;i<Model.metabolites.mets.size();i++){
    if(Model.metabolites.mets[i].isSecondary()){
      int j = FindExchange4Metabolite(Model.fullrxns.rxns,Model.metabolites.mets[i].id);
      int id = j;
      if(j!=-1){
	REACTION new_exchange = MagicExchange(Model.metabolites.mets[i],1000.0f, 0);
	Model.fullrxns.addReaction(new_exchange);
	Model.synrxns.addReaction(new_exchange);
	id = new_exchange.id;
      }
      printf("* %d %d %s*\n",Model.metabolites.mets[i].id,j,
	     Model.metabolites.mets[i].name);
      Model.fullrxns.change_Lb_and_Ub(id, -flux_bound, flux_bound);
    }
  }
  return;
}

REACTION MagicExchange(METABOLITE met,double flux_bound,int dir){
  REACTION rxn_add;
  char *temp = (char *) malloc(sizeof(char) * 64);
  STOICH stoich_add;

  rxn_add.init_likelihood = -3;
  rxn_add.isExchange = 1;
  rxn_add.net_reversible = dir;
  rxn_add.id = _db.MISSINGEXCHANGEFACTOR + met.id;

  /* Size of string could be a problem but for now I left this since we don't want 30 extra spaces in everything... */
  sprintf(temp,"%s%s","MagicEx_",met.name);
  strcpy(rxn_add.name,temp);
  stoich_add.met_id = met.id;
  stoich_add.rxn_coeff = -1;
  sprintf(stoich_add.met_name,"%s",met.name);
  rxn_add.stoich.push_back(stoich_add);
  rxn_add.stoich_part.push_back(stoich_add);

  if(dir <= 0) {  rxn_add.lb = -flux_bound; } else { rxn_add.lb = 0.0f; }
  if(dir >= 0) {  rxn_add.ub = flux_bound;  } else { rxn_add.ub = 0.0f; }
  //Problem here is that there are sometimes transporters that involve multiple metabolites.
  //So I want to force the flux to use those if it has to rather than rely on the MagicTrans reactions
  free(temp);
  return rxn_add;
}

/* Returns the ID for the reaction that transports metabolite "met_id" in the vector<REACTION> 
   reaction, if there is one in that vector, and -1 if there is not 

   rxn_direction is -1 if the metabolite should be imported and 1 if the metabolite should be exported, hence
   if rxn_direction = -1 we are looking for rxn_coeff*net_reversibility <= 0 for the external metabolite (because it is imported) , and
   if rxn_direction = 1 we are looking for rxn_coeff*net_reversibility >= 0 

   Change 3-25-11: Return the most likely transporter */
int FindTransport4Metabolite(const vector<REACTION> &reaction, const METSPACE &metspace, int met_id, 
			     int rxn_direction){
  unsigned int i,j,k;
  int pair_id;
  int whichExternal;
  int result(-1);
  double currentMin(-1.0f);
  pair_id = inOutPair(met_id, metspace);

  if(isExternalMet(metspace.metFromId(met_id).name, _db.E_tag)) { whichExternal = met_id;  } 
  else {  whichExternal = pair_id;  }

  for(int i=0;i<reaction.size();i++){
    if(reaction[i].transporter==1){      
      for(int j=0;j<reaction[i].stoich.size();j++){
	if(reaction[i].stoich[j].met_id == met_id){ 
	  /* Check directionality in cases where the external met ID is provided */
	  if(whichExternal == met_id) {
	    if(!(reaction[i].net_reversible == 0) && /* Accounts for rounding error - allow reversible transporters when they are present */
	       !(rxn_direction == -1 && reaction[i].stoich[j].rxn_coeff*(double)reaction[i].net_reversible < 0) &&
	       !(rxn_direction == 1 && reaction[i].stoich[j].rxn_coeff*(double)reaction[i].net_reversible > 0) ) {
	      break;
	    }
	  }
	  /* Look for transporters for the correct metabolite */
	  for(k=0;k<reaction[i].stoich.size();k++){
	    if(reaction[i].stoich[k].met_id == pair_id) {
	      /* Check that the likelihood is actually higher than the current maximum */
	      if(currentMin < reaction[i].current_likelihood) {
		/* Check directionality in cases where the pairID is the external metabolite */
		if(whichExternal == pair_id) {
		  if(!(reaction[i].net_reversible == 0) && 
		     !(rxn_direction == -1 && reaction[i].stoich[k].rxn_coeff*(double)reaction[i].net_reversible < 0) &&
		     !(rxn_direction == 1 && reaction[i].stoich[k].rxn_coeff*(double)reaction[i].net_reversible > 0) ) {
		    break;
		  }
		}
		/* All checked out */
		result = reaction[i].id;
		currentMin = reaction[i].current_likelihood;
	      } else { /* Likelihood is too low */
		break;
	      }
	    }
	  }
	} /* if(reaction[i].stoich[j].met_id == met_id */
      }
    }
  }
  return result;
}

REACTION MagicTransport(const vector<REACTION> &reaction, const METSPACE &metspace, int met_id, 
			char* name, int R){
  return MagicTransport(reaction,metspace,met_id,&name[0],R,1000.0f);
}

REACTION MagicTransport(const vector<REACTION> &reaction, const METSPACE &metspace, int met_id, 
			char* name, int R, double bound){
  REACTION rxn_add;
  char *temp = (char *) malloc(sizeof(char) * 64);
  STOICH stoich_add;

  rxn_add.init_likelihood = -3;
  rxn_add.transporter = 1;
  rxn_add.net_reversible = R;
  rxn_add.id = _db.MISSINGTRANSPORTFACTOR + met_id;

  /* Size of string could be a problem but for now I left this since we don't want 30 extra spaces in everything... */
  sprintf(temp,"%s%s","MagicTran_",name);
  strcpy(rxn_add.name,temp);
  stoich_add.met_id = met_id;
  stoich_add.rxn_coeff = -1;
  sprintf(stoich_add.met_name,"%s",name);
  rxn_add.stoich.push_back(stoich_add);
  rxn_add.stoich_part.push_back(stoich_add);
  stoich_add.met_id = inOutPair(met_id,metspace);
  if(stoich_add.met_id==-1){
    printf("MagicTransport: ERROR - NO COMPLEMENTARY METABOLITE: %d %s\n",met_id,name);}
  stoich_add.rxn_coeff = 1;
  sprintf(stoich_add.met_name, "%s", metspace.metFromId(stoich_add.met_id).name);
  rxn_add.stoich.push_back(stoich_add);
  rxn_add.stoich_part.push_back(stoich_add);
  if(R <= 0) { rxn_add.lb = -bound; } else { rxn_add.lb = 0.0f; }
  if(R >= 0) { rxn_add.ub =  bound; } else { rxn_add.ub = 0.0f; }

  //Problem here is that there are sometimes transporters that involve multiple metabolites.
  //So I want to force the flux to use those if it has to rather than rely on the MagicTrans reactions
  free(temp);
  return rxn_add;
}


/* Return a vector of transporters from synRxns for all byproducts or media in the defined GROWTH */
void GetTransportReactionsFromGrowth(const GROWTH &growth, vector<int> rxnDirections, 
				     const PROBLEM &ProblemSpace, vector<REACTION> &Transporters) {
  unsigned int i;
  Transporters.clear();
  vector<int> metIds;
  for(i=0; i<growth.media.size();i++) {
    metIds.push_back(growth.media[i].id);
  }
  for(i=0; i<growth.byproduct.size(); i++) {
    metIds.push_back(growth.byproduct[i].id);
  }
  GetTransportReactions(metIds, rxnDirections, ProblemSpace, Transporters);
  return;
}

/* General function to check for transport reactions in reaction vector reaction (or any reaction vector really)
   and returns the transporters if they are present in the vector, and magic transporters if they are not. 
   SKIPS OVER H for now - can define the bahavior later but for now the indexing of metIds will not be consistent
   with the indexing of the vector<REACTION> 

   rxnDirections: should be -1 for things that need to IMPORT something into the cell and +1 for things that need to EXPORT something

   It is NOT recommended to call this with synreaction because you will most likely end up with transporters with the wrong reversibility. */
void GetTransportReactions(vector<int> metIds, vector<int> rxnDirections, const PROBLEM &ProblemSpace,
		      vector<REACTION> &Transporters){
  Transporters.clear();
  vector<int> ExchangeRxns;
  REACTION TMPRXN;
  double tmpLikely = 0.0f;

  const vector<REACTION> &reaction = ProblemSpace.fullrxns.rxns;
  const vector<METABOLITE> &metabolite = ProblemSpace.metabolites.mets;
  const METSPACE &metspace = ProblemSpace.metabolites;

  for(int i=0;i<metIds.size();i++) {
    if(strcmp(metspace.metFromId(metIds[i]).name, _db.H_name) != 0) {
      /* FindTransport4Metabolite already looks for the most likely and also will account for direction. */
      int temp = FindTransport4Metabolite(reaction,metspace,metIds[i],rxnDirections[i]);
      if(temp == -1){
	TMPRXN = MagicTransport(reaction,metspace,metIds[i],metspace.metFromId(metIds[i]).name,0);
      } else {
	TMPRXN = ProblemSpace.fullrxns.rxnFromId(temp);
      }
    }
    Transporters.push_back(TMPRXN);
  }
  return;
}

/* For all of the metabolites with ID listed in metIdList, 
identifies if there is an exchange for that metabolite in the list of reactions
in ProblemSpace. If there is one, it puts that REACTION in exchanges. If there ISN'T one, 
it creates one (with the black magic flag to distinguish it) and then places it in
exchanges. 
*/
void GetExchangeReactions(vector<int> metIdList, vector<int> dirs, const PROBLEM &ProblemSpace, RXNSPACE &exchanges) {

  const METSPACE &metspace = ProblemSpace.metabolites;
  const RXNSPACE &rxnspace = ProblemSpace.fullrxns;

  exchanges.clear();

  char *name = (char*) malloc(sizeof(char)*64);

  for(int i=0;i<metIdList.size();i++) {
    /* Returns an existing reaction ID if possible */
    int temp = FindExchange4Metabolite(rxnspace.rxns, metIdList[i]);
    REACTION TMPRXN;
    METABOLITE tmpMet;

    if(temp==-1) {
      if(metspace.idIn(metIdList[i])) {
	sprintf(name, "%s", ProblemSpace.metabolites.metFromId(metIdList[i]).name);
	TMPRXN = GrowthExit(rxnspace.rxns, metIdList[i], dirs[i], 1000, name); 
      } else {
	/* Add the exchange anyway, but warn the user about possible perils... */
	printf("WARNING: In AddExchangeReactions a metID was passed that has no corresponding elements in the metabolite struct. Adding an entry for it and recalculating... \n");
	printf("You should replace this with the METABOLITE entry in the full metabolite vector later. If you are using the full metabolite vector this indicates a major problem. \n");
	printf("Offending metID: %d\n", metIdList[i]);
	tmpMet.id = metIdList[i];
	sprintf(tmpMet.name, "UNKNOWN_ID %d", tmpMet.id);
	sprintf(name, "%s", tmpMet.name);
	TMPRXN = GrowthExit(rxnspace.rxns, metIdList[i], dirs[i], 1000, name);
      }
    } else {
      TMPRXN = rxnspace.rxnFromId(temp);
      TMPRXN.net_reversible = dirs[i];
    }
    TMPRXN.init_likelihood = -2.0f;
    exchanges.addReaction(TMPRXN);
  }
  free(name);
  return;
}

/* Add Exchange Rxns for BYPRODUCTS from a SINGLE GROWTH (same syntax as oneMedia for consistency... */
void GetExchangeReactions_byproducts(const GROWTH &growth, const PROBLEM &ProblemSpace, RXNSPACE &eList){
  vector<int> metIdList;
  vector<int> dirs;
  eList.rxns.clear();
  for(int j=0;j<growth.byproduct.size();j++) {
    metIdList.push_back(growth.byproduct[j].id);
    dirs.push_back(1);
  }
  GetExchangeReactions(metIdList, dirs, ProblemSpace, eList);
  return;
}

/* Add Exchange Rxns for MEDIA from a SINGLE GROWTH */
void GetExchangeReactions_oneMedia(const GROWTH &growth, const PROBLEM &ProblemSpace, RXNSPACE &eList){
  vector<int> metIdList;
  vector<int> dirs;
  eList.rxns.clear();
  for(int j=0;j<growth.media.size();j++) {
    metIdList.push_back(growth.media[j].id);
    dirs.push_back(-1); /* Consumed only */
  }
  GetExchangeReactions(metIdList, dirs, ProblemSpace, eList);
  return;
}


/* Returns the ID for the exchange reaction for metabolite "met_id" in the vector<REACTION> 
   reaction, if there is one in that vector, and -1 if there is not 
   Does not depend on the "isExchange" property */
int FindExchange4Metabolite(const vector<REACTION> &reaction, int met_id){
  for(int i=0;i<reaction.size();i++){
     if(reaction[i].stoich[0].met_id==met_id && reaction[i].stoich.size()==1 && reaction[i].id!=_db.BIOMASS){ 
       return reaction[i].id;}
  }
  return -1; /* Exchange not found in database */
}
