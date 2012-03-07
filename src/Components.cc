#include "Exchanges.h"
#include "RunK.h"
#include "pathUtils.h"
#include "XML_loader.h"
#include "Printers.h"
#include "Paths2Model.h"
#include "genericLinprog.h"
#include "Grow.h"
#include "kShortest.h"
#include "Components.h"
#include "visual01.h"
#include "MyConstants.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <queue>
#include <vector>

/* Tests Growth via FBA component by component for a growth condition */
int ComponentTest(const PROBLEM &theModel){
  int Kq = 1;
  REACTION biomass = theModel.fullrxns.rxnFromId(_db.BIOMASS);
  PROBLEM baseModel = theModel;
  //Remove Full Biomass
  int bIdx = baseModel.fullrxns.idxFromId(_db.BIOMASS);
  baseModel.fullrxns.rxns.erase(baseModel.fullrxns.rxns.begin()+bIdx);
  baseModel.fullrxns.rxnMap();
  int sIdx = baseModel.synrxns.idxFromId(_db.BIOMASS);
  baseModel.synrxns.rxns.erase(baseModel.synrxns.rxns.begin()+sIdx);
  baseModel.synrxns.rxnMap();
  //Put in Exchanges for the Biomass Components so that things can flux out
  /*
  for(int i=0;i<biomass.stoich.size();i++){
    METABOLITE met = baseModel.metabolites.metFromId(biomass.stoich[i].met_id);
    REACTION new_exchange = MagicExchange(met,1000,-1);
    baseModel.fullrxns.addReaction(new_exchange);
    baseModel.synrxns.addReaction(new_exchange);
  }
  */

  vector<int> inputIds;
  vector<int> dirs;
  for(int j=0;j<baseModel.growth[0].media.size();j++){
    inputIds.push_back(baseModel.growth[0].media[j].id);
    dirs.push_back(0);
  }
  METSPACE inputs(baseModel.metabolites, inputIds);  
  checkExchangesAndTransporters(baseModel, theModel, inputIds, dirs);
  vector< vector< PATH > > allPaths;
  adjustLikelihoods(baseModel.synrxns.rxns, 1.0f, -1.0f, 1.1f, -3.0f, false);
  /*
  printf("baseModel: %d\n",(int)baseModel.synrxns.rxns.size());
  /*Print everything rxnModel
  for(int i=0;i<baseModel.synrxns.rxns.size();i++){
    printf("%d %d  %s\n",i,baseModel.synrxns.rxns[i].id, baseModel.synrxns.rxns[i].name);}
  printf("synrxn Loaded\n");
  */
  for(int i=0;i<biomass.stoich.size();i++){
    PROBLEM tempModel = baseModel;

    /*Load one component Model*/
    STOICH tmpS = biomass.stoich[i]; vector<STOICH> tmpStoich; tmpStoich.push_back(tmpS);
    int outputId = tmpS.met_id;
    REACTION tmpRxn = MakeObjRxn(tmpStoich);
    tempModel.fullrxns.addReaction(tmpRxn);
    tempModel.synrxns.addReaction(tmpRxn);
    //Load_Stoic_Part(tempModel.fullrxns.rxns,tempModel.metabolites.mets);
    //Load_Stoic_Part(tempModel.synrxns.rxns,tempModel.metabolites.mets);

    /*Run FBA*/
     char *matlab_str = (char *) malloc(sizeof(char) * 64);
    sprintf(matlab_str,"matlab_out_%d.mat",i);
    MATLAB_out(matlab_str,tempModel.fullrxns.rxns);

    printf("Running FBA on %s\n",biomass.stoich[i].met_name);
    vector<double> fbaResult = FBA_SOLVE(tempModel.fullrxns.rxns,tempModel.metabolites);
    //Print
    /*
    printf("FBA result vector for %s\n",biomass.stoich[i].met_name);
   for(int j=0;j<fbaResult.size();j++){
      printf("%s %f %d\n",tempModel.fullrxns.rxns[j].name,fbaResult[j],tempModel.fullrxns.rxns[j].net_reversible);
    }
    */
    PROBLEM fluxonly = reportFlux(tempModel,fbaResult);
    /*
    for(int j=0;j<fluxonly.fullrxns.rxns.size();j++){
      printf("%s %d\n",fluxonly.fullrxns.rxns[j].name,fluxonly.fullrxns.rxns[j].net_reversible);
    }
    */
    if(fluxonly.fullrxns.rxns.size()>0){
      VisualizeProblem(fluxonly,100+i);}
    else{
      //Try free cofactors
      FeedEnergy(tempModel,1.0f);
      for(int j=0;j<tempModel.fullrxns.rxns.size();j++){
	printf("== %d %d %4.0f %4.0f %s\n",j,tempModel.fullrxns.rxns[j].id,
	       tempModel.fullrxns.rxns[j].lb, tempModel.fullrxns.rxns[j].ub, tempModel.fullrxns.rxns[j].name);
      }
      vector<double> fbaResult3 = FBA_SOLVE(tempModel.fullrxns.rxns,tempModel.metabolites);
      PROBLEM fluxonly3 = reportFlux(tempModel,fbaResult3);
      if(fluxonly3.fullrxns.rxns.size()>0){
	VisualizeProblem(fluxonly3,i+300);}
    }
    
    vector<PATH> kpaths;
    adjustLikelihoods(tempModel.synrxns.rxns, 1.0f, -1.0f, 1.1f, -3.0f, false);

    printf("Running kShortest: %d %d %d %s %d\n",(int)tempModel.synrxns.rxns.size(),
	   (int)tempModel.metabolites.mets.size(),(int)inputs.mets.size(),
	   tempModel.metabolites.metFromId(outputId).name, Kq);
    /*
    printf("tempModel: %d\n",(int)tempModel.synrxns.rxns.size());
    /*Print everything rxnModel
    for(int i=0;i<tempModel.synrxns.rxns.size();i++){
      printf("%d %d  %s\n",i,tempModel.synrxns.rxns[i].id, tempModel.synrxns.rxns[i].name);}
    printf("synrxn Loaded\n");
    */   
    kShortest(kpaths,tempModel.synrxns,tempModel.metabolites, inputs,
    	      tempModel.metabolites.metFromId(outputId), Kq);
    printPathResults(kpaths);
    allPaths.push_back(kpaths);
  }

  VisualizeAllSolutions(theModel,allPaths, true);
  return 0; /*all is well*/
}

void Load_Stoich_Part_With_Cofactors(PROBLEM &Model){
  for(int i=0;i<Model.fullrxns.rxns.size();i++){
    Model.fullrxns.rxns[i].stoich_part.clear();
    for(int j=0;j<Model.fullrxns.rxns[i].stoich.size();j++){
      int metID = Model.fullrxns.rxns[i].stoich[j].met_id;
      if(Model.metabolites.metFromId(metID).secondary_lone!=1){
	Model.fullrxns.rxns[i].stoich_part.push_back(Model.fullrxns.rxns[i].stoich[j]);
      }
    }
  }
  return;
}

PROBLEM reportFlux(PROBLEM& Model, vector<double> &fbaSolution){
  PROBLEM temp = reportFlux(Model,fbaSolution,0.000001f);
  return temp;
}

PROBLEM reportFlux(PROBLEM& Model, vector<double> &fbaSolution, double cutoff){
  PROBLEM reporter;

  Load_Stoich_Part_With_Cofactors(Model);
  for(int i=0;i<fbaSolution.size();i++){
    if(abs(fbaSolution[i])>cutoff){
      REACTION tmpRxn = Model.fullrxns.rxns[i];
      //tmpRxn.stoich_part = tmpRxn.stoich; //includes everything
      if(fbaSolution[i]>0.0f){ tmpRxn.net_reversible =  1; tmpRxn.lb = 0.0f; }
      if(fbaSolution[i]<0.0f){ tmpRxn.net_reversible = -1; tmpRxn.ub = 0.0f; }
      if(tmpRxn.stoich_part.size()>1){ //This line isn't registering with ac either...
	reporter.fullrxns.addReaction(tmpRxn);
	//printf("%s %d ( %d )\n",tmpRxn.name,tmpRxn.id,(int)tmpRxn.stoich_part.size());
	//for(int j=0;j<tmpRxn.stoich_part.size();j++){
	//  printf("   %s %d\n",tmpRxn.stoich_part[j].met_name,tmpRxn.stoich_part[j].met_id);
	//}
      }
    }
  }
  if(reporter.fullrxns.rxns.size()>0){
    vector<int> metList;
    for(int i=0;i<reporter.fullrxns.rxns.size();i++){
      for(int j=0;j<reporter.fullrxns.rxns[i].stoich_part.size();j++){
	metList.push_back(reporter.fullrxns.rxns[i].stoich_part[j].met_id);
	//printf("* %s  %d\n",rxnspace.rxns[i].stoich_part[j].met_name,
	//       rxnspace.rxns[i].stoich_part[j].met_id);
      }
    }
    custom_unique(metList);
    //printf("custom_unique\n");
    for(int i=0;i<metList.size();i++){
      //printf("* *  %d\n",metList[i]);
      reporter.metabolites.addMetabolite(Model.metabolites.metFromId(metList[i]));
    }
    reporter.metabolites.metMap();
  }
  return reporter;
}
