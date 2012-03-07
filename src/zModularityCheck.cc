#include "Components.h"
#include "Exchanges.h"
#include "Grow.h"
#include "genericLinprog.h"
#include "Printers.h"
#include "Paths2Model.h"
#include "Modularity.h"
#include "RunK.h"
#include "XML_loader.h"

void HandWave(PROBLEM &theModel){
  /*
  int co2ex_idx = theModel.fullrxns.Ids2Idx[9406];
  theModel.fullrxns.rxns[co2ex_idx].net_reversible = 1;
  int GP4GH_idx = theModel.fullrxns.Ids2Idx[2247];
  theModel.fullrxns.rxns.erase(theModel.fullrxns.rxns.begin()+GP4GH_idx);
  theModel.fullrxns.rxnMap();
  */
  theModel.fullrxns.change_Lb_and_Ub(9326, 0, 777);
  theModel.fullrxns.change_Lb_and_Ub(9456, 0, 777);
  theModel.fullrxns.change_Lb_and_Ub(9460, 0, 777);
  theModel.fullrxns.change_Lb_and_Ub(9511, 0, 777);
  theModel.fullrxns.change_Lb_and_Ub(9674, 0, 777);
}

/* Indended to be used with the NEW files (ECOONLY) */
int main(int argc, char *argv[]){

  /*Load XML files*/
  printf("Starting program...\n");
  PROBLEM ProblemSpace;
  InputSetup(argc,argv,ProblemSpace);
  printf("XML files loaded\n");
  
  /*Load All Reactions in Model*/
  vector<REACTION> inModel;
  AllHardIncludes(ProblemSpace.fullrxns.rxns,inModel);
  printf("HardIncludes included\n");

  /*Load Exchanges (food) - shouldn't be needed, all rxns (included exchanges) were exported */

  LoadAllExchanges(ProblemSpace);
  for(int i=0;i<ProblemSpace.exchanges.rxns.size();i++){
    inModel.push_back(ProblemSpace.exchanges.rxns[i]);
    //printf("%d\n",ProblemSpace.exchanges.rxns[i].id);
  }
  printf("Exchanges Loaded\n");

  /*Load Biomass*/
  inModel.push_back(MakeBiomassFromGrowth(ProblemSpace.growth[0]));
  RXNSPACE rxnModel(inModel);
  //printf("vector<REACTION> --> RXNSPACE\n"); 

  /*Print everything rxnModel
  for(int i=0;i<rxnModel.rxns.size();i++){
    printf("%d\n",rxnModel.rxns[i].id);}
  printf("rxnModel Loaded\n");

  /*Alter the fluxes for exchanges*/
  printf("Fluxes Adjusted\n");
  PROBLEM theModel(rxnModel,ProblemSpace);
  FeedTheBeast(theModel.fullrxns,theModel.growth[0]);
  //HandWave(theModel);
  SynIncludes(theModel);
  theModel.synrxns.rxnMap();
  printf("Loaded Matrix ProblemSpace\n");

  for(int i=0;i<ProblemSpace.exchanges.rxns.size();i++){
    ProblemSpace.fullrxns.addReaction(ProblemSpace.exchanges.rxns[i]);
  }
  adjustLikelihoods(theModel.synrxns.rxns, 1.0f, -1.0f, 1.1f, -3.0f, false);   
  adjustLikelihoods(theModel.fullrxns.rxns, 1.0f, -1.0f, 1.1f, -3.0f, false);   

  /*
  printf("Everything: %d\n",(int)ProblemSpace.synrxns.rxns.size());
  printf("theModel: %d\n",(int)theModel.synrxns.rxns.size());
  //Print everything rxnModel

  /*Full FBA check*/
  //vector<double> fbaResult = FBA_SOLVE(ProblemSpace.fullrxns.rxns,ProblemSpace.metabolites);
  //printf("FBA check:\n"); for(int i=0;i<fbaResult.size();i++){printf("%f\n",fbaResult[i]);}
  /*Component by component FBA check*/
  //ComponentTest(ProblemSpace);

  /*Full inModel FBA check*/
  vector<double> fbaResult2 = FBA_SOLVE(theModel.fullrxns.rxns,theModel.metabolites);
  //printf("FBA check:\n"); for(int i=0;i<fbaResult2.size();i++){printf("%f\n",fbaResult2[i]);}
  /*Component by component FBA check*/
  ComponentTest(theModel);

  /*Loading Annotations*/
  //map<string,vector<int> > gene_map = GeneAnnotations();

  /*Reduce to Independent Pathways*/ 


  /*Report Pathways*/


}
