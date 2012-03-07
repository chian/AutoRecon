#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <vector>

#include "DataStructures.h"
#include "Printers.h"

using std::map;
using std::set;
using std::vector;

/**************** Generic printers **************/

void printIntMap(map<int, int> intMap) {
  map<int, int>::const_iterator i;
  if(intMap.empty()) {printf("EMPTY\n"); return; }
  for(i=intMap.begin(); i != intMap.end(); ++i) {
    printf("%d\t%d\n", i->first, i->second);
  }
  return;
}


void printIntVector(vector<int> intVector) {
  if(intVector.empty()) {printf("EMPTY\n"); return;}
  for(int i=0;i<intVector.size();i++)  { printf("%d ", intVector[i]);  }
  printf("\n");
  return;
}

void printDoubleVector(vector<double> doubleVector) {
  if(doubleVector.empty()) {printf("EMPTY\n"); return;}
  for(int i=0;i<doubleVector.size();i++)    {
    printf("%4.3f ", doubleVector[i]);
  }
  printf("\n");
  return;
}

void printIntSet(const set<int> &mySet) {
  if(mySet.empty()) { printf("EMPTY\n"); return; }
  for(set<int>::iterator it=mySet.begin(); it!=mySet.end(); it++) {
    printf("%d ", *it);
  }
  printf("\n");
  return;
}

/***************** Metablite / reaction printers ****************/

void printRxnFormula(const REACTION &rxn, char* rxnString, bool printStoichPart) {
  vector<STOICH> st;
  if(printStoichPart) {
    st = rxn.stoich_part;
  } else {
    st = rxn.stoich;
  }
  sort(st.begin(), st.end());
  bool arrowDone = false; /* Needed to make sure we always get an arrow including on exchanges */
  for(int i=0; i < st.size(); i++) {
    double rxnCoeff = st[i].rxn_coeff;
    if(rxnCoeff < 0.0f) {
      rxnCoeff *= -1;
    }
    char oneReactantString[96];
    sprintf(oneReactantString, "%1.4f %s", rxnCoeff, st[i].met_name);
    strcat(rxnString, oneReactantString);

    if(i != st.size() - 1) {
      if(st[i+1].rxn_coeff * st[i].rxn_coeff < 0.0f) {
        if(rxn.net_reversible == -1) {
          strcat(rxnString, " <-- ");
        } else if(rxn.net_reversible == 0) {
          strcat(rxnString, " <=> ");
        } else {
          strcat(rxnString, " --> ");
        }
        arrowDone = true;
      } else {
        strcat(rxnString, " + ");
      }
    }
  }
  if(!arrowDone) {
    if(rxn.net_reversible == -1) {
      strcat(rxnString, " <-- ");
    } else if(rxn.net_reversible == 0) {
      strcat(rxnString, " <=> ");
    } else {
      strcat(rxnString, " --> ");
    }
  }
  return;
}


void printMETABOLITEinputs(const METABOLITE &metabolite){
  printf("\tid: %05d\n",metabolite.id);
  printf("\tname: %s\n",metabolite.name);
  printf("\tcharge: %d\n",metabolite.charge);
  printf("\tinput: %d  ",metabolite.input);
  if(metabolite.input==0){ printf("(NO)\n");} else{printf("(YES)\n");}
  printf("\toutput: %d  ",metabolite.output);
  if(metabolite.output==0){ printf("(NO)\n");} else{printf("(YES)\n");}
  printf("\tbiomass: %d  ",metabolite.biomass);
  if(metabolite.biomass==0){ printf("(NO)\n");} else{printf("(YES)\n");}
  printf("\tsecondary_lone: %d  ",metabolite.secondary_lone);
  if(metabolite.secondary_lone==0){ printf("(NO)\n");} else{printf("(YES)\n");}
  if(metabolite.secondary_pair.size()>0){
    for(int j=0;j<metabolite.secondary_pair.size();j++){
      printf("\tsecondary_pair: %d  ",metabolite.secondary_pair[j]);
    }
  }
  else{printf("(YES - PARTNER ID# GIVEN)\n");}
  printf("\n");
}

void printGROWTHinputs(const GROWTH &growth){
  int i;
  printf("GROWTH MEDIA:\n");
  for(i=0;i<growth.media.size();i++){
    printf("\t%05d\t%1.3f\n",growth.media[i].id,growth.media[i].rate);
  }
  printf("\nBYPRODUCTS:\n");
  for(i=0;i<growth.byproduct.size();i++){
    printf("\t%05d\t%1.3f\n",growth.byproduct[i].id,
           growth.byproduct[i].rate);
  }
  printf("\nMUTATIONS:\n");
  for(i=0;i<growth.mutation.size();i++){
    printf("\t%s\t%1.3f\n",growth.mutation[i].genename.c_str(),
           growth.mutation[i].act_coef);
  }
  printf("\tnum_of_metabolites: %d\n",(int)growth.biomass.size());
  if(growth.biomass.size() > 0){
    printf("\n BIOMASS EQUATION:\n\n");
    printSTOICHIOMETRY_by_id(growth.biomass,growth.biomass.size(),1);}
  printf("\nGROWTH RATE: %f\n\n",growth.growth_rate);
}

void printSynRxns(const RXNSPACE &synrxns, const RXNSPACE &fullrxns) {
  for(int i=0; i<synrxns.rxns.size(); i++) {
    printf("Reaction %s (rev=%d) synrxns: ", synrxns.rxns[i].name, synrxns.rxns[i].net_reversible);
    for(int j=0; j<synrxns.rxns[i].syn.size(); j++) {
      printf("%s (rev=%d); ", fullrxns.rxnFromId(synrxns.rxns[i].syn[j]).name, fullrxns.rxnFromId(synrxns.rxns[i].syn[j]).net_reversible);
    }
    printf("\n");
  }
}

void printREACTIONintermediates(const REACTION &reaction, int print_type){
  printf("\tid: %05d\n",reaction.id);
  printf("\tname: %s\n",reaction.name);
  printf("\treversible: %d  ",reaction.net_reversible);
  if(reaction.net_reversible==0){ printf("(YES)\n");}
  if(reaction.net_reversible==-1){ printf("(NO - BACKWARDS)\n");}
  if(reaction.net_reversible==1){ printf("(NO - FORWARDS)\n");}
  printf("\tSynthesis reaction:");
  if(reaction.synthesis == -1) { printf("(NO)\n"); }
  else { printf("YES - metabolite %d\n", reaction.synthesis); }

  printf("\t [lb: %f; ub: %f]\n", reaction.lb, reaction.ub);
  printf("\ttransporter: %d  ",reaction.transporter);
  if(reaction.transporter==0){ printf("(NO)\n");} else{printf("(YES)\n");}
  printf("\tlikelihood_init: %f\n",reaction.init_likelihood);
  printf("\tlikelihood_curr: %f\n",reaction.current_likelihood);
  printf("\tnum_of_metabolites: %d\n",(int)reaction.stoich.size());
  if(reaction.stoich.size() > 0){
    if(print_type==1){
      printf("\n STOICHIOMETRY:\n\n");
      printSTOICHIOMETRY_by_id(reaction.stoich,reaction.stoich.size(),
                               reaction.init_reversible);}
  }
  printf("\n");
  printf("\tnum_of_non_secondary_metabolites: %d\n", (int)reaction.stoich_part.size());
  if(reaction.stoich_part.size() > 0){
    if(print_type==1){
      printf("\n STOICHIOMETRY:\n\n");
      printSTOICHIOMETRY_by_id(reaction.stoich_part,
                               reaction.stoich_part.size(),
                               reaction.init_reversible);}
  }
  printf("\n");
  if(reaction.syn.size()>1){
    printf("\tSYNONYMS: ");
    for(int i=0;i<reaction.syn.size();i++){
      printf("%05d, ",reaction.syn[i]);
    }
    printf("\n\n--------------------\n");
  }
  return;
}


void printSTOICHIOMETRY_by_id(const vector<STOICH> &stoich, int num_met, int rev){
  unsigned int i;
  STOICH temps;
  vector<STOICH> tempv = stoich;

  /* Quick little sort to make sure negatives are before positives */
  for(i=1;i<tempv.size();){
    if(tempv[i-1].rxn_coeff <= tempv[i].rxn_coeff){ i++; continue; }
    if(tempv[i-1].rxn_coeff > tempv[i].rxn_coeff){
      temps = tempv[i];
      tempv[i] = tempv[i-1];
      tempv[i-1] = temps;
      if(i>1){i--;}
      continue;
    }
  }

  /* Print reaction */
  i=0;
  while(tempv[i].rxn_coeff<0){
    printf("%2.2f %s ",tempv[i].rxn_coeff,tempv[i].met_name);
    if(tempv[i+1].rxn_coeff<0){ printf("+ ");}
    i++;}
  if(rev==-1){ printf(" <--  ");}
  if(rev==0) { printf(" <--> ");}
  if(rev==1) { printf("  --> ");}
  while(tempv[i].rxn_coeff>0 && i<tempv.size()){
    printf("%2.2f %s ",tempv[i].rxn_coeff,tempv[i].met_name);
    if(i!=(tempv.size()-1)){ printf("+ ");}
    i++;}
  printf("\n");
}

void printREACTIONvector(const vector<REACTION> &reaction, int print_type) {
  for(unsigned int i=0; i<reaction.size(); i++) { printREACTIONinputs(reaction[i], print_type);}
}

void printREACTIONinputs(const REACTION &reaction, int print_type){
  printf("\tid: %05d\n",reaction.id);
  printf("\tname: %s\n",reaction.name);
  printf("\tinit_reversible: %d  ",reaction.init_reversible);
  if(reaction.init_reversible==0){ printf("(YES)\n");}
  if(reaction.init_reversible==-1){ printf("(NO - BACKWARDS)\n");}
  if(reaction.init_reversible==1){ printf("(NO - FORWARDS)\n");}

  printf("\tFree make: %d ", reaction.freeMakeFlag);
  if(reaction.freeMakeFlag == 1) { printf("(YES)\n"); } else { printf("(NO)\n"); }

  printf("\tSynthesis reaction:");
  if(reaction.synthesis == -1) { printf("(NO)\n"); }
  else { printf("YES - metabolite %d\n", reaction.synthesis); }

  printf("\tnet_reversible: %d \n", reaction.net_reversible);
  printf("\tLower bound: %f\n", reaction.lb);
  printf("\tUpper bound: %f\n", reaction.ub);
  printf("\ttransporter: %d  ",reaction.transporter);
  if(reaction.transporter==0){ printf("(NO)\n");} else{printf("(YES)\n");}
  printf("\tlikelihood: %f",reaction.init_likelihood);
  if(reaction.init_likelihood < 0.0f){
    printf("(SPECIAL FLAG - ");
    if(reaction.init_likelihood > -1.1f && reaction.init_likelihood < -0.9f){
      printf("DO NOT INCLUDE IN NETWORK)\n");}
    if(reaction.init_likelihood > -2.1f && reaction.init_likelihood < -1.9f){
      printf("USER DEFINED INCLUDE IN NETWORK)\n");}
    if(reaction.init_likelihood > -3.1f && reaction.init_likelihood < -2.9f){
      printf("BLACK MAGIC)\n");}
    if(reaction.init_likelihood > -4.1f && reaction.init_likelihood < -3.9f){
      printf("SPONTANEOUS REACTION)\n");}
    if(reaction.init_likelihood > -5.1f && reaction.init_likelihood < -4.9f){
      printf("NOT FOUND)\n");}
  }
  else{ printf("\n");}
  printf("\tnum_of_metabolites: %d\n",(int)reaction.stoich.size());
  if(reaction.stoich.size() > 0){
    if(print_type==1){
      printf("\n STOICHIOMETRY:\n\n");
      printSTOICHIOMETRY_by_id(reaction.stoich,reaction.stoich.size(),
                               reaction.init_reversible);}
  }
  printf("\n");
  printf("\tnum_of_non_secondary_metabolites: %d\n", (int)reaction.stoich_part.size());
  if(reaction.stoich_part.size() > 0){
    if(print_type==1){
      printf("\n STOICHIOMETRY:\n\n");
      printSTOICHIOMETRY_by_id(reaction.stoich_part,
                               reaction.stoich_part.size(),
                               reaction.init_reversible);}
  }
  printf("\n");
}


void printRxnsFromIntVector(const vector<int> &intVector, const RXNSPACE &rxnspace) {
  if(intVector.empty()) {printf("EMPTY\n"); return;}
  for(int i=0;i<intVector.size();i++){
    printf("%s(%4.3f) ", rxnspace.rxnFromId(abs(intVector[i])).name, rxnspace.rxnFromId(abs(intVector[i])).init_likelihood);
  }
  printf("\n");
}

void printRxnsFromIntSet(const set<int> &intSet, const RXNSPACE &rxnspace) {
  if(intSet.empty()) { printf("EMPTY\n"); return; }
  for(set<int>::iterator it=intSet.begin(); it!=intSet.end(); it++) {
    printf("%s(%4.3f) ", rxnspace.rxnFromId(*it).name, rxnspace.rxnFromId(*it).init_likelihood);
  }
  printf("\n");
  return;
}

void printMetsFromIntVector(const vector<int> &intVector, const PROBLEM &ProblemSpace) {
  if(intVector.empty()) {printf("EMPTY\n"); return;}
  for(int i=0;i<intVector.size();i++){
    printf("%s ", ProblemSpace.metabolites.metFromId(intVector[i]).name);
  }
  printf("\n");
}

/* The metabolites in RXNSPACE should have the same order as the doubles in doubleVec
   Prints out the reaction names and then the double values in a nice format (omits zero values) */
void printDoubleVector_rxns(const RXNSPACE &rxnspace, const vector<double> &doubleVec) {
  if(rxnspace.rxns.size() != doubleVec.size()) {
    printf("WARNING: Attempt to call printDoubleVector_rxns with a RXNSPACE and double vector of different sizes\n");
    return;
  }
  for(int i=0; i<doubleVec.size(); i++) {
    if(doubleVec[i] > 0.0001 || doubleVec[i] < -0.0001) {
      printf("%s\t%1.3f\n", rxnspace.rxns[i].name, doubleVec[i]);
    }
  }
  return;
}

/* The metabolites in METSPACE should have the same order as the doubles in doubleVec
   Prints out the metbaolite names and then the double values in a nice format */
void printDoubleVector_mets(const METSPACE &metspace, const vector<double> &doubleVec) {
  if(metspace.mets.size() != doubleVec.size()) {
    printf("WARNING: Attempt to call printDoubleVector_mets with a METSPACE and double vector of different sizes\n");
    return;
  }
  for(int i=0; i<doubleVec.size(); i++) {
    if(doubleVec[i] > 0.0001 || doubleVec[i] < -0.0001) {
      printf("%s\t%1.3f\n", metspace.mets[i].name, doubleVec[i]);
    }
  }
  return;
}


void printMetsFromStoich(vector<STOICH> a){
  int i;
  for(i=0;i<a.size();i++){
    printf("%s ",a[i].met_name);
  }
  printf("\n");
  return;
}


/************* ETC printers ***********************/

void printNetReactionVector(const vector<NETREACTION> &netReactions, const PROBLEM &problemSpace) {
  for(int i=0; i<netReactions.size(); i++) {
    printf("ETC %d: ", i);
    for(int j=0; j<netReactions[i].rxnDirIds.size(); j++) {
      REACTION tmpRxn = problemSpace.fullrxns.rxnFromId(abs(netReactions[i].rxnDirIds[j]));
      printf("%s(%4.3f)\t", tmpRxn.name, tmpRxn.init_likelihood);
    }
    printf("\n");
  }
  return;
}


/*************** Path printers ******************/

void PrintPathSummary_verbose(const PATHSUMMARY &psum) {
  printf("----------------------------------------\n");
  printf("RxnDirIds:");
  printIntVector(psum.rxnDirIds);
  printf("deadEndIds:");
  printIntVector(psum.deadEndIds);
  printf("AllMetsConsumed:");
  printIntVector(psum.metsConsumed);
  printf("Output ID: %d\n", psum.outputId);
  printf("Growth Index: %d\n", psum.growthIdx[0]);
  printf("K number: %d\n", psum.k_number);
}

void PrintPathSummary(const vector<vector<vector<PATHSUMMARY> > > &psum){
  /* Standard Dummies */
  unsigned int i,j,k;
  for(i=0;i<psum.size();i++){
    for(j=0;j<psum[i].size();j++){
      for(k=0;k<psum[i][j].size();k++){
        printf("-%d-%d-%d- ",i,j,k);
        printIntVector(psum[i][j][k].rxnDirIds);
      }
    }
  }
  return;
}

void printPathResults(const vector<PATH> &path) {
  for(unsigned int i=0;i<path.size();i++) {
    printf("Path number: %d corresponding to output %d...\n", i, path[i].outputId);
    printf("Inputs required to reach output: ");
    printIntVector(path[i].inputIds);
    printf("Number of reactions: ");
    printf("%d\n", (int)path[i].rxnIds.size());
    printf("Reactions required: ");
    printIntVector(path[i].rxnIds);
    printf("Reaction directionality: ");
    printIntVector(path[i].rxnDirection);
    printf("Total likelihood: %4.3f \n\n", path[i].totalLikelihood);
    printf("Dead ends: ");
    printIntVector(path[i].deadEndIds);
    printf("---------------------------------------------------------------------\n\n");
  }
}

/* Prints the metabolite names in addition to just the IDs */
void printPathResults(const vector<PATH> &path, PROBLEM &ProblemSpace) {
  for(unsigned int i=0;i<path.size();i++) {
    printf("Path number: %d corresponding to output %d...\n", i, path[i].outputId);
    printf("Inputs required to reach output: "); 
    printMetsFromIntVector(path[i].inputIds,ProblemSpace);
    printf("Number of reactions: ");
    printf("%d\n", (int)path[i].rxnIds.size());
    printf("Reactions required: ");
    printRxnsFromIntVector(path[i].rxnIds,ProblemSpace.synrxns);
    printf("Reaction directionality: ");
    printIntVector(path[i].rxnDirection);
    printf("Total likelihood: %4.3f \n\n", path[i].totalLikelihood);
    printf("Dead ends: ");
    printMetsFromIntVector(path[i].deadEndIds,ProblemSpace);
    printf("---------------------------------------------------------------------\n\n");
  }
  return;
}

/* ALso prints reaction likelihoods from rxnspace (why do we need this if the rxnspace is also found in PROBLEMSPACE?) */
void printPathResults(const vector<PATH> &path, PROBLEM &ProblemSpace, RXNSPACE &rxnspace) {
  for(unsigned int i=0;i<path.size();i++) {
    printf("Path number: %d corresponding to output %s...\n", i, ProblemSpace.metabolites.metFromId(path[i].outputId).name);
    printf("Inputs required to reach output: "); 
    printMetsFromIntVector(path[i].inputIds,ProblemSpace);
    printf("Number of reactions: ");
    printf("%d\n", (int)path[i].rxnIds.size());
    printf("Reactions required: ");
    printRxnsFromIntVector(path[i].rxnIds,rxnspace);
    printf("Reaction likelihoods: ");
    for(int j=0;j<path[i].rxnIds.size();j++){
      printf("%1.2f  ",rxnspace.rxnFromId(path[i].rxnIds[j]).current_likelihood);}
    printf("\n");
    printf("Reaction directionality: ");
    printIntVector(path[i].rxnDirection);
    printf("Total likelihood: %4.3f \n\n", path[i].totalLikelihood);
    printf("Dead ends: ");
    printMetsFromIntVector(path[i].deadEndIds,ProblemSpace);
    printf("---------------------------------------------------------------------\n\n");
  }
  return;
}

void PrintPathSummary2(const vector<vector<vector<PATHSUMMARY> > > &psum, RXNSPACE &rxnspace){
  for(int i=0;i<psum.size();i++){
    for(int j=0;j<psum[i].size();j++){
      for(int k=0;k<psum[i][j].size();k++){
	printf("--------\n");
        printf("-%d-%d-%d- ",i,j,k);
        printRxnsFromIntVector(psum[i][j][k].rxnDirIds,rxnspace);
	printf("Priorities: \n");
	printRxnsFromIntVector(psum[i][j][k].rxnPriority, rxnspace);
      }
    }
  }
  return;
}

void PrintPathSummary2(const vector<vector<vector<PATHSUMMARY> > > &psum, PROBLEM &ProblemSpace){
  PrintPathSummary2(psum,ProblemSpace.synrxns);
  return;
}

void PrintGapfillResult(const vector<GAPFILLRESULT> &res, const PROBLEM &problemSpace, const vector<int> &kToPrint) {
  for(int j=0; j<res.size(); j++) {
    if(res[j].deadEndSolutions.size() == 0) { continue; }
    printf("%s (k=%d) --> ", problemSpace.metabolites.metFromId(res[j].deadMetId).name, kToPrint[j]);
    for(int n=0; n<res[j].deadEndSolutions[kToPrint[j]].size(); n++) {
      REACTION tmp = problemSpace.fullrxns.rxnFromId(res[j].deadEndSolutions[kToPrint[j]][n]);
      printf("%s(%4.3f)\t", tmp.name, tmp.init_likelihood);
    }
    printf("\n");
  }
}

/******************** Output files ***************/

/* Output reactions in InRxns to a text file (including name and direction) */
void MATLAB_out(const char* fileName, const vector<REACTION> &InRxns){
  unsigned int i,j,k;
  STOICH temps;
  vector<STOICH> stoich;
  FILE* output;
  output = fopen(fileName, "w");

  for(i=0;i<InRxns.size();i++){

    stoich = InRxns[i].stoich;
    /* Sort stoichs by coefficient */
    for(j=1;j<stoich.size();){
      if(stoich[j-1].rxn_coeff <= stoich[j].rxn_coeff){ j++; continue; }
      if(stoich[j-1].rxn_coeff > stoich[j].rxn_coeff){
        temps = stoich[j];
        stoich[j] = stoich[j-1];
        stoich[j-1] = temps;
        if(j>1){j--;}
        continue;
      }
    }

    bool allPositive(true);
    for(j=0; j<stoich.size(); j++) {
      if(stoich[j].rxn_coeff < 0.0f) {
        allPositive = false;
        break;
      }
    }

    /* 1st column: Print out name */
    fprintf(output,"%s\t",InRxns[i].name);

    /* 2nd column: Print out reaction */
    for(j=0,k=0;j<stoich.size();j++){
      if(stoich[j].rxn_coeff < 0.0f){
        fprintf(output,"%1.8f %s",-stoich[j].rxn_coeff,stoich[j].met_name);
        if((j+1) < stoich.size() && stoich[j+1].rxn_coeff < 0.0f){ fprintf(output," + ");}
      }
      /* If either 1) the sign changes between the current and next STOICH, or
         2) there's only positive coefficients [and thus we need to put an arrow --> before the first one] */
      if( ((j+1) < stoich.size() && stoich[j].rxn_coeff < 0.0f && stoich[j+1].rxn_coeff > 0.0f ) ||
          ( j==0 && allPositive) ){
        k++;
        if(InRxns[i].net_reversible==0){fprintf(output," <--> ");}
        if(InRxns[i].net_reversible==1){fprintf(output," --> ");}
        if(InRxns[i].net_reversible==-1){fprintf(output," <-- ");}
      }
      if(stoich[j].rxn_coeff > 0.0f){
        fprintf(output,"%1.8f %s",stoich[j].rxn_coeff,stoich[j].met_name);
        if((j+1) < stoich.size()){ fprintf(output," + ");}
      }
    }
    if(k==0){
      if(InRxns[i].net_reversible==0){fprintf(output," <--> ");}
      if(InRxns[i].net_reversible==1){fprintf(output," --> ");}
      if(InRxns[i].net_reversible==-1){fprintf(output," <-- ");}
    }
    fprintf(output, "\t");

    /* 3rd\4th column: Print out the flux bound */
    fprintf(output, "%1.5f\t%1.5f\t", InRxns[i].lb, InRxns[i].ub);

    /* 6th/6th column: Print out original and final reversibility */
    fprintf(output, "%d\t%d\t", InRxns[i].init_reversible, InRxns[i].net_reversible);

    /* 7th/8th column: Print out initial and current likelihood */
    fprintf(output, "%1.5f\t%1.5f", InRxns[i].init_likelihood, InRxns[i].current_likelihood);

    /* End of reaction entry */
    fprintf(output,"\n");
  }

  fflush(output);
  fclose(output);
  return;
}

/* Report reactions involved in psum (AFTER unsynonymizing - if you try this before unsyn it will die a horrible death) */
void PATHS_rxns_out(const char* fileName, const vector<PATHSUMMARY> &psum, const PROBLEM &problem) {

  const RXNSPACE &rxnspace = problem.fullrxns;
  const METSPACE &metspace = problem.metabolites;

  /* Reactions */
  FILE* output = fopen(fileName, "w");

  fprintf(output, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	  "PSUM ID", "GROWTH CONDITION", "TARGET METABOLITE", "K (#th shortest)", "REACTION", "RXNDIRID", "LIKELIHOOD");

  for(int i=0; i<psum.size(); i++) {
    for(int j=0; j<psum[i].rxnDirIds.size(); j++) {
      fprintf(output, "%ld\t%d\t%s\t%d\t%s\t%d\t%4.3f\n",
	      psum[i].id, psum[i].growthIdx[0], metspace.metFromId(psum[i].outputId).name, psum[i].k_number, rxnspace.rxnFromId(abs(psum[i].rxnDirIds[j])).name, psum[i].rxnDirIds[j], 
	      rxnspace.rxnFromId(abs(psum[i].rxnDirIds[j])).init_likelihood);
    }
  }

  fflush(output);
  fclose(output);
}

/* Report list of metabolites in each PATH. The reaction names in psum
 must correspond to the reactions in rxnspace (i.e. for un-synnonymized RXNSPACE, the RXNSPACE should be fullrxns,
 and for the synonymized version from Dijkstras it should be synrxns, etc.

This format can be easily sorted and checked e.g. for what reactions are foudn in what paths */
void PATHS_mets_out(const char* fileName, const vector<PATHSUMMARY> &psum, const PROBLEM &problem) {
  const RXNSPACE &rxnspace = problem.fullrxns;
  const METSPACE &metspace = problem.metabolites;
  /* Reactions */
  FILE* output = fopen(fileName, "w");

  fprintf(output, "%s\t%s\t%s\t%s\t%s\n",
	  "PSUM ID", "GROWTH CONDITION", "TARGET METABOLITE", "K (#th shortest)", "METABOLITE", "REACTION containing metabolite");

  for(int i=0; i<psum.size(); i++) {
    for(int j=0; j<psum[i].rxnDirIds.size(); j++) {
      REACTION tmp = rxnspace.rxnFromId(abs(psum[i].rxnDirIds[j]));
      for(int k=0; k<tmp.stoich.size(); k++) {
	fprintf(output, "%ld\t%d\t%s\t%d\t%s\t%s\n",
		psum[i].id, psum[i].growthIdx[0], metspace.metFromId(psum[i].outputId).name, psum[i].k_number, metspace.metFromId(tmp.stoich[k].met_id).name, tmp.name);
      }
    }
  }

  fflush(output);
  fclose(output);
}


void ANNOTATIONS_out(const char* filename, const vector<REACTION> &annotated_reaction_list) {
  FILE* output = fopen(filename, "w");
  fprintf(output, "%s\t%s\t%s\n", "REACTION", "ANNOTATION", "GENE_PROBABILITY");
  for(int i=0; i<annotated_reaction_list.size(); i++) {
    for(int j=0; j<annotated_reaction_list[i].annote.size(); j++) {
      fprintf(output, "%s\t%s\t%4.3f\n", annotated_reaction_list[i].name, annotated_reaction_list[i].annote[j].genename.c_str(), annotated_reaction_list[i].annote[j].probability);
    }
  }
  fclose(output);
}
