#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <libxml/parser.h>
#include <libxml/xmlmemory.h>

#include "DataStructures.h"
#include "MyConstants.h"
#include "pathUtils.h"
#include "Printers.h"
#include "XML_loader.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>

using std::vector;
using std::map;
using std::string;

void parseMETABOLITE(xmlDocPtr doc, xmlNodePtr cur, METSPACE &metspace);
void parseSTOICH(xmlDocPtr doc, xmlNodePtr cur, vector<STOICH> &stoich);
void parseREACTION(xmlDocPtr doc, xmlNodePtr cur, RXNSPACE &reaction);
void parseDoc(char *docname, RXNSPACE &rxnspace, METSPACE &metspace);
void parseMETAB(xmlDocPtr doc, xmlNodePtr cur, vector<MEDIA> &media);
void parseGROWTH(xmlDocPtr doc, xmlNodePtr cur, vector<GROWTH> &growth);
void parseMUTATION(xmlDocPtr doc, xmlNodePtr cur, vector<KNOCKOUT> &media);
void parseANNOTE (xmlDocPtr doc, xmlNodePtr cur, vector<ANNOTATION> &annote);

/* Parse both files (growth and rxn/metabolite) and do preliminary preparations 
(remove secondaries, calculate mets in each reaction, ...) */ 
void parseALL(char* file1, char* file2, PROBLEM &ProblemSpace){

  RXNSPACE &fullrxns = ProblemSpace.fullrxns;
  METSPACE &metspace = ProblemSpace.metabolites;
  vector<GROWTH> &growth = ProblemSpace.growth;

  parseDoc(file1, fullrxns, metspace);
  parseData(file2, growth);

  /* Mostly for reference, copy metabolite names over to the STOICHs from the METABOLITEs */
  addMetNameToStoichs(ProblemSpace.fullrxns.rxns, ProblemSpace.metabolites);
  addMetNameToStoichs(ProblemSpace.growth, ProblemSpace.metabolites);

  /* Now that the names are added... lets check and make sure we're not putting in insane inputs 
   and make sure that we actually have all the required special metabolites from DEBUGFLAGS in our metabolite list */
  checkConsistency(ProblemSpace);

  /* Remove secondaries and make magic bridges between cofactor pairs */
  Load_Stoic_Part(ProblemSpace.fullrxns.rxns, ProblemSpace.metabolites.mets);

  /* This is used to add "inputs" and "outputs" to our metabolite list, so that the path visualizer deals with them
     correctly.  For now, I only bother with the first growth condition, because I don't know how graphviz would respond to something having two different
     commands for how to format the circles. I imagine it would not be happy.*/
  vector<int> dum, dum2;
  loadInputsOutputs(ProblemSpace.growth, 0, ProblemSpace.metabolites, dum, dum2);

  /* Identify reactions that name "free" metabolites, i.e. reactions that are 1) NOT exchanges, and
     2) Only have secondary metabolites on one side. 
   
     If the metabolite is "free" the "freeMakeFlag" is set to 1. */
  identifyFreeReactions(ProblemSpace.fullrxns.rxns);

  //  addBridgeMetabolites(ProblemSpace);

  return;
}

void parseMETABOLITE (xmlDocPtr doc, xmlNodePtr cur, METSPACE &metspace) {
  xmlChar *key;
  METABOLITE tempm;
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"met_id"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("id: %s\n", key); */
      tempm.id = atoi((char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"charge"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("charge: %s\n", key); */
      tempm.charge = atoi((char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"name"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("name: %s\n", key); */
      strcpy(tempm.name,(char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"chemform"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("chemform: %s\n", key); */
      strcpy(tempm.chemform,(char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"secondary"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("secondary_lone: %s\n", key); */
      tempm.secondary_lone = atoi((char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"cofactor"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("secondary_pair: %s\n", key); */
      tempm.secondary_pair.push_back(atoi((char*)key));
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"noncentral"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("secondary_pair: %s\n", key); */
      tempm.noncentral = atoi((char*)key);
      xmlFree(key);
    }
    cur = cur->next;
  }

  metspace.addMetabolite(tempm);
  return;
}

void parseSTOICH (xmlDocPtr doc, xmlNodePtr cur, vector<STOICH> &stoich){
  xmlChar *key;
  cur = cur->xmlChildrenNode;
  STOICH temps;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"met_id"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("met_id: %s\n", key); */
      temps.met_id = atoi((char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"stoich"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("rxn_coeff: %s\n", key); */
      temps.rxn_coeff = atof((char*)key);
      xmlFree(key);
    }
    cur = cur->next;
  }

  if(temps.met_id == -1) { 
    printf("WARNING: Skipped over a STOICH with no met_id present in it - this is very likely an XML error!\n"); }
  else {      stoich.push_back(temps);  }
  return;
}

void parseANNOTE (xmlDocPtr doc, xmlNodePtr cur, vector<ANNOTATION> &annote){
  xmlChar *key;
  cur = cur->xmlChildrenNode;
  ANNOTATION temps;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"gene"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("met_id: %s\n", key); */
      temps.genename = (char*)key;
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"prob"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("rxn_coeff: %s\n", key); */
      temps.probability = atof((char*)key);
      xmlFree(key);
    }
    cur = cur->next;
  }
  annote.push_back(temps);
  return;
}

void parseREACTION (xmlDocPtr doc, xmlNodePtr cur, RXNSPACE &rxnspace) {
  xmlChar *key;
  REACTION tempr;
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"id"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("id: %s\n", key); */
      tempr.id = atoi((char*)key);
	xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"name"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("name: %s\n", key); */ 
      strcpy(tempr.name,(char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"s"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      parseSTOICH (doc, cur, tempr.stoich);
      //tempr.numMets += 1;
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"exchange"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      tempr.isExchange = atoi((char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"rev"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      tempr.init_reversible = atoi((char*)key);
      /* For now I set the net_reversible equal to init_reversible */
      tempr.net_reversible = tempr.init_reversible;
      if(tempr.init_reversible < 0)  { tempr.ub = 0.0f; } else { tempr.ub = 1000.0f; }
      if(tempr.init_reversible > 0)  { tempr.lb  = 0.0f;} else { tempr.lb = -1000.0f;}
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"transport"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      tempr.transporter = atoi((char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"likelihood"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      tempr.init_likelihood = atof((char*)key);
      /* Set current likelihood = net likelihood BUT if it is one of the special switches we have to do special things with them in Dijkstra's - that program will handle it */
      tempr.current_likelihood = tempr.init_likelihood;
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"annote"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      parseANNOTE(doc, cur, tempr.annote);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"synthesis"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      tempr.synthesis = atoi((char*)key);
      xmlFree(key);
    }
    cur = cur->next;
  }
  /* MATT CHANGE (6-27-11)
     Only bother including this reaction if there is anything present in S. Otherwise it's a dud... */
  if(!tempr.stoich.empty()) {
    rxnspace.addReaction(tempr);
  }
  return;
}

void parseDoc(char *docname, RXNSPACE &rxnspace, METSPACE &metspace) {
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile(docname);
  if (doc == NULL ) {
    fprintf(stderr,"Document not parsed successfully. \n");
    assert(doc != NULL);
  }
  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr,"empty document\n");
    xmlFreeDoc(doc);
    assert(cur != NULL);
  }
  if (xmlStrcmp(cur->name, (const xmlChar *) "model")) {
    fprintf(stderr,"document of the wrong type, root node != story");
    xmlFreeDoc(doc);
    assert(false);
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"metabolite"))){
      parseMETABOLITE (doc, cur, metspace);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"reaction"))){
      parseREACTION (doc, cur, rxnspace);
    }
    cur = cur->next;
  }
  xmlFreeDoc(doc);
  return;
}

/* Parse MEDIA */
void parseMETAB (xmlDocPtr doc, xmlNodePtr cur, vector<MEDIA> &media) {
  xmlChar *key;
  MEDIA tempme;
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"met_id"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      //printf("id: %s\n", key);
      tempme.id = atoi((char*)key);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"rate"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("name: %s\n", key); */
      tempme.rate = atof((char*)key);
      xmlFree(key);
    }
    if ( (!xmlStrcmp(cur->name, (const xmlChar *)"name"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      strcpy(tempme.name, (char*)key);
      xmlFree(key);
    }
    cur = cur->next;
  }

  /* Note - this could print garbage if it also had no name */
  if(tempme.id == -1) {  printf("WARNING: Metabolite %s had no met ID in the MEDIA file and therefore was skipped\n", tempme.name); }
  else {  media.push_back(tempme);  }

  return;
}

void parseMUTATION (xmlDocPtr doc, xmlNodePtr cur, vector<KNOCKOUT> &kos) {
  xmlChar *key;
  KNOCKOUT tempko;
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, (const xmlChar *)"gene_id")) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* How will we get the IDs? */
      tempko.genename = (char*)key;
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"act_coef"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      tempko.act_coef = atof((char*)key);
      xmlFree(key);
    }
    cur = cur->next;
  }
  kos.push_back(tempko);
  return;
}

void parseGROWTH (xmlDocPtr doc, xmlNodePtr cur, vector<GROWTH> &growth) {
  xmlChar *key;
  xmlNodePtr cur2;
  GROWTH tempg;
  //printf("in parseGROWTH\n");fflush(stdout);
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"media"))) {
      cur2 = cur->xmlChildrenNode;
      //printf("entering MEDIA\n");fflush(stdout);
      while (cur2 != NULL) {
	if ((!xmlStrcmp(cur2->name, (const xmlChar *)"metab"))) {
	  //printf("\tabout to parseMETAB...");fflush(stdout);
	  parseMETAB(doc,cur2,tempg.media);
	  //printf(" done\n");fflush(stdout);
	}
	cur2 = cur2->next;
      }
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"byproducts"))) {
      cur2 = cur->xmlChildrenNode;
      //printf("entering BYPRODUCTS\n");fflush(stdout);
      while (cur2 != NULL) {
	if ((!xmlStrcmp(cur2->name, (const xmlChar *)"metab"))) {
	  //printf("\tabout to parseMETAB...");fflush(stdout);
	  parseMETAB(doc,cur2,tempg.byproduct);
	  //printf(" done\n");fflush(stdout);
	}
	cur2 = cur2->next;
      }
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"biomass"))) {
      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL) {
	if ((!xmlStrcmp(cur2->name, (const xmlChar *)"s"))) {
	  parseSTOICH (doc, cur2, tempg.biomass);
	}
	cur2 = cur2->next;
      }
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"geneko"))) {
	  printf("about to parseMUTATION...");fflush(stdout);
	  parseMUTATION(doc, cur, tempg.mutation);
	  printf(" done\n");fflush(stdout);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"growth_rate"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      /* printf("gibbs: %s\n", key); */
      if(strcmp((char*)key,"NaN")==0 ||
	 strcmp((char*)key,"nan")==0 ||
	 strcmp((char*)key,"NAN")==0 ||
	 strcmp((char*)key,"Nan")==0 ) { tempg.growth_rate = 0.0f;}
      else{ tempg.growth_rate = atof((char*)key);}
      xmlFree(key);
    }
    cur = cur->next;
  }
  growth.push_back(tempg);
  return;
}

void parseData(char *docname, vector<GROWTH> &growth) {
  xmlDocPtr doc;
  xmlNodePtr cur;
  doc = xmlParseFile(docname);
  if (doc == NULL ) {
    fprintf(stderr,"Document not parsed successfully. \n");
    return;
  }
  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr,"empty document\n");
    xmlFreeDoc(doc);
    return;
  }
  if (xmlStrcmp(cur->name, (const xmlChar *) "data")) {
    fprintf(stderr,"document of the wrong type, root node != data");
    xmlFreeDoc(doc);
    return;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"growth"))){
      parseGROWTH (doc, cur, growth);
    }
    cur = cur->next;
  }
  xmlFreeDoc(doc);
  return;
}

void identifyFreeReactions(vector<REACTION> &reactions) {
  for(int i=0; i<reactions.size(); i++) {

    /* Exclude exchange reactions */
    if(reactions[i].isExchange == 1) { continue; }

    bool hasReactant = false;
    bool hasProduct = false;

    for(int j=0; j<reactions[i].stoich_part.size(); j++) {
      if(reactions[i].stoich_part[j].rxn_coeff < 0) { 
	hasReactant = true;
      }
      if(reactions[i].stoich_part[j].rxn_coeff > 0) {
	hasProduct = true;
      }
    }

    /* Empty reactions (all secondary) */
    if(!hasProduct && !hasReactant) { continue; }

    /* Has both (so let Dijkstras try to find a path) */
    if(hasProduct && hasReactant) { continue; }

    switch(reactions[i].init_reversible) {
    case -1:
      if(hasReactant) {
	reactions[i].freeMakeFlag = 1;
      }
      break;
    case 0:
      /* reversible: No matter which way, we consider it free */
      reactions[i].freeMakeFlag = 1;
      break;
    case 1:
      if(hasProduct) {
	reactions[i].freeMakeFlag = 1;
      }
      break;
    }
  }
}

/* Load stoich_part as stoichiometry without secondary metabolites 
ID of all secondaries now occurs internally in this function

1-16-12 - TEST of the "synthesis reaction" concept
Secondary_lones are not removed if they are considered "synthesis reactions"
for a given compound. 

FIXME: Do not remove things that we feed the cell - these can be assumed to be transported
and hence are legitimately "free".

*/
void Load_Stoic_Part(vector<REACTION> &reaction, 
		     const vector<METABOLITE> &metabolite){
  
  vector<int> singleSecondary_ID;
  vector<int> secondaryPair_ID1;
  vector<int> secondaryPair_ID2;

  for(int i=0;i<metabolite.size(); i++) {

    /* Singlet */
    if(metabolite[i].secondary_lone != 0) {
      singleSecondary_ID.push_back(metabolite[i].id);
    }

    /* We'll get duplicates doing it this way but that's actually good so we don't 
       need to check both directions explicitly */
    if(!metabolite[i].secondary_pair.empty()) {
      for(int j=0;j<metabolite[i].secondary_pair.size();j++){
	secondaryPair_ID1.push_back(metabolite[i].id);
	secondaryPair_ID2.push_back(metabolite[i].secondary_pair[j]);
      }
    }
  }
  
  for(int i=0;i<reaction.size();i++) {

    /* Build list of metabolites involved in the reaction */
    map<int, int> metId2StoichIdx;
    for(int j=0; j<reaction[i].stoich.size(); j++) {
      metId2StoichIdx[reaction[i].stoich[j].met_id] = j;
    }

    for(int j=0;j<reaction[i].stoich.size();j++) {
      STOICH tempPart;

      /* Test for double secondary and don't copy if it's part of a secondary pair in that reaction */
      if( secondaryPairPresent(secondaryPair_ID1, secondaryPair_ID2, reaction[i], reaction[i].stoich[j].met_id, metId2StoichIdx) ) { 
	continue; 
      }

      /* Single secondary - don't copy (do this second to allow for cases with both single and double secondaries) */
      vector<int>::iterator iter = find(singleSecondary_ID.begin(), singleSecondary_ID.end(), reaction[i].stoich[j].met_id);
      if(iter != singleSecondary_ID.end()) {
	/* NEW 1-17-2012 - skip over it if it is considered a synthesis reaction for that metabolite */
	if(_db.TEST_SYNTHESIS) {
	  if(reaction[i].synthesis != reaction[i].stoich[j].met_id) { continue; }
	} else {
	  continue; 
	}
      };
      
      /* Otherwise just copy it */
      tempPart.met_id = reaction[i].stoich[j].met_id;
      tempPart.rxn_coeff = reaction[i].stoich[j].rxn_coeff;
      sprintf(tempPart.met_name, "%s", reaction[i].stoich[j].met_name);
      reaction[i].stoich_part.push_back(tempPart);
    }
  }
}

/* Find out if metId is part of a cofactor pair (two values with the same index in secondaryList1 and secondaryList2 
 are cofactor pairs) - note that I only export the cofactor pair in ONE of the pair of metabolites now, in order to
 simplify making of magic bridges. */
bool secondaryPairPresent(const vector<int> &secondaryList1, const vector<int> &secondaryList2, const REACTION &reaction, int metId, map<int, int> &metId2StoichIdx) {

  vector<int>::const_iterator it = secondaryList1.begin();
  int idx;
  
  while(true) {

    /* We only worry about secondaryList1 because every pair is present in both set of metabolites.
       it has to be this way for downstream stuff tow ork. */
    it = find(it, secondaryList1.end(), metId);
    if(it==secondaryList1.end()) { break; }
    idx = distance(secondaryList1.begin(), it);

    map<int, int>::const_iterator metIter1, metIter2;

    metIter1 = metId2StoichIdx.find(secondaryList1[idx]);
    metIter2 = metId2StoichIdx.find(secondaryList2[idx]);

    /* Check that both mets are in the reaction and on oppposite sides */
    if(metIter1!=metId2StoichIdx.end() && metIter2!=metId2StoichIdx.end()) {
      if( ( reaction.stoich[metIter1->second].rxn_coeff *
	    reaction.stoich[metIter2->second].rxn_coeff ) < 0 ) {
	return true;
      }
    }
    it++;
  }

  return false;
}

/* Add metabolite names to the stoich structure */
void addMetNameToStoichs(vector<REACTION> &reaction, METSPACE &metspace) {

  for(int i=0;i<reaction.size();i++) {
    REACTION tmpRxn = reaction[i];
    /* Update stoich */
    for(int j=0;j<tmpRxn.stoich.size();j++) {
      sprintf(reaction[i].stoich[j].met_name, "%s",  metspace.metFromId(tmpRxn.stoich[j].met_id).name);
    }
    /* Update stoich_part */
    for(int j=0;j<tmpRxn.stoich_part.size();j++) {
      sprintf(reaction[i].stoich_part[j].met_name, "%s", metspace.metFromId(tmpRxn.stoich_part[j].met_id).name);
    }
  }
  return;
}

/* Add metabolite names to the stoich structure */
void addMetNameToStoichs(vector<GROWTH> &growth, METSPACE &metspace) {

  for(int j=0;j<growth.size();j++){
    for(int i=0;i<growth[j].biomass.size();i++) {
      sprintf(growth[j].biomass[i].met_name,"%s", metspace.metFromId(growth[j].biomass[i].met_id).name);
    }
  }
  return;
}

/* Fill the inputIds and outputIds vectors with the input and biomass (output) IDs declared in growth[growthIdx]
Also puts whether the metabolite is an input/output or not into the appropriate places in "metabolite".
The latter must be done for visualization only.
*/
void loadInputsOutputs(vector<GROWTH> &growth, int growthIdx, METSPACE &metspace, vector<int> &inputIds, vector<int> &outputIds) {

  if( (inputIds.size() != 0) || (outputIds.size() != 0 ) ) {
    printf("WARNING: either inputIds or outputIds provided to loadInputsOutputs was not empty! The function expects empty vectors and will clear them before filling back up.");
    inputIds.clear();
    outputIds.clear();
  }
  /* Load inputs and outputs from growth[growthIdx] */
  for(int j=0;j<growth[growthIdx].media.size();j++) {
    inputIds.push_back(growth[growthIdx].media[j].id);
    /* For visualization */
    metspace.metPtrFromId(growth[growthIdx].media[j].id)->input=1;
  }
  for(int j=0;j<growth[growthIdx].biomass.size();j++) {
    outputIds.push_back(growth[growthIdx].biomass[j].met_id);
    metspace.metPtrFromId(growth[growthIdx].biomass[j].met_id)->output=1;
  }
  /* For visualization, highlight outputs as "outputs" in addition to the biomass component */
  for(int j=0;j<growth[growthIdx].byproduct.size();j++) {
    metspace.metPtrFromId(growth[growthIdx].byproduct[j].id)->input=1;
  }
}

/* Makes "magic bridges" between cofactors so that Dijkstras can find paths to synthesize them (since the secondary metabolites are also potentially in the biomass equation, this will be necessary). The magic bridges have a likelihood of 1 before adjusting likelihoods and an ID that is 4000xx where xx is incremented for each cofactor pair. The reaction ID will be that of either the reactant or product that it starts with, plus 400,000. */
void makeMagicBridges(PROBLEM &ProblemSpace) {

  vector<REACTION> &reaction = ProblemSpace.synrxns.rxns;
  vector<METABOLITE> &metabolite = ProblemSpace.metabolites.mets;
  map<int, int> &metId2Idx = ProblemSpace.metabolites.Ids2Idx;
  
  /* Search for cofactor pairs */
  int ctr, flag;
  flag = _db.MAGICBRIDGEFACTOR;
  ctr = 0;

  vector< vector<int> > donePairs;
  

  for(int i=0;i<metabolite.size();i++) {
    REACTION tmp_rxn;
    STOICH tmp_stoich;
    vector<STOICH> tmp_stoich_vec;

    if(metabolite[i].secondary_pair.empty()) {
      continue;
    }
    for(int j=0;j<metabolite[i].secondary_pair.size();j++){

      vector<int> onePair(2, 0);
      onePair[0] = metabolite[i].id;
      onePair[1] = metabolite[i].secondary_pair[j];
      if(pairDone(donePairs, onePair)) { continue; }
      donePairs.push_back(onePair);

      /* Add a reaction from A to B (irreversible)  */
      tmp_rxn.id = ctr + flag;
      sprintf(tmp_rxn.name, "MAGICBRIDGE_%s_%s",  metabolite[i].name, 
	      metabolite[metId2Idx[metabolite[i].secondary_pair[j]]].name);
      tmp_rxn.init_reversible = 0; /* Reversible */
      tmp_rxn.net_reversible = tmp_rxn.init_reversible; 
      tmp_rxn.init_likelihood = 1.0f; /* Likely (net_likelihood will be gotten later) */
      tmp_rxn.current_likelihood = 0.1f; /* In case someone doesnt reset the likelihoods we'd rather this thing not die a horrible death */

      /* A */
      tmp_stoich.met_id = metabolite[i].id;
      tmp_stoich.rxn_coeff = -1;
      sprintf(tmp_stoich.met_name, "%s", metabolite[i].name);
      tmp_stoich_vec.push_back(tmp_stoich);
      /* B */
      tmp_stoich.reset();
      tmp_stoich.rxn_coeff = 1;
      tmp_stoich.met_id = metabolite[i].secondary_pair[j];
      sprintf(tmp_stoich.met_name, "%s", ProblemSpace.metabolites.metFromId(metabolite[i].secondary_pair[j]).name);
      tmp_stoich_vec.push_back(tmp_stoich);
      
      tmp_rxn.stoich = tmp_stoich_vec;
      // Allows it to actually be used in Dijkstras
      tmp_rxn.stoich_part = tmp_stoich_vec;
      
      reaction.push_back(tmp_rxn);
      ctr++;
      //printREACTIONintermediates(tmp_rxn, 1);
    }
  }
}

/* TRUE if onePair is one of the pairs within pairList (also would return true if onePair has two identical values) */
bool pairDone(const vector<vector<int> > &pairList, const vector<int> &onePair) {
  for(int i=0; i<pairList.size(); i++) {
    if(onePair[0] == pairList[i][0] || 
       onePair[0] == pairList[i][1] ) {
      if(onePair[1] == pairList[i][0] || 
	 onePair[1] == pairList[i][1] ) {
	return true;
      }
    }
  }
  return false;
}

/* Ensure consistency of NGAM with what we expect (ATP consumption), and
   make sure there is capacity to use GAM in the biomass equation (add ATP / ADP / Pi / H / H2O if they aren't already in there
   with 0 stoichiometry - and then we'll make a fucntion to modify them in series) */
void setUpMaintenanceReactions(PROBLEM &ProblemSpace) {
  DEBUGFLAGS db;

  /********* NGAM **************/
  int ATPM_id = Name2Ids(ProblemSpace.fullrxns.rxns, db.ATPM_name);
  if(ATPM_id == -1) {
    /* FIXME: What should we use as the ID? */
    /* For now I just make this an error */
    printf("ERROR: Specified Non-growth associated maintenance (ATPM) reaction not found!\n");
    assert(false);
  }
  /* Check that ATPM contains ATP, ADP, H2O, and PI (and possibly H) and nothing else */
  vector<int> mustIds;
  mustIds.push_back(Name2Ids(ProblemSpace.metabolites.mets, db.ATP_name));
  mustIds.push_back(Name2Ids(ProblemSpace.metabolites.mets, db.ADP_name));
  mustIds.push_back(Name2Ids(ProblemSpace.metabolites.mets, db.PI_name));
  mustIds.push_back(Name2Ids(ProblemSpace.metabolites.mets, db.H2O_name));
  /* H is optional because different databases may charge balance differently */
  int optId = Name2Ids(ProblemSpace.metabolites.mets, db.H_name);
  int numH(0);

  /* This is how we keep track of ones that HAVE to be there */
  vector<bool> mustOK(mustIds.size(), false);
  
  REACTION ATPM = ProblemSpace.fullrxns.rxnFromId(ATPM_id);
  for(int i=0; i<ATPM.stoich.size(); i++) {
    bool OK(false);
    for(int j=0; j<mustIds.size(); j++) {
      if(mustIds[j] == ATPM.stoich[i].met_id) { 
	OK = true; 
	mustOK[j] = true;
      }
    }
    if(optId == ATPM.stoich[i].met_id) { 
      OK = true; 
      numH = ATPM.stoich[i].rxn_coeff;
    }
    /* Not a recgnized reaction */
    if(!OK) { printf("ERROR: ATP maintenance reaction %s specified by DEBUGFLAGS is not the expected ATP consumption reaction!\n", db.ATPM_name); assert(false); }
  }

  /* Add NGAM to the biomass equation (take the number of H's from the provided ATPM reaction) 
   Any of those metabolites not present are added with a coefficient of 0 
  Note - I don't bother guessing the ATP maintenance value from teh biomass equation. Instead I just
  start at 0 and adjust it from there in the optimizer... 
  I'll just have to make sure everything stays on the correct side of the equation (i.e. no ATP generation!) */
  REACTION biomass = ProblemSpace.fullrxns.rxnFromId(db.BIOMASS);
  vector<bool> requiredPresent(mustIds.size(), false);
  bool hPresent(false);
  for(int i=0; i<biomass.stoich.size(); i++) { 
    for(int j=0; j<mustIds.size(); j++) {
      if(biomass.stoich[i].met_id == mustIds[j]) { requiredPresent[j] = true;  }
    }
    if(biomass.stoich[i].met_id == optId) { hPresent = true; }
  }
  for(int j=0; j<mustIds.size(); j++) {
    if(!requiredPresent[j]) { 
      STOICH tmpStoich;
      METABOLITE tmpMet = ProblemSpace.metabolites.metFromId(mustIds[j]);
      printf("WARNING: Biomass seems to be missing a metabolite %s that is required for NGAM calculations - adding to the biomass equation\n", tmpMet.name);
      sprintf(tmpStoich.met_name, "%s", tmpMet.name);
      tmpStoich.met_id = tmpMet.id;
      tmpStoich.rxn_coeff = 0.0f;
    }
  }

  return;

}

/* Consistency check for the names of metabolites from the two different XML files 
 This is to at least catch headaches in the act of happening rather than after they kill you 
 This is also the main reason I printed the names into the input files. */
void checkConsistency(const PROBLEM &ProblemSpace) {

  const METSPACE &metspace = ProblemSpace.metabolites;
  const RXNSPACE &rxnspace = ProblemSpace.fullrxns;
  const vector<GROWTH> &growth = ProblemSpace.growth;

  DEBUGFLAGS db;

  printf("Checking consistency of the two input XML files with each other...\n");

  for(int i=0; i<growth.size(); i++) {
    /* Check biomass */
    for(int j=0; j<growth[i].biomass.size(); j++) {
      int growthId = growth[i].biomass[j].met_id;
      if(!metspace.idIn(growthId)) {
	printf("ERROR: The provided InputData and Database XML files have inconsistent IDs, this program will now terminate\n");
	printf("Inconsistent name: %s in the biomass had no corresponding ID in the metabolite file \n", growth[i].biomass[j].met_name);
	assert(false);
      }
      if( strcmp( metspace.metFromId(growthId).name , growth[i].biomass[j].met_name ) != 0 ) {
	printf("ERROR: The provided InputData and Database XML files have inconsistent IDs, this program will now terminate\n");
	printf("Inconsistent name: %s in the biomass did not correspond with %s in the metabolite file despite both having the same ID \n", growth[i].biomass[j].met_name, metspace.metFromId(growthId).name);
	assert(false);
      }
    }

    /* Check byproducts */
    for(int j=0; j<growth[i].byproduct.size(); j++) {
      int growthId = growth[i].byproduct[j].id;
      if(!metspace.idIn(growthId)) {
	printf("ERROR: The provided InputData and Database XML files have inconsistent IDs, this program will now terminate\n");
	printf("Inconsistent name: %s in the byproducts had no corresponding ID in the metabolite file \n", growth[i].byproduct[j].name);
	assert(false);
      }
      if( strcmp( metspace.metFromId(growthId).name , growth[i].byproduct[j].name ) != 0 ) {
	printf("ERROR: The provided InputData and Database XML files have inconsistent IDs, this program will now terminate\n");
	printf("Inconsistent name: %s in the byproducts did not correspond with %s in the metabolite file despite both having the same ID \n", growth[i].byproduct[j].name, metspace.metFromId(growthId).name);
	assert(false);
      }
    }

    /* Check media */
    for(int j=0; j<growth[i].media.size(); j++) {
      int growthId = growth[i].media[j].id;
      if(!metspace.idIn(growthId)) {
        printf("ERROR: The provided InputData and Database XML files have inconsistent IDs, this program will now terminate\n");
	printf("Inconsistent name: %s in the media had no corresponding ID in the metabolite file \n", growth[i].media[j].name);
        assert(false);
      }
      if( strcmp( metspace.metFromId(growthId).name , growth[i].media[j].name ) != 0 ) {
        printf("ERROR: The provided InputData and Database XML files have inconsistent IDs, this program will now terminate\n");
        printf("Inconsistent name: %s in the media did not correspond with %s in the metabolite file despite both having the same ID \n", growth[i].media[j].name, metspace.metFromId(growthId).name);
        assert(false);
      }
    }
  }

  /* Test special metabolite and reaction names to make sure they exist */
  int h_id = Name2Ids(metspace.mets, db.H_name);
  if(h_id == -1) { printf("ERROR: DEBUGFLAGS claims metabolite %s is the name for H+ in your database, but it not found in the database\n", db.H_name); assert(false); }
  int h_e_id = Name2Ids(metspace.mets, db.H_plus_E);
  if(h_id == -1) { printf("ERROR: DEBUGFLAGS claims metabolite %s is the name for H+[e] in your database, but it not found in the database\n", db.H_plus_E); assert(false); }
  int na_id = Name2Ids(metspace.mets, db.Na_name);
  if(h_id == -1) { printf("ERROR: DEBUGFLAGS claims metabolite %s is the name for sodium (Na) in your database, but it not found in the database\n", db.Na_name); assert(false); }
  int na_e_id = Name2Ids(metspace.mets, db.Na_plus_E);
  if(h_id == -1) { printf("ERROR: DEBUGFLAGS claims metabolite %s is the name for external sodium (Na) in your database, but it not found in the database\n", db.Na_plus_E); assert(false); }
  int atp_id = Name2Ids(metspace.mets, db.ATP_name);
  if(h_id == -1) { printf("ERROR: DEBUGFLAGS claims metabolite %s is the name for ATP in your database, but it not found in the database\n", db.ATP_name); assert(false); }
  int adp_id = Name2Ids(metspace.mets, db.ADP_name);
  if(h_id == -1) { printf("ERROR: DEBUGFLAGS claims metabolite %s is the name for ADP in your database, but it not found in the database\n", db.ADP_name); assert(false); }
  int pi_id = Name2Ids(metspace.mets, db.PI_name);
  if(h_id == -1) { printf("ERROR: DEBUGFLAGS claims metabolite %s is the name for phosphate (pi) in your database, but it not found in the database\n", db.PI_name); assert(false); }
  int h2o_id = Name2Ids(metspace.mets, db.H2O_name);
  if(h_id == -1) { printf("ERROR: DEBUGFLAGS claims metabolite %s is the name for water (h2o) in your database, but it not found in the database\n", db.H2O_name); assert(false); }

  printf("Consistency check done - XML files have consistent names and IDs and both inputs are compatible with the DEBUGFLAGS information \n");

}
