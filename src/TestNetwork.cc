/* Series of functions that creates a small test network that looks like this:

          0             2
A  ------------> C   -------> D
      /                3     /
B----/       /-> G ---------/
        1   /
E--------------> F

Here the shortest path from A to D is to take A/B -> C -> D */

#include "TestNetwork.h"
#include <vector>
#include <map>
#include <stdio.h>

PROBLEM makeTestNetwork() {  
  vector<METABOLITE> metList;
  vector<REACTION> rxnList;
  vector<STOICH> stoichList;
  METABOLITE tmpMet;
  REACTION tmpRxn;
  STOICH tmpStoich;

  /* Construct metabolites */
  tmpMet = makeMetabolite(0, (char*)"A", 1, 0, 0, 0, 0, 0, 1.0f);   metList.push_back(tmpMet);
  tmpMet = makeMetabolite(1, (char*)"B", 1, 0, 0, 0, 0, 0, 1.0f);   metList.push_back(tmpMet);
  tmpMet = makeMetabolite(2, (char*)"C", 0, 0, 0, 0, 0, 0, 1.0f);   metList.push_back(tmpMet);
  tmpMet = makeMetabolite(3, (char*)"D", 0, 1, 0, 0, 0, 0, 1.0f);   metList.push_back(tmpMet);
  tmpMet = makeMetabolite(4, (char*)"E", 1, 0, 0, 0, 0, 0, 1.0f);   metList.push_back(tmpMet);
  tmpMet = makeMetabolite(5, (char*)"F", 0, 0, 0, 0, 0, 0, 1.0f);   metList.push_back(tmpMet);
  tmpMet = makeMetabolite(6, (char*)"G", 0, 0, 0, 0, 0, 0, 1.0f);   metList.push_back(tmpMet);

  /* Construct reactions; note that I don't bother with secondary vs not here but we can add that */

  /* Construct reaction 0 (A + B --> C) with likelihood 0 */
  stoichList.clear();
  tmpStoich = makeStoich((char*)"A", -1, 0); stoichList.push_back(tmpStoich);
  tmpStoich = makeStoich((char*)"B", -1, 1); stoichList.push_back(tmpStoich);
  tmpStoich = makeStoich((char*)"C", 1, 2);  stoichList.push_back(tmpStoich);
  tmpRxn = makeReaction(0, (char*)"RXN0", 0.0f, 0, 0, 0.0f, 1, 1, 1000.0f,stoichList);
  rxnList.push_back(tmpRxn);

  /* Construct reaction 1 (E --> G + F) with likelihood 1 */
  stoichList.clear();
  tmpStoich = makeStoich((char*)"E", -1, 4); stoichList.push_back(tmpStoich);
  tmpStoich = makeStoich((char*)"F", 1, 5);  stoichList.push_back(tmpStoich);
  tmpStoich = makeStoich((char*)"G", 1, 6);  stoichList.push_back(tmpStoich);
  tmpRxn = makeReaction(0, (char*)"RXN1", 0.0f, 0, 0, 1.0f, 1, 1, 1000.0f,stoichList);
  rxnList.push_back(tmpRxn);

  /* Construct reaction 2 (C --> D) with likelihood 2 */
  stoichList.clear();
  tmpStoich = makeStoich((char*)"C", -1, 2); stoichList.push_back(tmpStoich);
  tmpStoich = makeStoich((char*)"D", 1, 3);  stoichList.push_back(tmpStoich);
  tmpRxn = makeReaction(0, (char*)"RXN2", 0.0f, 0, 0, 2.0f, 1, 1, 1000.0f,stoichList);
  rxnList.push_back(tmpRxn);
  
  /* Construct reaction 3 (G --> D) with likelihood 3 */
  stoichList.clear();
  tmpStoich = makeStoich((char*)"G", -1, 6);  stoichList.push_back(tmpStoich);
  tmpStoich = makeStoich((char*)"D", 1, 3);  stoichList.push_back(tmpStoich);
  tmpRxn = makeReaction(0, (char*)"RXN3", 0.0f, 0, 0, 3.0f, 1, 1, 1000.0f,stoichList);
  rxnList.push_back(tmpRxn);

  RXNSPACE rxnspace(rxnList);
  METSPACE metspace(metList);

  /* Compute things:
     1. rxnsInvoled_nosec */
  calcMetRxnRelations_nosec(rxnspace, metspace);

  PROBLEM result;
  result.fullrxns = rxnspace;
  result.synrxns = rxnspace;

  result.metabolites = metspace;
}

METABOLITE makeMetabolite(int id, char* name, int input, int output, int biomass, 
			  int secondary_lone, int secondary_pair, int noncentral, 
			  double modifier) {

  METABOLITE met;
  met.id = id;
  sprintf(met.name, "%s", name);
  met.input = input;
  met.output = output;
  met.biomass = biomass;
  met.secondary_lone = secondary_lone;
  met.secondary_pair = secondary_pair;
  met.noncentral = noncentral;
  met.modifier = modifier;

  return met;
}

REACTION makeReaction(int id, char* name, double gibbs, double isExchange, double transporter, 
		      double init_likelihood, double init_reversible, int present, 
		      double flux_bound, vector<STOICH> &stoichList ) { 
  REACTION rxn;

  sprintf(rxn.name, "%s", name);

  rxn.gibbs = gibbs;
  rxn.isExchange = isExchange;
  rxn.transporter = transporter;
  rxn.init_likelihood = init_likelihood;
  rxn.current_likelihood = init_likelihood;
  rxn.init_reversible = init_reversible;
  rxn.net_reversible = init_reversible; /*Note - Dijkstras uses this one */
  rxn.present = present;
  rxn.flux_bound = flux_bound;
  rxn.stoich = stoichList;
  rxn.stoich_part = stoichList;

  return rxn;

}

STOICH makeStoich(char* metName, double stoichCoeff, int metId) {
  STOICH stoich;

  stoich.met_id = metId;
  stoich.rxn_coeff = stoichCoeff;
  sprintf(stoich.met_name, "%s", metName);

  return stoich;
}
