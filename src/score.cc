#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<assert.h>
#include<map>
#include<limits.h>
#include<gsl/gsl_min.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multimin.h>

#include"DataStructures.h"
#include"Exchanges.h"
#include"Grow.h"
#include"genericLinprog.h"
#include"MyConstants.h"

using namespace std;

ANSWER overwriteRXNS(vector<REACTION> &ETC_adjusted, const ANSWER &ans){
  ANSWER result(ans);
  for(int i=0;i<ETC_adjusted.size();i++){
    int rxn_idx = result.reactions.idxFromId(ETC_adjusted[i].id);
    REACTION &rxn_ptr = result.reactions.rxns[rxn_idx];
    assert(ETC_adjusted[i].stoich.size()==rxn_ptr.stoich.size());
    for(int j=0;j<rxn_ptr.stoich.size();j++){
      rxn_ptr.stoich[j] = ETC_adjusted[i].stoich[j];
    }
  } 
  return result;
}

map<int,double> getSecretionRates(const ANSWER &ans1, const GROWTH &growth1, 
				 const vector<double> &fba_solution){
  map<int,double> result; //met_id to secretion rate mapping

  //indexed by the order they appear in growth
  for(int i=0;i<growth1.byproduct.size();i++){
    int met_id = growth1.byproduct[i].id;
    int rxn_id = FindExchange4Metabolite(ans1.reactions.rxns,met_id);
    int rxn_idx= ans1.reactions.idxFromId(rxn_id);
    result[met_id] = fba_solution[rxn_idx];
  }

  return result;
}

//This weights growths equal to secretion rates. Only one growth at a time. 
double evalProtocal1(const vector<double> sim, const vector<double> exp){
  assert(sim.size()==exp.size());

  //measure weights so that the fitting on the absolute scale is roughly equal between
  //growth and secretion rates
  double secret_weight = 0.0f;
  double denominator = abs(sim[0]) + abs(exp[0]);
  if(denominator>100.0f || isinf(denominator)){denominator = 100.0f;}
  double growth_weight = 1.0f / (abs(denominator) * 0.5f);    //assumes growth is the 0th index
  for(int i=1;i<sim.size();i++){
    secret_weight += abs(sim[i]) + abs(exp[i]);
  }
  secret_weight *= 0.5f;
  if(secret_weight > 0.01){
    secret_weight = 1.0f/secret_weight;}
  else{
    secret_weight = 100;}

  double diffs[sim.size()];
  for(int i=0;i<sim.size();i++){
    diffs[i] = abs(sim[i] - exp[i]);
    if(exp[i] > 999.9f){
      if(sim[i] < 0.0001 && exp[i] < 0.0001){ 
	diffs[i] = 0.0;}
      else{ diffs[i] = 1.0;}
    }
    //printf("diffs %f\n",diffs[i]);
  }
  double weighted_diffs[sim.size()];
  //printf("growth_weight: %f\n",growth_weight);
  weighted_diffs[0] = diffs[0] * growth_weight;
  for(int i=1;i<sim.size();i++){
    weighted_diffs[i] = diffs[i] * secret_weight;}

  double score = 0.0f; 
  for(int i=0;i<sim.size();i++){
    score += weighted_diffs[i]; 
    //printf("weighted_diffs %f\n",weighted_diffs[i]);
}

  return score/2.0f;
}

double measureScore(SCORE1 &score1,vector<ANSWER> &ans, vector<GROWTH> &growth){

  vector<double> scores2;
  for(int i=0;i<ans.size();i++){
    //printf("modifying NGAM\n"); 
    modifyNGAM(ans[i],score1.ngam);
    //printf("modifying GAM\n"); 
    modifyGAM(ans[i],score1.gam);
    //printf("done\n\n");
    
    vector<vector<double> > sim_growth_rates, exp_growth_rates;
    vector<double> scores;
    //load the sim vs exp data
    for(int j=0;j<growth.size();j++){
      //printf("overwriting reactions\n");
      ANSWER temp_ans = overwriteRXNS(score1.ETC_adjusted, ans[i]);
      // printf("FBA solving\n");
      vector<double> g = FBA_SOLVE(temp_ans.reactions, temp_ans.metabolites);
      //growth rate
      score1.growthRate[j] = g[temp_ans.reactions.idxFromId(_db.BIOMASS)];
      //printf("Growth Rate: %f\n",g[0]);
      //printf("saving\n");
      vector<double> temp_vec_double;
      sim_growth_rates.push_back(temp_vec_double);
      exp_growth_rates.push_back(temp_vec_double);
      sim_growth_rates[j].push_back(g[temp_ans.reactions.idxFromId(_db.BIOMASS)]);
      exp_growth_rates[j].push_back(growth[j].growth_rate);
      //byproducts
      //printf("getting secretion rates\n");
      map<int,double> fba_byproduct_result = getSecretionRates(temp_ans,growth[j],g);
      //printf("going into loop\n");
      for(int k=0;k<growth[j].byproduct.size();k++){
	int met_id = growth[j].byproduct[k].id;
	sim_growth_rates[j].push_back(fba_byproduct_result[met_id]);
	exp_growth_rates[j].push_back(growth[j].byproduct[k].rate);
      }
      if(0){
	printf("#%d\t",j);
	for(int k=0;k<exp_growth_rates[j].size();k++){
	  printf("%#8g\t", exp_growth_rates[j][k]);}
	printf("\n# \t");
	for(int k=0;k<sim_growth_rates[j].size();k++){
	  printf("%#8g\t", sim_growth_rates[j][k]);}
	printf("\n");
      }
      //growth rates equal to secretion rates in fitting weight
      double weight = (sim_growth_rates[j][0] + exp_growth_rates[j][0]) * 0.5;
      double temp = evalProtocal1(sim_growth_rates[j],exp_growth_rates[j]) * weight;
      scores.push_back(temp);
      score1.growthScore[j] = temp;
    }
    //printf("loaded the sim vs exp data\n");

    //printf("going into scoring loop\n");
    double score = 0.0f;
    for(int j=0;j<scores.size();j++){ score += scores[j];}
    //printf("done\n");
    score = score / (int)scores.size();
    //printf("%f\n",score);
    scores2.push_back(score);
  }

  double score = 0.0f;
  for(int i=0;i<scores2.size();i++){ score += scores2[i];}
  score = score / (int)scores2.size();
  score1.score = score;
  //printf("made score score \n");

  return score1.score;
}

class PARAMS{
public:
  vector<ANSWER> *ans;
  PROBLEM *ProblemSpace;
};

double my_f(const gsl_vector * x, void * params){
  PARAMS *P = (PARAMS *) params;
  SCORE1 score1;
  score1.ngam = gsl_vector_get(x,0); //size_t i is the type here...
  score1.gam = gsl_vector_get(x,1);
  double result = measureScore(score1,*P->ans,P->ProblemSpace->growth);
  return result;
}

void optimizeAMs(vector<ANSWER> &ans, PROBLEM &ProblemSpace, double init_AM){
  printf("optimizeAMs: 0\n");
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss;
  gsl_multimin_function minex_func;
  printf("optimizeAMs: 1\n");

  size_t iter = 0;
  int status;
  double size;
  printf("optimizeAMs: 2\n");

  gsl_vector *AMs = (gsl_vector *) gsl_vector_alloc(sizeof(gsl_vector) * 2);
  gsl_vector_set_all(AMs,init_AM);
  printf("optimizeAMs: 3\n");

  ss = gsl_vector_alloc(2);
  gsl_vector_set_all (ss,1.0);
  printf("optimizeAMs: 4\n");

  /* We have to allocate some memory here */
  PARAMS *par = new PARAMS;
  printf("optimizeAMs: 5a\n");

  par->ans = &ans;
  printf("optimizeAMs: 5b\n");
  par->ProblemSpace = &ProblemSpace;
  printf("optimizeAMs: 5c\n");

  minex_func.n = 2;
  minex_func.f = my_f;
  minex_func.params = par;
  printf("optimizeAMs: 6\n");

  s = gsl_multimin_fminimizer_alloc (T,2);
  gsl_multimin_fminimizer_set (s, &minex_func, AMs, ss);
  printf("optimizeAMs: 7\n");

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if(status){ break;}

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      if(status == GSL_SUCCESS){
	printf("converged to minimum at\n");}
      printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
	      iter,
	      gsl_vector_get (s->x, 0),
	      gsl_vector_get (s->x, 1),
	      s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_vector_free(AMs);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  delete par;

  return;
}
