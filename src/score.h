#ifndef _SCORE_H_
#define _SCORE_H_

#include<map>
#include<vector> 
#include"DataStructures.h"

using std::vector;
using std::map;

ANSWER overwriteRXNS(vector<REACTION> &ETC_adjusted, const ANSWER &ans);
map<int,double> getSecretionRates(const ANSWER &ans1, const GROWTH &growth1,
				  const vector<double> &fba_solution);
double evalProtocal1(const vector<double> sim, const vector<double> exp);
double measureScore(SCORE1 &score1, vector<ANSWER> &ans, vector<GROWTH> &growth);

void optimizeAMs(vector<ANSWER> &ans, PROBLEM &ProblemSpace, double init_AM);

#endif
