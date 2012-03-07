#ifndef _KSHORTESTH
#define _KSHORTESTH

#include "DataStructures.h"
#include "pathUtils.h"
#include "shortestPath.h"
#include <map>
#include <vector>

/* Main K-shortest algorithm */
void kShortest(vector<vector<PATH> > &result, RXNSPACE &rxnspace, METSPACE &metspace, METSPACE &inputs, 
	       METSPACE &outputs, int K);
void kShortest(vector<PATH> &result, RXNSPACE &rxnspace, METSPACE &metspace, METSPACE &inputs, 
	       METABOLITE output, int K);
void kShortest2(vector<vector<PATH> > &result, RXNSPACE &rxnspace, METSPACE &metspace, METSPACE &inputs, 
	        METSPACE &outputs, int K, RXNSPACE &truedir);
void kShortest2(vector<PATH> &result, RXNSPACE &rxnspace, METSPACE &metspace, METSPACE &inputs, 
	        METABOLITE output, int K, RXNSPACE &truedir);
#endif
