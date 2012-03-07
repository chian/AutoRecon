#ifndef COMPONENTS_H
#define COMPONENTS_H

#include "shortestPath.h"
#include "pathUtils.h"
#include "visual01.h"
#include "MersenneTwister.h"
#include "DataStructures.h"
#include "RunK.h"
#include "Exchanges.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>

/*Component by component FBA testing*/
int ComponentTest(const PROBLEM &theModel);
PROBLEM reportFlux(PROBLEM &model, vector<double> &fbaSolution, double cutoff);
PROBLEM reportFlux(PROBLEM &model, vector<double> &fbaSolution);

#endif
