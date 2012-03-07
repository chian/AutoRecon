#ifndef _VISUAL
#define _VISUAL

#include <cstdio>
#include <vector>
#include "DataStructures.h"
#include "shortestPath.h"

void Paths2Dot(FILE *dotput, const METSPACE &metspace, const RXNSPACE &rxnspace, const char *label);
void convertPathForVisuals(const PATH &path, const PROBLEM &ProblemSpace, RXNSPACE &rxnspace, METSPACE &metspace, int useSyns);
void VisualizePath2File(const char *fileName, const char *label, const PATH &path, 
			const PROBLEM &ProblemSpace, int useSyns);
void visualizePathSummary2File(const char* fileBase, const char* label, const vector<PATHSUMMARY> &psum, const PROBLEM &ProblemSpace, int useSyns);

void VisualizeProblem(PROBLEM &small_model, const int i);
void VisualizeAllSolutions(const PROBLEM &ProblemSpace, const vector< vector< PATH> > &paths, int useSyns);
bool hasValidSecondaryPair(const REACTION &rxn, const METSPACE &workingMets, int met_id);
#endif
