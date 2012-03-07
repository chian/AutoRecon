#ifndef ANNOTATIONS_H_
#define ANNOTATIONS_H_

#include <map>
#include <vector>
#include "pathUtils.h"
#include "DataStructures.h"

using std::vector;
using std::map;
using std::string;

map< string, vector<VALUESTORE>  > GeneAnnotations(RXNSPACE &small_list, double cut, double cut2);
void removeAnnotation(RXNSPACE &rxnspace, const string &genename, const VALUESTORE &val);

#endif
