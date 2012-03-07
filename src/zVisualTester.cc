#include "kShortest.h"
#include "shortestPath.h"
#include "genericLinprog.h"
#include "pathUtils.h"
#include "visual01.h"
#include "XML_loader.h"
#include "Printers.h"
#include <vector>
#include <map>

void testKShortest_visuals();

int main() {
	testKShortest_visuals();
	return 0;
}

void testKShortest_visuals() {
	int K(4);
	int maxPathLength(5000);
	
	vector<REACTION> reaction;
    vector<METABOLITE> metabolite;
	vector<GROWTH> growth;

	vector<int> inputIds;
	vector<int> outputIds;
	vector< vector<PATH> > result;
	int i, j,k;
	map<int, int> metIds2Idx, metIdx2Ids;

	char* docName = "model_likelihood_eco.xml";
	char* docName2 = "inputdata_eco.xml";
	/* IMPORTANT - need to use this instead of parseDOC and parseGROWTH if you want to use visualization */
	parseALL(docName, docName2, reaction, metabolite, growth);

	/* Get linkages between metabolites and reactions */
	ProblemSpace.metabolites.metMap();
	
	/*Deal with special flags */
	adjustLikelihoods(reaction, 1.0f, -1.0f, 1.1f, 0.0f);
	
	/* Load inputs and outputs */
	/*for(i=0; i<1;i++) {
		for(j=0;j<growth[i].media.size();j++) {
			inputIds.push_back(growth[i].media[j].id);
		}
		//For testing purpose, I only want a couple of them...
		for(j=3;j<6;j++) {
			outputIds.push_back(growth[i].biomass[j].met_id);
		}
	} */
	
	loadInputsOutputs(growth, 0, metabolite, metIds2Idx, inputIds, outputIds);
	result = kShortest(reaction, metabolite, inputIds, outputIds, maxPathLength, 1, K);
	
	vector<REACTION> newReaction;
	vector<METABOLITE> newMetabolite;
	vector<PATH> path;
	char* label = "TEST";
	char fileName[32];
	FILE *dotput;
	int status;
	
	for(i=0;i<result.size();i++) {
		path = result[i];
		printPathResults(reaction, result[i]);
		for(j=0;j<path.size();j++) {
			status = sprintf(fileName, "TESTPATH_%d_%d.dot", i,j);
			if(status == -1) { throw; }
			dotput = fopen(fileName, "w");
			convertPathForVisuals(path[j], reaction, metabolite, newReaction, newMetabolite);
			printf("METABOLITE LIST #%d:\n",j);
			for(k=0;k<newMetabolite.size();k++) {
				printMETABOLITEinputs(newMetabolite[k]);
			}
			printf("REACTION LIST #%d:\n",j);
			for(k=0;k<newReaction.size();k++) {
				printREACTIONinputs(newReaction[k], 1);
			}
			Paths2Dot(dotput, newMetabolite, newReaction, label);
			fclose(dotput);
		}
	}
}
