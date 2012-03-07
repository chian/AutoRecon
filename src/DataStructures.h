// Data structures used throughout the project

#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

#include <algorithm>
#include <climits>
#include <map>
#include <string>
#include <vector>

using std::vector;
using std::map;
using std::string;

/* Input data classes */
struct MEDIA;
struct KNOCKOUT;
struct ANNOTATION;
class STOICH;
class REACTION;
class GROWTH;
class METABOLITE;
class RXNSPACE;
class METSPACE;

class PATH;
class PATHSUMMARY;

class VALUESTORE;
class GRAPHSTORE;
class BADIDSTORE;

class GAPFILLRESULT;
class ANSWER;
class SIMULATIONRESULT;
class KORESULT;
class SCORE;

class RXNSPACE{
 public:
  vector<REACTION> rxns;

  RXNSPACE();
  RXNSPACE(const vector<REACTION> &rxnVec);
  RXNSPACE(const RXNSPACE& existingSpace, const vector<int> &idSubset);

  void clear();

  void changeReversibility(int id, int new_rev);
  void change_Lb(int id, double new_lb);
  void change_Ub(int id, double new_ub);
  void change_Lb_and_Ub(int id, double new_lb, double new_ub);

  void changeId(int oldId, int newId);

  void addReactionVector(const vector<REACTION> &rxnVec);
  void addReaction(const REACTION &rxn);
  void removeRxnFromBack();
  REACTION rxnFromId(int id) const;
  REACTION* rxnPtrFromId(int id);
  const REACTION* rxnPtrFromId(int id) const;
  int idxFromId(int id) const;
  bool idIn(int id) const;
  void rxnMap();

  RXNSPACE operator=(const RXNSPACE& init);
  REACTION & operator[](int idx);
  bool operator==(const RXNSPACE &rhs);

 private:
  /* Note - this is just a REFERENCE POINT - it always starts at 0 and all changes to ATPM are relative to whatever the user inputs */
  double currentAtpm;
  map<int,int> Ids2Idx;
  int numRxns;

};

class METSPACE{
 public:
  vector<METABOLITE> mets;

  METSPACE();
  METSPACE(const vector<METABOLITE> &metVec);
  METSPACE(const METSPACE &existingSpace, const vector<int> &idSubset);
  METSPACE(const RXNSPACE &rxnspace, const METSPACE &largeMetSpace);

  void clear();
  void removeMetFromBack();
  void addMetabolite(const METABOLITE &met);
  METABOLITE metFromId(int id) const;
  METABOLITE* metPtrFromId(int id);
  int idxFromId(int id) const;
  bool idIn(int id) const;
  void metMap();

  METSPACE operator=(const METSPACE& init);
  METABOLITE & operator[](int idx);
  map<int, int> Ids2Idx;

 private:
  int numMets;
};

class PROBLEM{
 public:
  RXNSPACE fullrxns;
  RXNSPACE synrxns;
  RXNSPACE synrxnsR;
  RXNSPACE exchanges;

  METSPACE secondaryLones;
  METSPACE metabolites;
  METSPACE cofactors;

  vector<GROWTH> growth;

  PROBLEM(const RXNSPACE &x, const PROBLEM &ProblemSpace);
  PROBLEM();
  void clear();
};

class GROWTH{
 public:
  /* Note that these should be ID's not Indexes */
  vector<MEDIA> media;
  vector<MEDIA> byproduct;
  /* FIXME: Needs secretion rates here (default would be 1000) */
  vector<STOICH> biomass;

  vector<KNOCKOUT> mutation;
  double growth_rate;

  GROWTH();
  void reset();
};

struct MEDIA{
  int id; /* Metabolite id */
  char name[64]; /* Metabolie name */
  double rate; /* mmol/gDW/hr 
		  (uptake rate for media, secretion rate for byproducts */
  bool operator==(const MEDIA &rhs) const;
  bool operator>(const MEDIA &rhs) const;
  bool operator<(const MEDIA &rhs) const;
  
  MEDIA();
};

struct KNOCKOUT{
  int id;
  string genename;
  double act_coef; /* percentage of maximum growth rate achieved with knockout (0 = lethal knockout, 1 = wild type growth rate) */
};

class METABOLITE{
 public:
  /* Externally (XML/User) defined parameters */
  int id; /* Has to be big matrix row index */
  char name[64];
  int charge;
  int input; /* Is it an input?: 0 for no, 1 for yes */
  int output; /* Is it an output? 0 for no, 1 for yes */
  int biomass; /* Is it a biomass component? 0 for no, 1 for yes */
  int secondary_lone; /* 0 for no, 1 for yes */
  vector<int> secondary_pair; /* Int of the ID for each possible secondary pair */
  int noncentral; /* 0 = central (not used for ETC), 1 = noncentral (used for ETC), -1 = undefined */
  char chemform[64];
  double modifier; /* Reserved for things that can be used to modify metabolite cost in Dijkstras algorithm, such as metabolomics data or other thigns we calculate */
  
  /* Connected reactions */
  vector<int> rxnsInvolved_nosec;
  
  METABOLITE();
  void reset();
  int isSecondary() const;
};

class STOICH{
  public:
  int met_id;
  double rxn_coeff;
  char met_name[64];
  
  bool operator==(const STOICH &rhs) const;
  bool operator>(const STOICH &rhs) const;
  bool operator<(const STOICH &rhs) const;
  STOICH();
  void reset();
};

class REACTION{
 public:
  
  int id;
  int synthesis; /* ID for metabolite that the REACTION synthesizes - if any (-1 otherwise) */
  char name[64];
  vector<STOICH> stoich; /* Full chemical reaction */

  /* stoich, but without the secondary metabolties (Load_Stoic_Part deals with this) */
  /* This is the one used by Dijkstras - if you want to use fullrxns with dijkstras you need to copy it over first */
  vector<STOICH> stoich_part;

  int transporter; /* 0 for no, 1 for yes */

  /* Reactions with either no reactants or no products [including secondaries] */
  int isExchange;

  /* Special flag for reactions that go from secondaries ONLY to a metabolite
     that is not a secondary.   ex. R00005 2 CO2 + 2 NH3 <=> 	Urea-1-carboxylate + H2O */
  int freeMakeFlag;

  double init_likelihood; /* based on user input:
			     1>x>0 PROBABILITY INCLUSION
			     -1 DO NOT INCLUDE (100% sure)
			     -2 HARD INCLUDE (100% sure)
			     -3 BLACK MAGIC
			     -4 SPONTANEOUS REACTION
			     -5 NOT ON HIT LIST
			  */

  double current_likelihood;   /* Likelihood used by Dijkstras - after modifying them, put them here */
  double old_likelihood; /* Temporary storage used by K-shortest to save old likelihood when setting new one to -1 */
  int init_reversible; /* -1 backwards only, 0 reversible, 1 forward only */
  int net_reversible;   /* Reversibility used by Dijkstras and Simplex */
  double lb;  double ub; /* LB and UB (should be consistent with net_reversible - net_reversible < 0 means lb < 0 and ub = 0) */

  int revPair; /* The iD of the reaction gonig the other way in the network (if there is one) and -1 if none */
  vector<int> syn;    /* Special vector for storing reactions in the same reaction-class */
  vector<ANNOTATION> annote; /* List of gene annotations */
  
  bool operator==(const REACTION &rhs) const;
  bool operator>(const REACTION &rhs) const;
  bool operator<(const REACTION &rhs) const;
  REACTION();
  void reset();
};

struct ANNOTATION{
  double probability;
  string genename;
  /* THis compares > because we want to go in opposite order.. when sorting these. */
  bool operator<(const ANNOTATION &rhs) const;
};

class NETREACTION{
 public:
  /* Reaction IDs that create the NETREACTION (negative numbers are used in reverse) */
  vector<int> rxnDirIds;
  /* The NETREACTION itself */
  REACTION rxn;
  bool operator==(const NETREACTION &rhs) const;
};

class PATH{
 public:
  int outputId; /* The Path by definition is trying to reach a specific output. This lists that 
		   output */
  vector<int> inputIds; /* Inputs that were needed to reach the output */
  vector<int> rxnIds; /* Reactions in the path */
  vector<int> rxnDirection; /* Dijkstras will return this... which direction is the rxn going in 
			       the particular path? */
  vector<int> metsConsumedIds; /* Metabolites that are nodes on the path tree for this particular 
				  path */
  vector<int> deadEndIds; /* Dead end metabolites on the way to reach the specified output 
			     (EMPTY if no path is found) */
  vector<int> rxnPriority; /* Priority for particular reaction [higher priority means the reaction was closer to the beginning of the pathway] - useful for gapfind */
  double totalLikelihood; /* Sum of provided likelihoods for a given path */
  PATH();
};

class PATHSUMMARY{
 public:
  long int id;
  int outputId; /* Output (target component) for the given path - needed to run FBA / magic 
		   exits */
  vector<int> growthIdx; /* Which growth does it correspond to? - MB - did we decide on an Id or Idx for 
		    this? I guess since we'll number the growth ourselves and probably just 
		    enumerate it will be the same thing */
  int k_number; /* I assume this is whether this is the 1st, 2nd or kth path in a list */
  vector<int> rxnDirIds;
  vector<int> rxnPriority;
  vector<int> deadEndIds;
  vector<int> metsConsumed;
  double likelihood;
  PATHSUMMARY();
  bool metIn(int metId) const;
  bool operator==(const PATHSUMMARY &rhs1);
  PATHSUMMARY clear();
  PATHSUMMARY operator=(PATH &onepath);
  PATHSUMMARY operator/(PATH &onepath);
};

/* This is to be used only for the purposes of the priority queue. I put it this way so that
   we can keep track of index but still use the priority queue to make minimization more efficient. */
class VALUESTORE{
 public:
  int id;
  double value;
  /* The > is not a typo - we are looking for the MINIMUM and not the MAXIMUM so we need to switch the sign
     to make it work with priority_queue... 
   Incidentally this also makes it so that a vector<VALUESTORE> sorts max to min instead of min to max ... */
  bool operator<(const VALUESTORE &rhs) const { return this[0].value > rhs.value; }
};

/* Setup priority queue for K-shortest (again, we want the minimum and not the maximum)
   We exclude specific reactions from a particular PATH and then pass those onto the next iteration */
class GRAPHSTORE{
 public:
  vector<int> excludedRxnIds;
  PATH path;  
  bool operator<(const GRAPHSTORE &rhs) const { return this[0].path.totalLikelihood > rhs.path.totalLikelihood; }
};

class GAPFILLRESULT{
 public:  
  int deadMetId; /* ID of any essential magic exits given the specified combination of PATHSUMMARY */
  vector< vector<int> > deadEndSolutions; /* Dijkstras solutions for metabolite deadMetId - deadEndSolutions[i] is the i'th shortest gapfill solution */
};

/* Setup for storage of ID's of reactions that cannot be used in Dijkstras algorithm and which metablites block them.
   Since we use a SET to store them we need a "<" operator - I order them first by reaction ID and then by length of the met ID vectors,
   and finally by the values in the met ID vector themselves (if they are the same)

   To remove the values I also need an == operator. */
class BADIDSTORE {
 public:
  int badRxnId;
  vector<int> badMetIds;
  BADIDSTORE();
  bool operator<(const BADIDSTORE &rhs) const;
  bool operator==(const BADIDSTORE &rhs) const;
};

class INNERPOPSTORE {
 public:
  double score;
  vector<int> whichK;
  vector<int> essentialExits;
  INNERPOPSTORE();
  bool operator<(const INNERPOPSTORE &rhs) const;
  bool operator==(const INNERPOPSTORE &rhs) const;
};

class SCORE1{
 public:
  double score;

  //parameter
  double gam;
  double ngam;
  vector<REACTION> ETC_adjusted;
  map<int,double> growthRate;
  map<int,double> growthScore;
};

class ANSWER {
 public:
  /*  Reaction IDs in the model
      Path numbers
      Bottlenecks fixed
       Parameters (ETC/ATPM)
       Growth rates (Experimental)
       Growth rates (Predicted by algorithm)
       Secretion rates (Experimental) for each media in the PROBLEM
       Secretion rates (Predicted by algorithm) for each media in the PROBLEM
       Total network likelihood */

  static map<string, vector<VALUESTORE> > annoteToRxns;
  static vector<NETREACTION> etc;
  static vector<vector<vector<PATHSUMMARY> > > pList;

  RXNSPACE reactions;
  METSPACE metabolites;
  /* Set of path IDs used in model (the genetic algorithm will shuffle these) */
  vector<const PATHSUMMARY*> passedPaths;

  vector<const METABOLITE*> fixedGaps;
  vector<vector<const REACTION*> > fixes; //indexed by the fixedGaps they solve 

  vector<const REACTION*> essentialMagicExits;

  /* ETC IDs (not sure how I'll do this) for ETCs connected to our network via ETC_CONNECT */
  vector<const REACTION*> etcIds; 
  vector<const REACTION*> etcConnect;
  vector<KORESULT> knockoutResults;

  //optimum scoring answer
  SCORE1 optLocalScore;

  ANSWER();

  //copy constructor
  ANSWER(const ANSWER &source){
    /* These aren't necessary because these variables are declared as STATIC so they have the same value automatically
       as any other ANSWER */
    /*    this->annoteToRxns = source.annoteToRxns;
    this->etc = source.etc;
    this->pList = source.pList; */

    this->reactions = source.reactions;
    this->metabolites = source.metabolites;
    this->knockoutResults = source.knockoutResults;
    this->optLocalScore = source.optLocalScore;

    /* This should be OK since passedPaths is a STATIC it should always be in the same memory location... */
    this->passedPaths = source.passedPaths;

    /* Since the reactions and metabolites are NOT static (by design)
       we need to move the pointers to the new answer */
    for(int i=0; i<source.fixedGaps.size(); i++) { this->fixedGaps.push_back( this->metabolites.metPtrFromId( source.fixedGaps[i]->id )); }
    for(int i=0; i<source.fixes.size(); i++) {
      vector<const REACTION*> tmp;
      for(int j=0; j<source.fixes[i].size(); j++) {
	tmp.push_back(this->reactions.rxnPtrFromId( source.fixes[i][j]->id ));
      }
      this->fixes.push_back(tmp);
    }
    for(int i=0; i<source.essentialMagicExits.size(); i++) { this->essentialMagicExits.push_back( this->reactions.rxnPtrFromId( source.essentialMagicExits[i]->id )); }
    for(int i=0; i<source.etcIds.size(); i++) { this->etcIds.push_back( this->reactions.rxnPtrFromId( source.etcIds[i]->id )); }
    for(int i=0; i<source.etcConnect.size(); i++) { this->etcConnect.push_back( this->reactions.rxnPtrFromId( source.etcConnect[i]->id )); }
  }

};

class KORESULT{

};

/* Initialize static members of the ANSWER class */
void initializeAnswer(const vector<NETREACTION> &myEtc, const vector<vector<vector<PATHSUMMARY> > > &myPlist, const map<string, vector<VALUESTORE> > &annoteToRxn);

#endif // _DATASTRUCTURES_H
