#include <cstdio>
#include <cstdlib>

#ifndef MYCONST_H
#define MYCONST_H

/* See MyConstants.cc for definitions and values for all of these switches.
   If you want to change a value you must re-compile for it to take effect */
class DEBUGFLAGS{
 public:

  bool TEST_SYNTHESIS;

  int GAPFILL_K;
  int INITIAL_K;
  double ANNOTE_CUTOFF_1;
  double ANNOTE_CUTOFF_2;

  bool DEBUGETC;
  bool DEBUGSYN;
  bool DEBUGPATHS;
  bool DEBUGRUNK;
  bool DEBUGGAPFILL;
  bool DEBUGBRIDGES;
  bool DEBUGFVA;
  bool DEBUGFBA;
  bool PARSIMONY;
  bool OUTPUTPATHRESULTS;
  bool VISUALIZEPATHS;
  bool PRINTGAPFILLRESULTS;
  bool PRINTETCRESULTS;
  bool PRINTANNOTATIONS;
  bool PRINTSHOULDGROW;
  bool SAVEANNOTATIONS;
  int SYNFACTOR;
  int REVFACTOR;
  int MAGICBRIDGEFACTOR;
  int BLACKMAGICFACTOR;
  int MISSINGEXCHANGEFACTOR;
  int MISSINGTRANSPORTFACTOR;
  int BRIDGEFIXFACTOR;
  int ETCFACTOR;
  int MINFACTORSPACING;
  int BRIDGEMETFACTOR;
  int BIOMASS;
  /* Database conventions */
  char E_tag[8];
  char H_name[8];
  char Na_name[8];
  char Na_plus_E[16];
  char H_plus_E[16];
  char ATPM_name[16];
  char ATP_name[8];
  char ADP_name[8];
  char H2O_name[8];
  char PI_name[8];

  double FLUX_CUTOFF;
  double GROWTH_CUTOFF;

  DEBUGFLAGS();
};

/* Global debugflags variable - should be used in all files in the project instead of each one defning its own copy  */
const DEBUGFLAGS _db;

#endif
