#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "MyConstants.h"

DEBUGFLAGS::DEBUGFLAGS() {

  /* Test new methods for synthesis reactions- if FALSE we ignore the synthesis tag */
  TEST_SYNTHESIS = true;

  /******************* K-shortest parameters *********/
  INITIAL_K = 1;
  GAPFILL_K = 3;

  /* Annotation cutoff - remove annotations that are less than this
     probability compared to the maximum  [one reaction --> multiple genes cutoff]
   Originally 0.7 */
  ANNOTE_CUTOFF_1 = 0.7;

  /* Annotation cutoff 2 - remove annotations that are less than this 
     probability compared to the maximum [one gene --> multiple reactions cutoff] */
  ANNOTE_CUTOFF_2 = 0.2;

  /******************** Algorithmic switches *********/
  /* True if you want to use maximum-parsimony instaed of maximum-likelihood */
  PARSIMONY = false;

  /************ Output file switches ******************/
  /* True if you want to output visualization of paths */
  VISUALIZEPATHS = true;
  /* True if you want to output some tab-delimited text files with path information */
  OUTPUTPATHRESULTS = true;
  /* True if you want to print "SHOULDGROW", which is a text file containing all the reactions after the LP but before
     gapfilling (and which seems to be a perrenial pain in my ass to get it to work correctly) */
  PRINTSHOULDGROW = true;
  /* True to print a tab-delimited list of reaction annotations to a file */
  SAVEANNOTATIONS = true;

  /************ Printer switches **********************/
  /* True if you want to print debugging information for the electron transport chain finder / linker */
  DEBUGETC = false;
  /* True if you want a printout of the reaction synonyms */
  DEBUGSYN = false;
  /* True if you want a printout of all of the paths found after running K shortest [suggested true unless you get bothered by the ETC chain dead end filling printouts] */
  DEBUGPATHS = true;
  /* True if you want a list of reactions passed into K-shortest every time runK is run (very large, suggested false) */
  DEBUGRUNK = false;
  /* True if you want to print out information about gapfind / gapfill (more detailed than PRINTGAPFILLRESULTS - suggested false) */
  DEBUGGAPFILL = false;
  /* True if you want to print detailed information about the magic bridge fixing */
  DEBUGBRIDGES = false;
  /* True to get printout of FVA results */
  DEBUGFVA = false;
  /* True to get printouts of FBA results (including GLPK printouts) */
  DEBUGFBA = false;
  /* True if you want to get gapfill results (suggested true) */
  PRINTGAPFILLRESULTS = true;
  /* True to get a list of ETCs */
  PRINTETCRESULTS = false;
  /* True to print annotations for each reaction */
  PRINTANNOTATIONS = false;

  /************ FLUX and GROWTH RATE precision ( DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING!) ************/
  FLUX_CUTOFF = 1E-7;
  GROWTH_CUTOFF = 1E-5;

  /************ Internal constants (DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING!) ********/
  MISSINGTRANSPORTFACTOR = 2000000;/* ID for transports added because none were present for media components */
  MISSINGEXCHANGEFACTOR = 3000000; /* For reactions added because no exchanges were present for media components */
  MAGICBRIDGEFACTOR = 4000000;     /* ID for magic bridges */
  BLACKMAGICFACTOR = 5000000;      /* ID for magic entrances and exits (used for gapfilling) */
  ETCFACTOR = 6000000;             /* Not sure what this was for */
  REVFACTOR = 7000000;             /* ID for reactions that are reversed versions of existing reactions in the network */
  SYNFACTOR = 8000000;             /* ID for the synonymized version of reactions */
  BRIDGEFIXFACTOR = 9000000;       /* ID for reactions that were used to fix magic bridges ("MB" versions) */
  BIOMASS = 99999;
  MINFACTORSPACING = 500000; /* This is the minimum distance between two special flags -
                                It's needed in order to extract later which reactions have which tag.
				It must be larger than num_mets + BRIDGEMETFACTOR and also larger than num_metabolites but must be less than the difference of any of the two factors above. */

  BRIDGEMETFACTOR = 100000; /* For bridge (_MB) metabolites [used after gapfilling] */

  /* Database conventions for important metabolites and reactions, and for compartment names (no compartment = no _x_ tag) */
  strcpy(E_tag, "_e_");
  strcpy(H_name, "h");
  strcpy(Na_name, "na1");
  sprintf(H_plus_E, "%s%s", H_name, E_tag);
  sprintf(Na_plus_E,"%s%s", Na_name, E_tag);

  /* Non-growth associated ATP maintenance reaction */
  strcpy(ATPM_name, "ATPM");
  /* Growth-associated maintenance - we need these in order to modify the ATPM value of the biomass equation - along with H_name above */
  strcpy(ATP_name, "atp");
  strcpy(ADP_name, "adp");
  strcpy(H2O_name, "h2o");
  strcpy(PI_name, "pi");
}

