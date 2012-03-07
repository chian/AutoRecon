/* Contains functions that I deleted but don't know if I will ever need them again */

void modifyRxnsForMagicBridge(PROBLEM &problemSpace, int bridgeId, vector<int> &rxnIds, vector<bool> &toErase) {
  REACTION bridge = problemSpace.synrxns.rxnFromId(bridgeId);
  RXNSPACE &fullrxns = problemSpace.fullrxns;

  int origSize = fullrxns.rxns.size(); /* Because the size of fullrxns will change */
  for(int j=0; j<origSize; j++) {
    for(int i=0; i<bridge.stoich_part.size(); i++) {
      /* Change 11-8-11 - add the biomass equation to the list of things that need to be given _MB */
      if( isDanglingSecondary(fullrxns.rxns[j], problemSpace.metabolites, bridge.stoich_part[i].met_id) ||
          fullrxns.rxns[j].id == _db.BIOMASS) {
        if(_db.DEBUGBRIDGES) {printf("Changing metabolite %s to _MB version in reaction %s\n", bridge.stoich_part[i].met_name, fullrxns.rxns[j].name); }
        changeMetToMB(fullrxns.rxns[j], bridge.stoich_part[i].met_id, bridge.stoich_part[i].met_name, true);
      } else {
        if(fixedBridge(fullrxns.rxns[j], bridge)) {
          /* CHANGE 10-5-11: We need to actually have two copies of these reactions, one with _MB (to connect to the path) and one without
             (to allow cofactor generation). Some of these reactions, e.g. FTHFL, are used as part of cofactor regeneration loops and also
             as part of the synthesis pathways of other metabolites, so we need both. However, this could help us, because it effectively separates
             the cofactor internal usage from its synthesis. */
          REACTION tmpRxn = fullrxns.rxns[j];
          tmpRxn.magicBridgeFix.clear(); /* Prevents us from getting _MB_MB_MB_MB_... multiple copies of rxns */
          tmpRxn.id += _db.BRIDGEFIXFACTOR;
          sprintf(tmpRxn.name, "%s%s", tmpRxn.name, "_MB");

          changeMetToMB(tmpRxn, bridge.stoich_part[0].met_id, bridge.stoich_part[0].met_name, false);
          changeMetToMB(tmpRxn, bridge.stoich_part[1].met_id, bridge.stoich_part[1].met_name, false);
          fullrxns.addReaction(tmpRxn);
          rxnIds.push_back(tmpRxn.id);
          toErase.push_back(false);

          break; /* Prevents modifying the same reaction twice (once for each of the two magic bridge metabolites) */
        }
      }
    }
  }
}

/* Returns true if the reaction "rxn" was used to fix the bridge "bridge" and false otherwise */
bool fixedBridge(const REACTION &rxn, const REACTION &bridge) {
  for(int i=0; i<rxn.magicBridgeFix.size(); i++) {
    if(rxn.magicBridgeFix[i] == bridge.id) { return true; }
  }
  return false;
}

/* Change a metabolite in the REACTION to the MB version [magic bridge version], being careful to only remove ONE of them
 if the stoichiometry is not 1:1 (e.g. GLUSx has 2 glu-L on the right side and we only one ONE to be converted to _MB and one to remain as glu_L

If it's a dangling secondary (isDangling = true), then we want to replace all of them regardless of stoichiometry...
[this is to prevent the biomass equation from becoming wonky]

If the met_id is not present in rxn this function does nothing to the reaction. */
void changeMetToMB(REACTION &rxn, int met_id, const char* metName, bool isDangling) {
  int origSize = rxn.stoich.size();
  for(int i=0; i<origSize; i++) {
    if(rxn.stoich[i].met_id == met_id) {
      if(isDangling || rougheq(rxn.stoich[i].rxn_coeff, 1, 1E-5) == 1 || rougheq(rxn.stoich[i].rxn_coeff, -1, 1E-5) == 1 ) {
        // Just replace the name and ID, and don't bother adding a new component to the STOICH.
        rxn.stoich[i].met_id = met_id + _db.BRIDGEMETFACTOR;
        sprintf(rxn.stoich[i].met_name, "%s%s", metName, "_MB");
      } else {
        STOICH tmpStoich;
        // We always want one less to be reacted in total, regardless of which side of the reaction it's on (I think)
        if(rxn.stoich[i].rxn_coeff > 0) {
          rxn.stoich[i].rxn_coeff -= 1.0f;
          tmpStoich.rxn_coeff = 1.0f;
        } else {
          rxn.stoich[i].rxn_coeff += 1.0f;
          tmpStoich.rxn_coeff = -1.0f;
        }
        tmpStoich.met_id = met_id + _db.BRIDGEMETFACTOR;
        sprintf(tmpStoich.met_name, "%s%s", metName, "_MB");
        rxn.stoich.push_back(tmpStoich);
        break;
      }
    }
  }
  return;
 }

/* Returns TRUE if metabolite with ID "met_id" is a dangling secondary in reaction "reaction".

Dangling secondary:  if metabolite  in REACTION that has secondary pairs, but NONE of them
are present in the reaction.

Use this function to identify if we should modify the STOICH to contain the "magic bridge" version
of the secondary metabolite instead of the "normal" version

This reaction is intended to be used on fullrxns (so it uses STOICH, not STOICH_PART) */
bool isDanglingSecondary(const REACTION &reaction, const METSPACE &metspace, int met_id) {
  bool hasDesiredMet(false);

  METABOLITE tmpMet = metspace.metFromId(met_id);
  if(tmpMet.secondary_pair.empty()) {
    printf("WARNING: Met ID passed to isDanglingSecondary that has no secondary metabolites. This is most likely an error.");
    assert(false);
  }


  /* If the met_id isn't in the REACTION, we definitely shouldn't try to change the stoich for met_id */
  for(int i=0; i<reaction.stoich.size(); i++) {
    if(reaction.stoich[i].met_id == met_id) { hasDesiredMet = true; break; }
  }
  if(!hasDesiredMet) { return false; }

  /* If the met_id is in there, look for any of the secondary_pairs. If any is in there return FALSE (not a dangling secondary) otherwise return TRUE */
  bool hasSecondaryPair(false);
  for(int i=0; i<reaction.stoich.size(); i++) {
    for(int j=0; j<tmpMet.secondary_pair.size(); j++) {
      if(tmpMet.secondary_pair[j] == reaction.stoich[i].met_id) { return false; }
    }
  }

  return true;

}


/* adds "Bridge Metabolites" to the end of the METSPACE.                                                                                                                                                       
   The "Bridge Metabolites" are meant to be place-holders for                                                                                                                                                  
   metabolites that are normally part of a cofactor pair, but are                                                                                                                                              
   utilized as part of the normal path (e.g. 10fthf and thf in the path for                                                                                                                                      synthesis of 5mthf in E. coli)                                                                                                                                                                                                                                                                                                                                                                                            The placeholder is necessary to avoid conflicts with reactions in which they ARE a secondary pair                                                                                                           
   after doing the gapfilling step... */
void addBridgeMetabolites(PROBLEM &problemSpace) {
  DEBUGFLAGS db;
  METSPACE &metspace = problemSpace.metabolites;

  int initSize = metspace.mets.size(); // Since the size changes we have to declare this before the loop.                                                                                                      
  for(int i=0; i<initSize; i++) {
    if(metspace.mets[i].secondary_pair.empty()) { continue; }
    METABOLITE tmpMet = metspace.mets[i];
    sprintf(tmpMet.name, "%s%s", metspace.mets[i].name, "_MB");
    tmpMet.id = tmpMet.id + db.BRIDGEMETFACTOR;
    /* Should NOT be considered a cofactor or cofactor pair (so we don't make magic entrances for these) */
    tmpMet.secondary_lone = 0;
    tmpMet.secondary_pair.clear();

    metspace.addMetabolite(tmpMet);
  }

  return;
}
