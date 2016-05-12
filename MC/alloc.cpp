#include "alloc.h"

namespace MC {

MCBasic* allocMC(const std::string& potTag) {
  MCBasic* pmc = NULL;
  if(potTag == "REBO") {
    pmc = new MC::MC_REBO;
  }
  else if(potTag == "LCBOPII") {
    pmc =  new MC::MC_LCBOPII;
  }
  return pmc;
}

MCAgent* allocMCAgent(const std::string& agentTag, std::istream& fin) {
  std::string keyword;
  MCAgent* pmcAgent = NULL;
  if(agentTag == "Default") {
    pmcAgent = new MC::MCAgent;
  }
  else if(agentTag == "CtrlZ") {
    MC::MCAgent_CtrlZ* p = new MC::MCAgent_CtrlZ;
    fin >> keyword >> p->zmin()
	>> keyword >> p->zmax();
    pmcAgent = static_cast<MCAgent*>(p);
  }
  else if(agentTag == "STube") {
    pmcAgent = new MC::MCAgent_STube;
  }
  return pmcAgent;
}

}/* MC */
