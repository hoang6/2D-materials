#ifndef _MC_ALLOC_H
#define _MC_ALLOC_H

#include "MCBasic.h"
#include "MC_REBO.h"
#include "MC_LCBOPII.h"
#include "MCAgent.h"
#include "MCAgent_CtrlZ.h"
#include "MCAgent_STube.h"
#include <string>
#include <iostream>

namespace MC {

MCBasic* allocMC(const std::string& potTag);

MCAgent* allocMCAgent(const std::string& agentTag, std::istream& fin);

}/* MC */

#endif
