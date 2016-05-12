#ifndef _UTIL_CONSTANTS_H
#define _UTIL_CONSTANTS_H

#include <limits>

namespace Util {

/*************** math constants ***************/  
const double PI = 3.14159265358979323846;
const double SIN60 = 0.86602540378443864676;
const double MAX_DOUBLE = std::numeric_limits<double>::max();
const double MIN_DOUBLE = std::numeric_limits<double>::min();
const double EPS_DOUBLE = std::numeric_limits<double>::epsilon();
  
/*************** physical constants ***************/  
const double BOLTZMANN_CONSTANT = 8.61734315e-5;               /* eV/K */
const double PLANCK_CONSTANT = 4.1356673310e-15;               /* eV*s */
const double REDUCED_PLANCK_CONSTANT = PLANCK_CONSTANT/2.0/PI; /* eV*s */
const double AVOGADRO_CONSTANT = 6.0221417930e23;              /* mol^(-1) */
const double ELECTRON_CHARGE = 1.60217648740e-19;              /* C */

}/* Util */

#endif
