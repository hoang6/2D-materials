#ifndef _UTIL_RNG_H
#define _UTIL_RNG_H

#include "../RNG/SFMT.h"
#include <cstdlib>

namespace Util {

namespace RNG {

void seed(unsigned int _seed);

double uniform_SFMT(); //[0, 1)

double uniform_std();  //[0, 1)

/************************************************************/
inline void seed(unsigned int _seed) { 
  init_gen_rand(_seed); 
  srand(_seed); 
}

inline double uniform_SFMT() { 
  return genrand_res53(); 
}

inline double uniform_std() { 
  return double(std::rand())/(double(RAND_MAX)+1.0); 
}

}/* RNG */

}/* Util */

#endif
