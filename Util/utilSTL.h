#ifndef _UTIL_UTIL_STL_H
#define _UTIL_UTIL_STL_H

#include <algorithm>

namespace Util {

template<class CONT> void purge(CONT& cont);

template<class CONT, class T> bool push_back_unique(CONT& cont, const T& t);

/************************************************************/
template<class CONT>
void purge(CONT& cont) {
  typename CONT::iterator it;
  for(it = cont.begin(); it != cont.end(); it++) {
    delete (*it);
    *it = NULL;
  }
}

template<class CONT, class T>
bool push_back_unique(CONT& cont, const T& t) {
  if(std::find(cont.begin(),cont.end(),t) != cont.end()) return false;
  cont.push_back(t);
  return true;
}

}/* Util */

#endif
