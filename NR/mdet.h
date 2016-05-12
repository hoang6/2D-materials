#ifndef _NR_MDET_H
#define _NR_MDET_H

#include "ludcmp.h"

namespace NR {

/************************************************************/
//notice: the input matrix a will be destroyed
template<class T>
T mdet(T **a, int n);

/************************************************************/
template<class T>
T mdet(T **a, int n) {
  T d;
  int j, *indx;
  
  Vector<int> vindx(1,n);
  indx=vindx.pointer();
  ludcmp(a,n,indx,&d);
  for(j=1;j<=n;j++) d *= a[j][j];

  return d;
}

}/* NR */

#endif
