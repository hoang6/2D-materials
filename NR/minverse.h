#ifndef _NR_MINVERSE_H
#define _NR_MINVERSE_H

#include "ludcmp.h"
#include "lubksb.h"

namespace NR {

/************************************************************/
//notice: the input matrix a will be destroyed
template<class T>
void minverse(T **a, T **result, int n);

/************************************************************/
template<class T>
void minverse(T **a, T **result, int n) {
  T d,*col;
  int i,j,*indx;

  Vector<T> vcol(1,n);
  Vector<int> vindx(1,n);
  col=vcol.pointer();
  indx=vindx.pointer();
  ludcmp(a,n,indx,&d);
  for(j=1;j<=n;j++) {
    for(i=1;i<=n;i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(a,n,indx,col);
    for(i=1;i<=n;i++) result[i][j]=col[i];
  }
}

}/* NR */

#endif
