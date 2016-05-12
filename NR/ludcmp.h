#ifndef _NR_LUDCMP_H
#define _NR_LUDCMP_H

#include "util.h"

namespace NR {

/************************************************************/
template<class T>
void ludcmp(T **a, int n, int *indx, T *d);

class Error_ludcmp {};

/************************************************************/
template<class T>
void ludcmp(T **a, int n, int *indx, T *d) {
  const T TINY(0.0);
  int i,imax,j,k;
  T big,dum,sum,temp;
  T *vv;
  
  Vector<T> vvv(1,n);
  vv = vvv.pointer();
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) {
      throw Error_ludcmp();
    }
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    imax = 0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
}

}/* NR */

#endif
