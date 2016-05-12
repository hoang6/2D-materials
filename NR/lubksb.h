#ifndef _NR_LUBKSB_H
#define _NR_LUBKSB_H

/************************************************************/
template<class T>
void lubksb(T **a, int n, int *indx, T b[]);

/************************************************************/
template<class T>
void lubksb(T **a, int n, int *indx, T b[]) {
  int i,ii=0,ip,j;
  T sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

#endif
