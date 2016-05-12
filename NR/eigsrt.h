#ifndef _NR_ERGSRT_H
#define _NR_ERGSRT_H

namespace NR {

/********** prototypes **********/
template<class T>
void eigsrt(T d[], T **v, int n);

/********** definitions **********/
template<class T>
void eigsrt(T d[], T **v, int n) {
  //sort eigenvalues into descending order
  int k,j,i;
  T p;
  for (i=1;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<=n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=1;j<=n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}

}

#endif
