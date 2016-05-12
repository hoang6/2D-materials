#ifndef _NR_JACOBI_H
#define _NR_JACOBI_H

#include "util.h"

namespace NR {

/********** prototypes **********/
template<class T>
void jacobi(T **a, int n, T d[], T **v, int *nrot);

class Error_jacobi {};

/********** definitions **********/
#define NR_JACOBI_ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

template<class T>
void jacobi(T **a, int n, T d[], T **v, int *nrot) {
  //input matrix a[1...n][1...n]
  //eigen values d[1...n]
  //columns of v[1...n][1...n] contain normalized eigenvectors
  //nrot number of rotation needed
  int j,iq,ip,i;
  T tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  Vector<T> vb(1,n);
  Vector<T> vz(1,n);
  b=vb.pointer();
  z=vz.pointer();
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++)
	sm += ABS(a[ip][iq]);
    }
    if (sm == 0.0) {
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
	g=100.0*ABS(a[ip][iq]);
	if (i > 4 && (T)(ABS(d[ip])+g) == (T)ABS(d[ip])
	    && (T)(ABS(d[iq])+g) == (T)ABS(d[iq]))
	  a[ip][iq]=0.0;
	else if (ABS(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((T)(ABS(h)+g) == (T)ABS(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(ABS(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    NR_JACOBI_ROTATE(a,j,ip,j,iq)
	  }
	  for (j=ip+1;j<=iq-1;j++) {
	    NR_JACOBI_ROTATE(a,ip,j,j,iq)
	  }
	  for (j=iq+1;j<=n;j++) {
	    NR_JACOBI_ROTATE(a,ip,j,iq,j)
	  }
	  for (j=1;j<=n;j++) {
	    NR_JACOBI_ROTATE(v,j,ip,j,iq)
	  }
	  ++(*nrot);
	}
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  throw Error_jacobi();
}

#undef NR_JACOBI_ROTATE

}

#endif
