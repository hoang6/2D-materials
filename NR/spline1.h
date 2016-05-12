#ifndef _NR_SPLINE1_H
#define _NR_SPLINE1_H

#include "util.h"
#include <algorithm>

namespace NR {

/********** prototypes **********/
template<class T>
void spline(T x[], T y[], int n, T yp1, T ypn, T y2[]);

template<class T>
void splint(T xa[], T ya[], T y2a[], int n, T x, T *y);

class Error_splint {};

template<class T>
void splint_max(T xa[], T ya[], T y2a[], int n, T *xmax, T *ymax);

template<class T>
void splint_min(T xa[], T ya[], T y2a[], int n, T *xmin, T *ymin);

/********** definitions **********/
template<class T>
void spline(T x[], T y[], int n, T yp1, T ypn, T y2[]) {
  int i,k;
  T p,qn,sig,un,*u;

  Vector<T> vu(1,n-1);
  u=vu.pointer();
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
}

template<class T>
void splint(T xa[], T ya[], T y2a[], int n, T x, T *y) {
  int klo,khi,k;
  T h,b,a;
  
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) {
    throw Error_splint();
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

template<class T>
void splint_max(T xa[], T ya[], T y2a[], int n, T *xmax, T *ymax) {
  int j, k, max_i;
  double A, B, C, D; //Ax^2+Bx+C = 0
  T xroot[2];
  T* xcand;
  T* ycand;
  int count;

  Vector<T> vxcand(1,n+(n-1)*2);
  Vector<T> vycand(1,n+(n-1)*2);
  xcand=vxcand.pointer();
  ycand=vycand.pointer();
  for(k = 1; k <= n; k++) {
    xcand[k] = xa[k];
    ycand[k] = ya[k];
  }
  count = n;

  for(k = 1; k <= n-1; k++) {
    A = 3.0*(y2a[k+1]-y2a[k]);
    B = 2.0*(-3.0*xa[k]*y2a[k+1]+3*xa[k+1]*y2a[k]);
    C = 
      2*xa[k+1]*xa[k]*(y2a[k+1]-y2a[k])+
      xa[k]*xa[k]*(2.0*y2a[k+1]+y2a[k])-
      xa[k+1]*xa[k+1]*(y2a[k+1]+2*y2a[k])+
      6.0*(ya[k+1]-ya[k]);
    D = B*B-4*A*C;
    if(D < 0.0) continue;
    D = sqrt(D);
    xroot[0] = T((-B+D)/2.0/A);
    xroot[1] = T((-B-D)/2.0/A);
    for(j = 0; j < 2; j++) {
      if(xroot[j] > xa[k] && xroot[j] < xa[k+1]) {
	count++;
	xcand[count] = xroot[j];
	splint(xa,ya,y2a,n,xcand[count],&ycand[count]);
      }
    }//end for(j...
  }//end for(k...

  max_i = std::max_element(ycand+1,ycand+count+1)-ycand;
  *xmax = xcand[max_i];
  *ymax = ycand[max_i];
}

template<class T>
void splint_line_max(T xa[], T ya[], T y2a[], int n, 
		      T slope, T b, T *xmax, T *ymax) {
  int j, k, max_i;
  double A, B, C, D; //Ax^2+Bx+C = 0
  T xroot[2];
  T* xcand;
  T* ycand;
  int count;

  Vector<T> vxcand(1,n+(n-1)*2);
  Vector<T> vycand(1,n+(n-1)*2);
  xcand=vxcand.pointer();
  ycand=vycand.pointer();
  for(k = 1; k <= n; k++) {
    xcand[k] = xa[k];
    ycand[k] = ya[k];
    ycand[k] += slope*xcand[k]+b;
  }
  count = n;

  for(k = 1; k <= n-1; k++) {
    A = 3.0*(y2a[k+1]-y2a[k]);
    B = 2.0*(-3.0*xa[k]*y2a[k+1]+3*xa[k+1]*y2a[k]);
    C = 
      2*xa[k+1]*xa[k]*(y2a[k+1]-y2a[k])+
      xa[k]*xa[k]*(2.0*y2a[k+1]+y2a[k])-
      xa[k+1]*xa[k+1]*(y2a[k+1]+2*y2a[k])+
      6.0*(ya[k+1]-ya[k])+slope*(xa[k+1]-xa[k])*6.0;
    D = B*B-4*A*C;
    if(D < 0.0) continue;
    D = sqrt(D);
    xroot[0] = T((-B+D)/2.0/A);
    xroot[1] = T((-B-D)/2.0/A);
    for(j = 0; j < 2; j++) {
      if(xroot[j] > xa[k] && xroot[j] < xa[k+1]) {
	count++;
	xcand[count] = xroot[j];
	splint(xa,ya,y2a,n,xcand[count],&ycand[count]);
	ycand[count] += slope*xcand[count]+b;
      }
    }//end for(j...
  }//end for(k...

  max_i = std::max_element(ycand+1,ycand+count+1)-ycand;
  *xmax = xcand[max_i];
  *ymax = ycand[max_i];
}

} /* NR */

#endif
