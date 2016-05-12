/*************************************************************************
Copyright (c) 1980-2007, Jorge Nocedal.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

This software is freely available for educational or commercial  purposes.
We expect that all publications describing work using this software  quote
at least one of the references given below:
    * J. Nocedal. Updating  Quasi-Newton  Matrices  with  Limited  Storage
      (1980), Mathematics of Computation 35, pp. 773-782.
    * D.C. Liu and J. Nocedal. On the  Limited  Memory  Method  for  Large
      Scale  Optimization  (1989),  Mathematical  Programming  B,  45,  3,
      pp. 503-528.
*************************************************************************/

#ifndef _NR_LBFGS_H
#define _NR_LBFGS_H

#include "util.h"

/*************************************************************************
        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
                          JORGE NOCEDAL

The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
of memory.

The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm  (instead of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.

This subroutine uses the FuncGrad subroutine which calculates the value of
the function F and gradient G in point X. The programmer should define the
FuncGrad subroutine by himself.  It  should  be  noted that the subroutine
doesn't need to waste time for memory allocation of array G,  because  the
memory is allocated in calling the subroutine. Setting a dimension of array
G  each  time  when  calling  a  subroutine  will excessively slow down an
algorithm.

The programmer could also redefine the LBFGSNewIteration subroutine  which
is called on each new step. The current point X, the function value F  and
the  gradient  G  are  passed  into  this  subroutine. It is reasonable to
redefine the subroutine for better debugging, for  example,  to  visualize
the solution process.

Input parameters:
    N   -   problem dimension. N>0
    M   -   number of corrections in the BFGS scheme of Hessian
            approximation update. Recommended value:  3<=M<=7. The smaller
            value causes worse convergence, the bigger will  not  cause  a
            considerably better convergence, but will cause a fall in  the
            performance. M<=N.
    X   -   initial solution approximation.
            Array whose index ranges from 1 to N.
    EpsG -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if the condition ||G|| < EpsG  is
            satisfied, where ||.|| means Euclidian norm, G - gradient, X -
            current approximation.
    EpsF -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if on iteration  number  k+1  the
            condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}    is
            satisfied.
    EpsX -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if on iteration number k+1    the
            condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts- maximum number of iterations. If MaxIts=0, the number of
            iterations is unlimited.

Output parameters:
    X   -   solution approximation. Array whose index ranges from 1 to N.
    Info-   a return code:
                    * -1 wrong parameters were specified,
                    * 0 interrupted by user,
                    * 1 relative function decreasing is less or equal to EpsF,
                    * 2 step is less or equal EpsX,
                    * 4 gradient norm is less or equal to EpsG,
                    * 5 number of iterations exceeds MaxIts.

FuncGrad routine description. User-defined.
Input parameters:
    X   -   array whose index ranges from 1 to N.
Output parameters:
    F   -   function value at X.
    G   -   function gradient.
            Array whose index ranges from 1 to N.
The memory for array G has already been allocated in the calling subroutine,
and it isn't necessary to allocate it in the FuncGrad subroutine.
*************************************************************************/
namespace NR {

/********** prototypes **********/
template<class F, class T>
void lbfgs(T x[],
	   const int& n,
	   const int& m,
	   const T& epsg,
	   const T& epsf,
	   const T& epsx,
	   const int& maxits,
	   int* return_info,
	   F& funcd);

/********** definitions **********/

template<class F, class T>
void lbfgs(T x[],
	   const int& n,
	   const int& m,
	   const T& epsg,
	   const T& epsf,
	   const T& epsx,
	   const int& maxits,
	   int* return_info,
	   F& funcd)
{
  T* w;
  T f;
  T fold;
  T tf;
  T txnorm;
  T v;
  T* xold;
  T* tx;
  T* g;
  T* diag;
  T* ta;
  bool finish;
  T gnorm;
  T stp1;
  T ftol;
  T stp;
  T ys;
  T yy;
  T sq;
  T yr;
  T beta;
  T xnorm;
  int iter;
  int nfun;
  int point;
  int ispt;
  int iypt;
  int maxfev;
  int bound;
  int npt;
  int cp;
  int i;
  int nfev;
  int inmc;
  int iycn;
  int iscn;
  T xtol;
  T gtol;
  T stpmin;
  T stpmax;
  int& info = *return_info;

  Vector<T> vw(1, n*(2*m+1)+2*m);
  Vector<T> vg(1, n);
  Vector<T> vxold(1, n);
  Vector<T> vtx(1, n);
  Vector<T> vdiag(1, n);
  w = vw.pointer();
  g = vg.pointer();
  xold = vxold.pointer();
  tx = vtx.pointer();
  diag = vdiag.pointer();

  x = x-1;
  fold = f;
  iter = 0;
  info = 0;
  if( n<=0||m<=0||m>n||epsg<0||epsf<0||epsx<0||maxits<0 )
    {
      info = -1;
      return;
    }
  nfun = 1;
  point = 0;
  npt = point*n;
  finish = false;
  for(i = 1; i <= n; i++) 
    {
      diag[i] = 1;
    }
  
  //---------- define consts ----------
  gtol = 0.01; //0.9;
  ftol = 0.0001;
  xtol = 100*5.0e-16;
  maxfev = 20;
  stpmin = 0.001;//pow(T(10), T(-20));
  stpmax = 0.01;//pow(T(10), T(20));
  //---------- end define consts ------

  f = funcd.func(&x[1]);
  funcd.dfunc(&x[1], &g[1]);
  //debug
  funcd.callBack(&x[1], &g[1], T(0));  

  ispt = n+2*m;
  iypt = ispt+n*m;
  for(i = 1; i <= n; i++)
    {
      w[ispt+i] = -g[i]*diag[i];
    }
  gnorm = sqrt(vdotv(&g[1], &g[1], n));
  if(gnorm<=epsg )
    {
      info = 4;
      return;
    }
  stp1 = 1/gnorm;
  while(true)
    {
      vcopy(&x[1], &xold[1], n);
      iter = iter+1;
      info = 0;
      bound = iter-1;
      if( iter!=1 )
	{
	  if( iter>m )
	    {
	      bound = m;
	    }
	  ys = vdotv(&w[iypt+npt+1], &w[ispt+npt+1], n);
	  yy = vdotv(&w[iypt+npt+1], &w[iypt+npt+1], n);
	  if(ys == 0.0) ys = 1.0;
	  if(yy == 0.0) yy = 1.0;
	  for(i = 1; i <= n; i++)
	    {
	      diag[i] = ys/yy;
	    }
	  cp = point;
	  if( point==0 )
	    {
	      cp = m;
	    }
	  w[n+cp] = 1/ys;
	  for(i = 1; i <= n; i++)
	    {
	      w[i] = -g[i];
	    }
	  cp = point;
	  for(i = 1; i <= bound; i++)
	    {
	      cp = cp-1;
	      if( cp==-1 )
		{
		  cp = m-1;
		}
	      sq = vdotv(&w[ispt+cp*n+1], &w[1], n);
	      inmc = n+m+cp+1;
	      iycn = iypt+cp*n;
	      w[inmc] = w[n+cp+1]*sq;
	      kdotvadd(-w[inmc], &w[iycn+1], 1.0, &w[1], &w[1], n);
	    }
	  for(i = 1; i <= n; i++)
	    {
	      w[i] = diag[i]*w[i];
	    }
	  for(i = 1; i <= bound; i++)
	    {
	      yr = vdotv(&w[iypt+cp*n+1], &w[1], n);
	      beta = w[n+cp+1]*yr;
	      inmc = n+m+cp+1;
	      beta = w[inmc]-beta;
	      iscn = ispt+cp*n;
	      kdotvadd(beta, &w[iscn+1], 1.0, &w[1], &w[1], n);
	      cp = cp+1;
	      if( cp==m )
		{
		  cp = 0;
		}
	    }
	  for(i = 1; i <= n; i++)
	    {
	      w[ispt+point*n+i] = w[i];
	    }
	}
      nfev = 0;
      stp = 1;
      if( iter==1 )
	{
	  stp = stp1;
	}
      for(i = 1; i <= n; i++)
	{
	  w[i] = g[i];
	}

      //mc search ...
      lbfgsmcsrch(n, x, f, g, w, ispt+point*n+1, stp, ftol, xtol, maxfev, info, nfev, diag, gtol, stpmin, stpmax, funcd);
      if( info!=1 )
	{
	  if( info==0 )
	    {
	      info = -1;
	      return;
	    }
	}
      nfun = nfun+nfev;
      npt = point*n;
      for(i = 1; i <= n; i++)
	{
	  w[ispt+npt+i] = stp*w[ispt+npt+i];
	  w[iypt+npt+i] = g[i]-w[i];
	}
      point = point+1;
      if( point==m )
	{
	  point = 0;
	}
      if( iter>maxits&&maxits>0 )
	{
	  info = 5;
	  return;
	}
      
      //debug
      funcd.callBack(&x[1], &g[1], T(0));

      gnorm = sqrt(vdotv(&g[1], &g[1], n));
      if( gnorm<=epsg )
	{
	  info = 4;
	  return;
	}
      tf = MAX(fabs(fold), MAX(fabs(f), 1.0));
      if( fold-f<=epsf*tf )
	{
	  info = 1;
	}
      vcopy(&xold[1], &tx[1], n);
      vsub(&tx[1], &x[1], &tx[1], n);
      xnorm = sqrt(vdotv(&x[1], &x[1], n));
      txnorm = MAX(xnorm, sqrt(vdotv(&xold[1], &xold[1], n)));
      txnorm = MAX(txnorm, 1.0);
      v = sqrt(vdotv(&tx[1], &tx[1], n));
      if( v<=epsx )
	{
	  info = 2;
	  return;
	}
      fold = f;
      vcopy(&x[1], &xold[1], n);
    }
}

template<class F, class T>
void lbfgsmcsrch(const int& n,
		 T* x,
		 T& f,
		 T* g,
		 T* s,
		 int sstart,
		 T& stp,
		 const T& ftol,
		 const T& xtol,
		 const int& maxfev,
		 int& info,
		 int& nfev,
		 T* wa,
		 const T& gtol,
		 const T& stpmin,
		 const T& stpmax,
		 F& funcd)
{
  int infoc;
  int j;
  bool brackt;
  bool stage1;
  T dg;
  T dgm;
  T dginit;
  T dgtest;
  T dgx;
  T dgxm;
  T dgy;
  T dgym;
  T finit;
  T ftest1;
  T fm;
  T fx;
  T fxm;
  T fy;
  T fym;
  T p5;
  T p66;
  T stx;
  T sty;
  T stmin;
  T stmax;
  T width;
  T width1;
  T xtrapf;
  T zero;
  T mytemp;
    
  sstart = sstart-1;
  p5 = 0.5;
  p66 = 0.66;
  xtrapf = 4.0;
  zero = 0;
  f = funcd.func(&x[1]);
  funcd.dfunc(&x[1], &g[1]);
  infoc = 1;
  info = 0;
  if( n<=0||stp<=0||ftol<0||gtol<zero||xtol<zero||stpmin<zero||stpmax<stpmin||maxfev<=0 )
    {
      return;
    }
  dginit = 0;
  for(j = 1; j <= n; j++)
    {
      dginit = dginit+g[j]*s[j+sstart];
    }
    if( dginit>=0 )
    {
      return;
    }
    brackt = false;
    stage1 = true;
    nfev = 0;
    finit = f;
    dgtest = ftol*dginit;
    width = stpmax-stpmin;
    width1 = width/p5;
    for(j = 1; j <= n; j++)
    {
      wa[j] = x[j];
    }
    stx = 0;
    fx = finit;
    dgx = dginit;
    sty = 0;
    fy = finit;
    dgy = dginit;
    while(true)
    {
      if( brackt )
        {
	  if( stx<sty )
            {
	      stmin = stx;
	      stmax = sty;
            }
	  else
            {
	      stmin = sty;
	      stmax = stx;
            }
        }
      else
        {
	  stmin = stx;
	  stmax = stp+xtrapf*(stp-stx);
        }
      if( stp>stpmax )
        {
	  stp = stpmax;
        }
      if( stp<stpmin )
        {
	  stp = stpmin;
        }
      if( brackt&&(stp<=stmin||stp>=stmax)||nfev>=maxfev-1||infoc==0||brackt&&stmax-stmin<=xtol*stmax )
        {
	  stp = stx;
        }
      for(j = 1; j <= n; j++)
        {
	  x[j] = wa[j]+stp*s[j+sstart];
        }
      f = funcd.func(&x[1]);
      funcd.dfunc(&x[1], &g[1]);
      info = 0;
      nfev = nfev+1;
      dg = 0;
      for(j = 1; j <= n; j++)
        {
	  dg = dg+g[j]*s[j+sstart];
        }
      ftest1 = finit+stp*dgtest;
      if( brackt&&(stp<=stmin||stp>=stmax)||infoc==0 )
        {
	  info = 6;
        }
      if( stp==stpmax&&f<=ftest1&&dg<=dgtest )
        {
	  info = 5;
        }
      if( stp==stpmin&&(f>ftest1||dg>=dgtest) )
        {
	  info = 4;
        }
      if( nfev>=maxfev )
	{
	  info = 3;
	}
      if( brackt&&stmax-stmin<=xtol*stmax )
        {
	  info = 2;
        }
      if( f<=ftest1&&fabs(dg)<=-gtol*dginit )
        {
	  info = 1;
        }
      if( info!=0 )
        {
	  return;
        }
      mytemp = ftol;
      if( gtol<ftol )
        {
	  mytemp = gtol;
        }
      if( stage1&&f<=ftest1&&dg>=mytemp*dginit )
        {
	  stage1 = false;
        }
        if( stage1&&f<=fx&&f>ftest1 )
        {
	  fm = f-stp*dgtest;
	  fxm = fx-stx*dgtest;
	  fym = fy-sty*dgtest;
	  dgm = dg-dgtest;
	  dgxm = dgx-dgtest;
	  dgym = dgy-dgtest;
	  //mcstep ...
	  lbfgsmcstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc, funcd);
	  fx = fxm+stx*dgtest;
	  fy = fym+sty*dgtest;
	  dgx = dgxm+dgtest;
	  dgy = dgym+dgtest;
        }
        else
        {
	  //mcstep...
	  lbfgsmcstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc, funcd);
        }
        if( brackt )
        {
	  if( fabs(sty-stx)>=p66*width1 )
            {
	      stp = stx+p5*(sty-stx);
            }
	  width1 = width;
	  width = fabs(sty-stx);
        }
    }
}

template<class T>
void lbfgsmcstep(T& stx,
		 T& fx,
		 T& dx,
		 T& sty,
		 T& fy,
		 T& dy,
		 T& stp,
		 const T& fp,
		 const T& dp,
		 bool& brackt,
		 const T& stmin,
		 const T& stmax,
		 int& info,
		 T* pT)
{
  bool bound;
  T gamma;
  T p;
  T q;
  T r;
  T s;
  T sgnd;
  T stpc;
  T stpf;
  T stpq;
  T theta;
  
  info = 0;
  if( brackt&&(stp<=MIN(stx, sty)||stp>=MAX(stx, sty))||dx*(stp-stx)>=0||stmax<stmin )
    {
      return;
    }
  sgnd = dp*(dx/fabs(dx));
  if( fp>fx )
    {
      info = 1;
      bound = true;
      theta = 3*(fx-fp)/(stp-stx)+dx+dp;
      s = MAX(fabs(theta), MAX(fabs(dx), fabs(dp)));
      gamma = s*sqrt(pow(theta/s, 2)-dx/s*(dp/s));
      if( stp<stx )
        {
	  gamma = -gamma;
        }
      p = gamma-dx+theta;
      q = gamma-dx+gamma+dp;
      r = p/q;
      stpc = stx+r*(stp-stx);
      stpq = stx+dx/((fx-fp)/(stp-stx)+dx)/2*(stp-stx);
      if( fabs(stpc-stx)<fabs(stpq-stx) )
        {
	  stpf = stpc;
        }
      else
        {
	  stpf = stpc+(stpq-stpc)/2;
        }
      brackt = true;
    }
  else
    {
      if( sgnd<0 )
        {
	  info = 2;
	  bound = false;
	  theta = 3*(fx-fp)/(stp-stx)+dx+dp;
	  s = MAX(fabs(theta), MAX(fabs(dx), fabs(dp)));
	  gamma = s*sqrt(pow(theta/s, 2)-dx/s*(dp/s));
	  if( stp>stx )
            {
	      gamma = -gamma;
            }
	  p = gamma-dp+theta;
	  q = gamma-dp+gamma+dx;
	  r = p/q;
	  stpc = stp+r*(stx-stp);
	  stpq = stp+dp/(dp-dx)*(stx-stp);
	  if( fabs(stpc-stp)>fabs(stpq-stp) )
            {
	      stpf = stpc;
            }
	  else
            {
	      stpf = stpq;
            }
	  brackt = true;
        }
      else
	{
	  if( fabs(dp)<fabs(dx) )
	    {
	      info = 3;
	      bound = true;
	      theta = 3*(fx-fp)/(stp-stx)+dx+dp;
	      s = MAX(fabs(theta), MAX(fabs(dx), fabs(dp)));
	      gamma = s*sqrt(MAX(T(0), pow(theta/s, 2)-dx/s*(dp/s)));
	      if( stp>stx )
                {
		  gamma = -gamma;
                }
	      p = gamma-dp+theta;
	      q = gamma+(dx-dp)+gamma;
	      r = p/q;
	      if( r<0&&gamma!=0 )
                {
		  stpc = stp+r*(stx-stp);
                }
	      else
                {
		  if( stp>stx )
                    {
		      stpc = stmax;
                    }
		  else
                    {
		      stpc = stmin;
                    }
                }
	      stpq = stp+dp/(dp-dx)*(stx-stp);
	      if( brackt )
                {
		  if( fabs(stp-stpc)<fabs(stp-stpq) )
                    {
		      stpf = stpc;
                    }
		  else
		    {
		      stpf = stpq;
		    }
                }
	      else
                {
		  if( fabs(stp-stpc)>fabs(stp-stpq) )
                    {
		      stpf = stpc;
                    }
		  else
                    {
		      stpf = stpq;
                    }
                }
            }
	  else
	    {
	      info = 4;
	      bound = false;
	      if( brackt )
                {
		  theta = 3*(fp-fy)/(sty-stp)+dy+dp;
		  s = MAX(fabs(theta), MAX(fabs(dy), fabs(dp)));
		  gamma = s*sqrt(pow(theta/s, 2)-dy/s*(dp/s));
		  if( stp>sty )
                    {
		      gamma = -gamma;
                    }
		  p = gamma-dp+theta;
		  q = gamma-dp+gamma+dy;
		  r = p/q;
		  stpc = stp+r*(sty-stp);
		  stpf = stpc;
                }
	      else
                {
		  if( stp>stx )
                    {
		      stpf = stmax;
                    }
		  else
                    {
		      stpf = stmin;
                    }
                }
            }
        }
    }
  if( fp>fx )
    {
      sty = stp;
      fy = fp;
      dy = dp;
    }
  else
    {
      if( sgnd<0.0 )
        {
	  sty = stx;
	  fy = fx;
	  dy = dx;
        }
      stx = stp;
      fx = fp;
      dx = dp;
    }
  stpf = MIN(stmax, stpf);
  stpf = MAX(stmin, stpf);
  stp = stpf;
  if( brackt&&bound )
    {
      if( sty>stx )
	{
	  stp = MIN(stx+0.66*(sty-stx), stp);
        }
      else
        {
	  stp = MAX(stx+0.66*(sty-stx), stp);
        }
    }
}

}/* NR */

#endif
