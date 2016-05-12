/*----------------------------radic.c---------------------------------
  This file contains routines for calculating Tij and Fij as defined
  in DW Brenner J. Phys Cond Mat 14 (2002) 783-802.  This is
  accomplished through the use of tricubic interpolation.  This
  software was written by Stephen Jordan and modeled on code from
  Zyvex as well as DW Brenner's original Fortran.
  ---------------------------------------------------------------------*/

#include "radic.h"
#include "in3.h"
#include "globalConstants.h"
#include <cmath>
#include <cstdio>
#include <iostream>

/* The first coordinate of in3 is an index into CLMN or TLMN.
   The second coordinate is either 1, 2, or 3 for L (XINT1), M (XINT2), or N
   (CONJUG). The value od in3 is the power to raise XINT1, XINT2, or CONJUG
   to before multiplying by the corresponding element of CLMN or TLMN and
   then adding to evaluate the cubic spline.  CLMN contains the coefficients
   for the polynomial interpolation in RADIC.  Similarly, TLMN contains the
   coefficients for polynomial interpolation in TOR.  In3_initialized keeps
   track of whether the coefficients, etc. been loaded yet.*/

namespace tbtools {

//store computed Fij and Tij
#define TYPE_MAX 3
#define N_MAX 11
//value
static bool FijFlag[TYPE_MAX][TYPE_MAX][N_MAX][N_MAX][N_MAX];
static bool TijFlag[N_MAX][N_MAX][N_MAX];
static double FijVal[TYPE_MAX][TYPE_MAX][N_MAX][N_MAX][N_MAX];
static double TijVal[N_MAX][N_MAX][N_MAX];
//derivative
static bool dFijdXNT1Flag[TYPE_MAX][TYPE_MAX][N_MAX][N_MAX][N_MAX];
static bool dFijdXNT2Flag[TYPE_MAX][TYPE_MAX][N_MAX][N_MAX][N_MAX];
static bool dFijdCONJUGFlag[TYPE_MAX][TYPE_MAX][N_MAX][N_MAX][N_MAX];
static bool dTijdXNT1Flag[N_MAX][N_MAX][N_MAX];
static bool dTijdXNT2Flag[N_MAX][N_MAX][N_MAX];
static bool dTijdCONJUGFlag[N_MAX][N_MAX][N_MAX];
static double dFijdXNT1Val[TYPE_MAX][TYPE_MAX][N_MAX][N_MAX][N_MAX];
static double dFijdXNT2Val[TYPE_MAX][TYPE_MAX][N_MAX][N_MAX][N_MAX];
static double dFijdCONJUGVal[TYPE_MAX][TYPE_MAX][N_MAX][N_MAX][N_MAX];
static double dTijdXNT1Val[N_MAX][N_MAX][N_MAX];
static double dTijdXNT2Val[N_MAX][N_MAX][N_MAX];
static double dTijdCONJUGVal[N_MAX][N_MAX][N_MAX];

/* Fij(NiT,NjT,Nijconj) */
/* See equation 14 in the paper */
double Fij(int KI, int KJ, double XNT1, double XNT2, double CONJUG) {
  int J,L,M,N,KIKJ;
  double X, RAD;
  const double *clmn_ptr;
  double xnt1_pow[4];
  double xnt2_pow[4];
  double conj_pow[4];

  //for fast computation
  int iXNT1=int(XNT1);
  int iXNT2=int(XNT2);
  int iCONJUG=int(CONJUG);
  bool iflag=0;
  if(KI<TYPE_MAX && KJ<TYPE_MAX &&
     fabs(iXNT1-XNT1)<EPS_DOUBLE && iXNT1<N_MAX && 
     fabs(iXNT2-XNT2)<EPS_DOUBLE && iXNT2<N_MAX &&
     fabs(iCONJUG-CONJUG)<EPS_DOUBLE && iCONJUG<N_MAX) iflag = 1;
  if(iflag && FijFlag[KI][KJ][iXNT1][iXNT2][iCONJUG]) {
    return FijVal[KI][KJ][iXNT1][iXNT2][iCONJUG];
  }

  //  printf("Fij: %i\t%i\t%.20e\t%.20e\t%.20e\t", KI, KJ, XNT1, XNT2, CONJUG);
  XNT1++; //SPJ: to make it consistent with the paper
  XNT2++; //SPJ: likewise
  if(!in3_initialized) init_in3();
  L=(int)floor(XNT1);
  M=(int)floor(XNT2);
  N=(int)floor(CONJUG);
  RAD = 0;
  KIKJ = KI + KJ - 1;
  if(L >= 4) {
    L = 4;
    XNT1 = 4;
  }
  if(M >= 4) {
    M = 4;
    XNT2 = 4;
  }
  if(M < 1) {
    M = 1;
    XNT2 = 1;
  }
  if(L < 1) {
    L = 1;
    XNT1 = 1;
  }
  if(N >= 9) {
    N = 9;
    CONJUG = 9;
  }
  xnt1_pow[0] = 1.0;
  xnt1_pow[1] = XNT1;
  xnt1_pow[2] = XNT1*XNT1;
  xnt1_pow[3] = XNT1*xnt1_pow[2];
  xnt2_pow[0] = 1.0;
  xnt2_pow[1] = XNT2;
  xnt2_pow[2] = XNT2*XNT2;
  xnt2_pow[3] = XNT2*xnt2_pow[2];
  conj_pow[0] = 1.0;
  conj_pow[1] = CONJUG; 
  conj_pow[2] = CONJUG*CONJUG;
  conj_pow[3] = CONJUG*conj_pow[2];
  clmn_ptr = &CLMN[KIKJ][L][M][N][0];
  for (J=1; J<65; ++J) {
    X=*++clmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
    RAD += X;
  }
  //for symmetry
  clmn_ptr = &CLMN[KIKJ][M][L][N][0];
  for (J=1; J<65; ++J) {
    X=*++clmn_ptr*xnt2_pow[IN3[J][1]]*xnt1_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
    RAD += X;
  }
  RAD = RAD/2.0;
  //  printf("%.20e\n", RAD);
  if(iflag) {
    FijFlag[KI][KJ][iXNT1][iXNT2][iCONJUG] = 1;
    FijVal[KI][KJ][iXNT1][iXNT2][iCONJUG] = RAD;
  }
  return RAD;
}

/*Tij(NiT, NjT, Nijconj))*/
/* Tricubic spline for the torsional rotation.  It's Tij as in equation
   18 of the paper. */
double Tij(double XNT1, double XNT2, double CONJUG) {
  int J,L,M,N;
  double X;
  double ATOR;

  //for fast computation
  int iXNT1=int(XNT1);
  int iXNT2=int(XNT2);
  int iCONJUG=int(CONJUG);
  bool iflag=0;
  if(fabs(iXNT1-XNT1)<EPS_DOUBLE && iXNT1<N_MAX && 
     fabs(iXNT2-XNT2)<EPS_DOUBLE && iXNT2<N_MAX &&
     fabs(iCONJUG-CONJUG)<EPS_DOUBLE && iCONJUG<N_MAX) iflag = 1;
  if(iflag && TijFlag[iXNT1][iXNT2][iCONJUG]) 
    return TijVal[iXNT1][iXNT2][iCONJUG];

  //  printf("Tij:%.20e\t%.20e\t%.20e\t", XNT1, XNT2, CONJUG);
  XNT1++;  //SPJ: to make it consistent with the paper
  XNT2++;  //SPJ: likewise
  if(!in3_initialized) init_in3();
  ATOR = 0;
  if(XNT1<4 || XNT2<4) {
    const double *tlmn_ptr;
    double xnt1_pow[4];
    double xnt2_pow[4];
    double conj_pow[4];
    xnt1_pow[0] = 1.0;
    xnt1_pow[1] = XNT1;
    xnt1_pow[2] = XNT1*XNT1;
    xnt1_pow[3] = XNT1*xnt1_pow[2];
    xnt2_pow[0] = 1.0;
    xnt2_pow[1] = XNT2;
    xnt2_pow[2] = XNT2*XNT2;
    xnt2_pow[3] = XNT2*xnt2_pow[2];
    conj_pow[0] = 1.0;
    conj_pow[1] = CONJUG;
    conj_pow[2] = CONJUG*CONJUG;
    conj_pow[3] = CONJUG*conj_pow[2];
    L = (int)floor(XNT1);
    M = (int)floor(XNT2);
    N = (int)floor(CONJUG);
    tlmn_ptr = &TLMN[L][M][N][0];
    for(J=1; J<65; J++) {
      X = *++tlmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
      ATOR += X;
    }
    //for symmetry
    tlmn_ptr = &TLMN[M][L][N][0];
    for (J=1; J<65; ++J) {
      X=*++tlmn_ptr*xnt2_pow[IN3[J][1]]*xnt1_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
      ATOR += X;
    }
    ATOR = ATOR/2.0;
  }
  //  printf("%.20e\n", ATOR);
  if(iflag) {
    TijFlag[iXNT1][iXNT2][iCONJUG] = 1;
    TijVal[iXNT1][iXNT2][iCONJUG] = ATOR;
  }
  return ATOR;
}

double dFij_dXNT1(int KI, int KJ, double XNT1, double XNT2, double CONJUG) {
  int J,L,M,N,KIKJ;
  double X, RAD;
  const double *clmn_ptr;
  double xnt1_pow[4];
  double xnt2_pow[4];
  double conj_pow[4];
  
  //for fast computation
  int iXNT1=int(XNT1);
  int iXNT2=int(XNT2);
  int iCONJUG=int(CONJUG);
  bool iflag=0;
  if(KI<TYPE_MAX && KJ<TYPE_MAX &&
     fabs(iXNT1-XNT1)<EPS_DOUBLE && iXNT1<N_MAX && 
     fabs(iXNT2-XNT2)<EPS_DOUBLE && iXNT2<N_MAX &&
     fabs(iCONJUG-CONJUG)<EPS_DOUBLE && iCONJUG<N_MAX) iflag = 1;
  if(iflag && dFijdXNT1Flag[KI][KJ][iXNT1][iXNT2][iCONJUG]) {
    return dFijdXNT1Val[KI][KJ][iXNT1][iXNT2][iCONJUG];
  }

  XNT1++; //SPJ: to make it consistent with the paper
  XNT2++; //SPJ: likewise
  if(!in3_initialized) init_in3();
  L=(int)floor(XNT1);
  M=(int)floor(XNT2);
  N=(int)floor(CONJUG);
  RAD = 0;
  KIKJ = KI + KJ - 1;
  if(L >= 4) return 0.0;
  if(M >= 4) {
    M = 4;
    XNT2 = 4;
  }
  if(M < 1) {
    M = 1;
    XNT2 = 1;
  }
  if(L < 1) return 0.0;
  if(N >= 9) {
    N = 9;
    CONJUG = 9;
  }
  xnt1_pow[0] = 0.0;
  xnt1_pow[1] = 1.0;
  xnt1_pow[2] = 2*XNT1;
  xnt1_pow[3] = 3*XNT1*XNT1;
  xnt2_pow[0] = 1.0;
  xnt2_pow[1] = XNT2;
  xnt2_pow[2] = XNT2*XNT2;
  xnt2_pow[3] = XNT2*xnt2_pow[2];
  conj_pow[0] = 1.0;
  conj_pow[1] = CONJUG; 
  conj_pow[2] = CONJUG*CONJUG;
  conj_pow[3] = CONJUG*conj_pow[2];
  clmn_ptr = &CLMN[KIKJ][L][M][N][0];
  for (J=1; J<65; ++J) {
    X=*++clmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
    RAD += X;
  }
  //for symmetry
  clmn_ptr = &CLMN[KIKJ][M][L][N][0];
  for (J=1; J<65; ++J) {
    X=*++clmn_ptr*xnt2_pow[IN3[J][1]]*xnt1_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
    RAD += X;
  }
  RAD = RAD/2.0;
  if(iflag) {
    dFijdXNT1Flag[KI][KJ][iXNT1][iXNT2][iCONJUG] = 1;
    dFijdXNT1Val[KI][KJ][iXNT1][iXNT2][iCONJUG] = RAD;
  }
  return RAD;
}
  
double dFij_dXNT2(int KI, int KJ, double XNT1, double XNT2, double CONJUG) {
  int J,L,M,N,KIKJ;
  double X, RAD;
  const double *clmn_ptr;
  double xnt1_pow[4];
  double xnt2_pow[4];
  double conj_pow[4];

  //for fast computation
  int iXNT1=int(XNT1);
  int iXNT2=int(XNT2);
  int iCONJUG=int(CONJUG);
  bool iflag=0;
  if(KI<TYPE_MAX && KJ<TYPE_MAX &&
     fabs(iXNT1-XNT1)<EPS_DOUBLE && iXNT1<N_MAX && 
     fabs(iXNT2-XNT2)<EPS_DOUBLE && iXNT2<N_MAX &&
     fabs(iCONJUG-CONJUG)<EPS_DOUBLE && iCONJUG<N_MAX) iflag = 1;
  if(iflag && dFijdXNT2Flag[KI][KJ][iXNT1][iXNT2][iCONJUG]) {
    return dFijdXNT2Val[KI][KJ][iXNT1][iXNT2][iCONJUG];
  }
  
  XNT1++; //SPJ: to make it consistent with the paper
  XNT2++; //SPJ: likewise
  if(!in3_initialized) init_in3();
  L=(int)floor(XNT1);
  M=(int)floor(XNT2);
  N=(int)floor(CONJUG);
  RAD = 0;
  KIKJ = KI + KJ - 1;
  if(L >= 4) {
    L = 4;
    XNT1 = 4;
  }
  if(M >= 4) return 0.0;
  if(M < 1) return 0.0;
  if(L < 1) {
    L = 1;
    XNT1 = 1;
  }
  if(N >= 9) {
    N = 9;
    CONJUG = 9;
  }
  xnt1_pow[0] = 1.0;
  xnt1_pow[1] = XNT1;
  xnt1_pow[2] = XNT1*XNT1;
  xnt1_pow[3] = XNT1*xnt1_pow[2];
  xnt2_pow[0] = 0.0;
  xnt2_pow[1] = 1.0;
  xnt2_pow[2] = 2*XNT2;
  xnt2_pow[3] = 3*XNT2*XNT2;
  conj_pow[0] = 1.0;
  conj_pow[1] = CONJUG; 
  conj_pow[2] = CONJUG*CONJUG;
  conj_pow[3] = CONJUG*conj_pow[2];
  clmn_ptr = &CLMN[KIKJ][L][M][N][0];
  for (J=1; J<65; ++J) {
    X=*++clmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
    RAD += X;
  }
  //for symmetry
  clmn_ptr = &CLMN[KIKJ][M][L][N][0];
  for (J=1; J<65; ++J) {
    X=*++clmn_ptr*xnt2_pow[IN3[J][1]]*xnt1_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
    RAD += X;
  }
  RAD = RAD/2.0;
  if(iflag) {
    dFijdXNT2Flag[KI][KJ][iXNT1][iXNT2][iCONJUG] = 1;
    dFijdXNT2Val[KI][KJ][iXNT1][iXNT2][iCONJUG] = RAD;
  }
  return RAD;
}
  
double dFij_dCONJUG(int KI, int KJ, double XNT1, double XNT2, double CONJUG) {
  int J,L,M,N,KIKJ;
  double X, RAD;
  const double *clmn_ptr;
  double xnt1_pow[4];
  double xnt2_pow[4];
  double conj_pow[4];

  //for fast computation
  int iXNT1=int(XNT1);
  int iXNT2=int(XNT2);
  int iCONJUG=int(CONJUG);
  bool iflag=0;
  if(KI<TYPE_MAX && KJ<TYPE_MAX &&
     fabs(iXNT1-XNT1)<EPS_DOUBLE && iXNT1<N_MAX && 
     fabs(iXNT2-XNT2)<EPS_DOUBLE && iXNT2<N_MAX &&
     fabs(iCONJUG-CONJUG)<EPS_DOUBLE && iCONJUG<N_MAX) iflag = 1;
  if(iflag && dFijdCONJUGFlag[KI][KJ][iXNT1][iXNT2][iCONJUG]) {
    return dFijdCONJUGVal[KI][KJ][iXNT1][iXNT2][iCONJUG];
  }

  XNT1++; //SPJ: to make it consistent with the paper
  XNT2++; //SPJ: likewise
  if(!in3_initialized) init_in3();
  L=(int)floor(XNT1);
  M=(int)floor(XNT2);
  N=(int)floor(CONJUG);
  RAD = 0;
  KIKJ = KI + KJ - 1;
  if(L >= 4) {
    L = 4;
    XNT1 = 4;
  }
  if(M >= 4) {
    M = 4;
    XNT2 = 4;
  }
  if(M < 1) {
    M = 1;
    XNT2 = 1;
  }
  if(L < 1) {
    L = 1;
    XNT1 = 1;
  }
  if(N >= 9) return 0.0;
  xnt1_pow[0] = 1.0;
  xnt1_pow[1] = XNT1;
  xnt1_pow[2] = XNT1*XNT1;
  xnt1_pow[3] = XNT1*xnt1_pow[2];
  xnt2_pow[0] = 1.0;
  xnt2_pow[1] = XNT2;
  xnt2_pow[2] = XNT2*XNT2;
  xnt2_pow[3] = XNT2*xnt2_pow[2];
  conj_pow[0] = 0.0;
  conj_pow[1] = 1.0;
  conj_pow[2] = 2*CONJUG;
  conj_pow[3] = 3*CONJUG*CONJUG;
  clmn_ptr = &CLMN[KIKJ][L][M][N][0];
  for (J=1; J<65; ++J) {
    X=*++clmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
    RAD += X;
  }
  //for symmetry
  clmn_ptr = &CLMN[KIKJ][M][L][N][0];
  for (J=1; J<65; ++J) {
    X=*++clmn_ptr*xnt2_pow[IN3[J][1]]*xnt1_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
    RAD += X;
  }
  RAD = RAD/2.0;
  if(iflag) {
    dFijdCONJUGFlag[KI][KJ][iXNT1][iXNT2][iCONJUG] = 1;
    dFijdCONJUGVal[KI][KJ][iXNT1][iXNT2][iCONJUG] = RAD;
  }
  return RAD;
}

double dTij_dXNT1(double XNT1, double XNT2, double CONJUG) {
  int J,L,M,N;
  double X;
  double ATOR;
  
  //for fast computation
  int iXNT1=int(XNT1);
  int iXNT2=int(XNT2);
  int iCONJUG=int(CONJUG);
  bool iflag=0;
  if(fabs(iXNT1-XNT1)<EPS_DOUBLE && iXNT1<N_MAX && 
     fabs(iXNT2-XNT2)<EPS_DOUBLE && iXNT2<N_MAX &&
     fabs(iCONJUG-CONJUG)<EPS_DOUBLE && iCONJUG<N_MAX) iflag = 1;
  if(iflag && dTijdXNT1Flag[iXNT1][iXNT2][iCONJUG]) 
    return dTijdXNT1Val[iXNT1][iXNT2][iCONJUG];

  XNT1++;  //SPJ: to make it consistent with the paper
  XNT2++;  //SPJ: likewise
  if(!in3_initialized) init_in3();
  ATOR = 0;
  if(XNT1<4 || XNT2<4) {
    const double *tlmn_ptr;
    double xnt1_pow[4];
    double xnt2_pow[4];
    double conj_pow[4];
    xnt1_pow[0] = 0.0;
    xnt1_pow[1] = 1.0;
    xnt1_pow[2] = 2*XNT1;
    xnt1_pow[3] = 3*XNT1*XNT1;
    xnt2_pow[0] = 1.0;
    xnt2_pow[1] = XNT2;
    xnt2_pow[2] = XNT2*XNT2;
    xnt2_pow[3] = XNT2*xnt2_pow[2];
    conj_pow[0] = 1.0;
    conj_pow[1] = CONJUG;
    conj_pow[2] = CONJUG*CONJUG;
    conj_pow[3] = CONJUG*conj_pow[2];
    L = (int)floor(XNT1);
    M = (int)floor(XNT2);
    N = (int)floor(CONJUG);
    tlmn_ptr = &TLMN[L][M][N][0];
    for(J=1; J<65; J++) {
      X = *++tlmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
      ATOR += X;
    }
    //for symmetry
    tlmn_ptr = &TLMN[M][L][N][0];
    for (J=1; J<65; ++J) {
      X=*++tlmn_ptr*xnt2_pow[IN3[J][1]]*xnt1_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
      ATOR += X;
    }
    ATOR = ATOR/2.0;
  }
  if(iflag) {
    dTijdXNT1Flag[iXNT1][iXNT2][iCONJUG] = 1;
    dTijdXNT1Val[iXNT1][iXNT2][iCONJUG] = ATOR;
  }
  return ATOR;
}

double dTij_dXNT2(double XNT1, double XNT2, double CONJUG) {
  int J,L,M,N;
  double X;
  double ATOR;

  //for fast computation
  int iXNT1=int(XNT1);
  int iXNT2=int(XNT2);
  int iCONJUG=int(CONJUG);
  bool iflag=0;
  if(fabs(iXNT1-XNT1)<EPS_DOUBLE && iXNT1<N_MAX && 
     fabs(iXNT2-XNT2)<EPS_DOUBLE && iXNT2<N_MAX &&
     fabs(iCONJUG-CONJUG)<EPS_DOUBLE && iCONJUG<N_MAX) iflag = 1;
  if(iflag && dTijdXNT2Flag[iXNT1][iXNT2][iCONJUG]) 
    return dTijdXNT2Val[iXNT1][iXNT2][iCONJUG];

  XNT1++;  //SPJ: to make it consistent with the paper
  XNT2++;  //SPJ: likewise
  if(!in3_initialized) init_in3();
  ATOR = 0;
  if(XNT1<4 || XNT2<4) {
    const double *tlmn_ptr;
    double xnt1_pow[4];
    double xnt2_pow[4];
    double conj_pow[4];
    xnt1_pow[0] = 1.0;
    xnt1_pow[1] = XNT1;
    xnt1_pow[2] = XNT1*XNT1;
    xnt1_pow[3] = XNT1*xnt1_pow[2];
    xnt2_pow[0] = 0.0;
    xnt2_pow[1] = 1;
    xnt2_pow[2] = 2*XNT2;
    xnt2_pow[3] = 3*XNT2*XNT2;
    conj_pow[0] = 1.0;
    conj_pow[1] = CONJUG;
    conj_pow[2] = CONJUG*CONJUG;
    conj_pow[3] = CONJUG*conj_pow[2];
    L = (int)floor(XNT1);
    M = (int)floor(XNT2);
    N = (int)floor(CONJUG);
    tlmn_ptr = &TLMN[L][M][N][0];
    for(J=1; J<65; J++) {
      X = *++tlmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
      ATOR += X;
    }
    //for symmetry
    tlmn_ptr = &TLMN[M][L][N][0];
    for (J=1; J<65; ++J) {
      X=*++tlmn_ptr*xnt2_pow[IN3[J][1]]*xnt1_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
      ATOR += X;
    }
    ATOR = ATOR/2.0;
  }
  if(iflag) {
    dTijdXNT2Flag[iXNT1][iXNT2][iCONJUG] = 1;
    dTijdXNT2Val[iXNT1][iXNT2][iCONJUG] = ATOR;
  }
  return ATOR;
}

double dTij_dCONJUG(double XNT1, double XNT2, double CONJUG) {
  int J,L,M,N;
  double X;
  double ATOR;

  //for fast computation
  int iXNT1=int(XNT1);
  int iXNT2=int(XNT2);
  int iCONJUG=int(CONJUG);
  bool iflag=0;
  if(fabs(iXNT1-XNT1)<EPS_DOUBLE && iXNT1<N_MAX && 
     fabs(iXNT2-XNT2)<EPS_DOUBLE && iXNT2<N_MAX &&
     fabs(iCONJUG-CONJUG)<EPS_DOUBLE && iCONJUG<N_MAX) iflag = 1;
  if(iflag && dTijdCONJUGFlag[iXNT1][iXNT2][iCONJUG]) 
    return dTijdCONJUGVal[iXNT1][iXNT2][iCONJUG];

  XNT1++;  //SPJ: to make it consistent with the paper
  XNT2++;  //SPJ: likewise
  if(!in3_initialized) init_in3();
  ATOR = 0;
  if(XNT1<4 || XNT2<4) {
    const double *tlmn_ptr;
    double xnt1_pow[4];
    double xnt2_pow[4];
    double conj_pow[4];
    xnt1_pow[0] = 1.0;
    xnt1_pow[1] = XNT1;
    xnt1_pow[2] = XNT1*XNT1;
    xnt1_pow[3] = XNT1*xnt1_pow[2];
    xnt2_pow[0] = 1.0;
    xnt2_pow[1] = XNT2;
    xnt2_pow[2] = XNT2*XNT2;
    xnt2_pow[3] = XNT2*xnt2_pow[2];
    conj_pow[0] = 0.0;
    conj_pow[1] = 1.0;
    conj_pow[2] = 2*CONJUG;
    conj_pow[3] = 3*CONJUG*CONJUG;
    L = (int)floor(XNT1);
    M = (int)floor(XNT2);
    N = (int)floor(CONJUG);
    tlmn_ptr = &TLMN[L][M][N][0];
    for(J=1; J<65; J++) {
      X = *++tlmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
      ATOR += X;
    }
    //for symmetry
    tlmn_ptr = &TLMN[M][L][N][0];
    for (J=1; J<65; ++J) {
      X=*++tlmn_ptr*xnt2_pow[IN3[J][1]]*xnt1_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
      ATOR += X;
    }
    ATOR = ATOR/2.0;
  }
  if(iflag) {
    dTijdCONJUGFlag[iXNT1][iXNT2][iCONJUG] = 1;
    dTijdCONJUGVal[iXNT1][iXNT2][iCONJUG] = ATOR;
  }
  return ATOR;  
}

#undef TYPE_MAX
#undef N_MAX

}/* tbtools */
