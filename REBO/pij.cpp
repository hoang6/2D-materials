/*---------------------------Pij.c---------------------------
  This software uses polynomial interpolation to calculate
  values of Pij, as in equation 8 of DW Brenner J Cond Mat
  14 (2002) 783-802.  It was written by Stephen Jordan and
  modeled on code by Zyvex as well as D.W. Brenner's original
  Fortran.
  -----------------------------------------------------------*/

#include "pij.h"
#include "in2.h"
#include "globalConstants.h"
#include <cmath>
#include <cstdio>
#include <iostream>

namespace tbtools {

/*IN2 stores exponents for the 16 terms of a 2D cubic polynomial.
  CLM stores the coefficients for the polynomial interpolation.
  xh stores the values of Pij where i and j are integers. */

//store computed Pij
#define TYPE_MAX 3
#define N_MAX 5
//value
static bool PijFlag[TYPE_MAX][N_MAX][N_MAX];
static double PijVal[TYPE_MAX][N_MAX][N_MAX];
//derivative
static bool dPijdNihFlag[TYPE_MAX][N_MAX][N_MAX];
static bool dPijdNicFlag[TYPE_MAX][N_MAX][N_MAX];
static double dPijdNihVal[TYPE_MAX][N_MAX][N_MAX];
static double dPijdNicVal[TYPE_MAX][N_MAX][N_MAX];


//follows normal rounding rules
double round(double input) {
  double a;
  a = floor(input);
  if((input - a) > 0.5) return a + 1.0;
  //otherwise
  return a;
}


//interpolates Pij using a bicubic spline
double spline(int KJ, double Nih, double Nic) {
  int J, NH, NC;
  double X, ANSY;
  double xx1_pow[4];
  double xx2_pow[4];
  double XX1, XX2;
  const double *clm_ptr;
  XX1 = Nih + 1.0;  //to make it match the paper
  XX2 = Nic + 1.0;  //likewise
  //bicubic interpolation or perhaps bicubic spline
  if (XX1 < 1) XX1 = 1;
  if (XX2 < 1) XX2 = 1;
  NH=(int)floor(XX1);
  NC=(int)floor(XX2);
  if(NH > MAX_BC_NEIGHBORS || NC > MAX_BC_NEIGHBORS) {
    //cerr << "Invalid neighbor numbers in Pij.\n";
    return(0.0);
  }
  if (KJ==0) {
    std::cerr << "Invalid atom type in Pij.\n";
    return(0.0);
  }
  xx1_pow[0] = 1.0;
  xx1_pow[1] = XX1;
  xx1_pow[2] = XX1*XX1;
  xx1_pow[3] = XX1*XX1*XX1;
  xx2_pow[0] = 1.0;
  xx2_pow[1] = XX2;
  xx2_pow[2] = XX2*XX2;
  xx2_pow[3] = XX2*XX2*XX2;
  clm_ptr = &CLM[KJ][NH][NC][0];
  ANSY=0;
  for (J=1; J<17; J++) {
    X = *++clm_ptr*xx1_pow[IN2[J][1]]*xx2_pow[IN2[J][2]];
    ANSY += X;
  }
  return ANSY;
}


/* Compute P sub i j, which is referenced in equation 8 on page 789.
   (of J. Phys. Cond. Mat vol 14 (2002) 783-802)  
   KJ is the identity of the "j" atom in P sub i j.
   The i atom must always be a carbon.
   If i is hydrogen don't call Pij, just use 0.
   Nih is N super H sub i and Nic is N super C sub i.
   For KJ, 1 = carbon, 2 = hydrogen.*/

double Pij(int KJ, double Nih, double Nic) {
  double rNih, rNic;
  int ih, ic;
  double return_val;

  //for fast computation
  int iNih = int(Nih);
  int iNic = int(Nic);
  bool iflag = 0;
  if(KJ < TYPE_MAX && 
     double(iNih)==Nih && iNih<N_MAX &&
     double(iNic)==Nic && iNic<N_MAX) iflag = 1;
  if(iflag && PijFlag[KJ][iNih][iNic])
    return PijVal[KJ][iNih][iNic];

  if(in2_initialized == 0) init_in2();
  rNih = round(Nih);
  rNic = round(Nic);
  ih = (int)rNih;
  ic = (int)rNic;
  if(fabs(rNih - Nih) < EPS_DOUBLE && fabs(rNic - Nic) < EPS_DOUBLE) {
    return_val = xh[KJ][ih+1][ic+1];
  }
  else
    return_val = spline(KJ, Nih, Nic);

  if(iflag) {
    PijFlag[KJ][iNih][iNic] = 1;
    PijVal[KJ][iNih][iNic] = return_val; 
  }
  return return_val;
}

double dspline_dNih(int KJ, double Nih, double Nic) {
  int J, NH, NC;
  double X, ANSY;
  double xx1_pow[4];
  double xx2_pow[4];
  double XX1, XX2;
  const double *clm_ptr;
  XX1 = Nih + 1.0;  //to make it match the paper
  XX2 = Nic + 1.0;  //likewise
  //bicubic interpolation or perhaps bicubic spline
  if (XX1 < 1) return 0.0;
  if (XX2 < 1) XX2 = 1;
  NH=(int)floor(XX1);
  NC=(int)floor(XX2);
  if(NH > MAX_BC_NEIGHBORS || NC > MAX_BC_NEIGHBORS) {
    //std::cerr << "Invalid neighbor numbers in Pij.\n";
    return(0.0);
  }
  if (KJ==0) {
    std::cerr << "Invalid atom type in Pij.\n";
    return(0.0);
  }
  xx1_pow[0] = 0.0;
  xx1_pow[1] = 1;
  xx1_pow[2] = 2*XX1;
  xx1_pow[3] = 3*XX1*XX1;
  xx2_pow[0] = 1.0;
  xx2_pow[1] = XX2;
  xx2_pow[2] = XX2*XX2;
  xx2_pow[3] = XX2*XX2*XX2;
  clm_ptr = &CLM[KJ][NH][NC][0];
  ANSY=0;
  for (J=1; J<17; J++) {
    X = *++clm_ptr*xx1_pow[IN2[J][1]]*xx2_pow[IN2[J][2]];
    ANSY += X;
  }
  return ANSY;
}

double dspline_dNic(int KJ, double Nih, double Nic) {
  int J, NH, NC;
  double X, ANSY;
  double xx1_pow[4];
  double xx2_pow[4];
  double XX1, XX2;
  const double *clm_ptr;
  XX1 = Nih + 1.0;  //to make it match the paper
  XX2 = Nic + 1.0;  //likewise
  //bicubic interpolation or perhaps bicubic spline
  if (XX1 < 1) XX1 = 1;
  if (XX2 < 1) return 0.0;
  NH=(int)floor(XX1);
  NC=(int)floor(XX2);
  if(NH > MAX_BC_NEIGHBORS || NC > MAX_BC_NEIGHBORS) {
    //std::cerr << "Invalid neighbor numbers in Pij.\n";
    return(0.0);
  }
  if (KJ==0) {
    std::cerr << "Invalid atom type in Pij.\n";
    return(0.0);
  }
  xx1_pow[0] = 1.0;
  xx1_pow[1] = XX1;
  xx1_pow[2] = XX1*XX1;
  xx1_pow[3] = XX1*XX1*XX1;
  xx2_pow[0] = 0.0;
  xx2_pow[1] = 1.0;
  xx2_pow[2] = 2*XX2;
  xx2_pow[3] = 3*XX2*XX2;
  clm_ptr = &CLM[KJ][NH][NC][0];
  ANSY=0;
  for (J=1; J<17; J++) {
    X = *++clm_ptr*xx1_pow[IN2[J][1]]*xx2_pow[IN2[J][2]];
    ANSY += X;
  }
  return ANSY;
}  

double dPij_dNih(int KJ, double Nih, double Nic) {
  if(in2_initialized == 0) init_in2();
  //for fast computation
  int iNih = int(Nih);
  int iNic = int(Nic);
  bool iflag = 0;
  if(KJ < TYPE_MAX && 
     double(iNih)==Nih && iNih<N_MAX &&
     double(iNic)==Nic && iNic<N_MAX) iflag = 1;
  if(iflag && dPijdNihFlag[KJ][iNih][iNic])
    return dPijdNihVal[KJ][iNih][iNic];

  double return_val = dspline_dNih(KJ, Nih, Nic);
  
  if(iflag) {
    dPijdNihFlag[KJ][iNih][iNic] = 1;
    dPijdNihVal[KJ][iNih][iNic] = return_val; 
  }
  return return_val;
}
  
double dPij_dNic(int KJ, double Nih, double Nic) {
  if(in2_initialized == 0) init_in2();
  //for fast computation
  int iNih = int(Nih);
  int iNic = int(Nic);
  bool iflag = 0;
  if(KJ < TYPE_MAX && 
     double(iNih)==Nih && iNih<N_MAX &&
     double(iNic)==Nic && iNic<N_MAX) iflag = 1;
  if(iflag && dPijdNicFlag[KJ][iNih][iNic])
    return dPijdNicVal[KJ][iNih][iNic];

  double return_val = dspline_dNic(KJ, Nih, Nic);
  
  if(iflag) {
    dPijdNicFlag[KJ][iNih][iNic] = 1;
    dPijdNicVal[KJ][iNih][iNic] = return_val; 
  }
  return return_val;
}

#undef MAX_BC_NEIGHBS
#undef TYPE_MAX
#undef N_MAX

}/* tbtools */
