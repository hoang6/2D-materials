/*------------------------------constants.c------------------------------------
  This file contains the constant parameters of Brenner's potential. This
  software was written by Stephen Jordan and modeled on code from Zyvex as well
  as DW Brenner's original Fortran.
  ----------------------------------------------------------------------------*/

#include "constants.h"
#include "globalConstants.h"
#include <cmath>      //for exp
#include <iostream>

namespace tbtools {

//Dmax
static double DmaxVal[1+2][1+2] = {
  {0.0, 0.0, 0.0},
  {0.0, 2.0, 1.8},
  {0.0, 1.8, 1.7}
};
//Dmin
static double DminVal[1+2][1+2] = {
  {0.0, 0.0, 0.0},
  {0.0, 1.7, 1.3},
  {0.0, 1.3, 1.1}
};  
//Q as in equation 5 of Brenner's new paper
static double QVal[1+2][1+2] = {
  {0.0, 0.0, 0.0},
  {0.0, 0.3134602960833, 0.340775728},
  {0.0, 0.340775728, 0.370471487045}
};  
//A as in equation 5 of Brenner's new paper
static double AVal[1+2][1+2] = {
  {0.0, 0.0, 0.0},
  {0.0, 10953.544162170, 149.94098723},
  {0.0, 149.94098723, 32.817355747}
};
//alpha as in equation 5 of Brenner's new paper
//(not the same as alpha in the old paper)
static double alphaVal[1+2][1+2] = {
  {0.0, 0.0, 0.0},
  {0.0, 4.7465390606595, 4.10254983},
  {0.0, 4.10254983, 3.536298648}
};
static double XDBVal[1+2][1+2][1+2];
static double REGVal[1+2][1+2][1+2];
static double BVal[1+2][1+2][1+3];
static double betaVal[1+2][1+2][1+3]; 
static int initialized = 0;

/*This is taken from Brenner's software, it has no counterpart in
  any of his papers, but REG*exp(XDB) serves a purpose analogous
  to the exponential term exp(lambda sub ijk) in equation 8 of
  Brenner's new paper.*/
double XDBGen(int type1, int type2, int type3) {
  if(type1 == 2 && type2 == 2 && type3 == 2) return 4.0;
  if(type1 == 2 && type2 == 1 && type3 == 2) return 4.0;
  if(type1 == 2 && type2 == 2 && type3 == 1) return 4.0;
  if(type1 == 2 && type2 == 1 && type3 == 1) return 4.0;
  //otherwise:
  return 0.0;
}


double REGGen(int type1, int type2, int type3) {
  const double RHH = 0.7415886997;
  const double RCH = 1.09;
  if(type1 == 2 && type2 == 1 && type3 == 2) 
    return std::exp(XDBGen(type1, type2, type3)*(RHH-RCH));
  if(type1 == 2 && type2 == 2 && type3 == 1) 
    return std::exp(XDBGen(type1, type2, type3)*(RCH-RHH));
  //otherwise:
  return 1.0;
}


//this is B as used in equation 6 of Brenner's new paper
double BGen(int type1, int type2, int n) {
  if(type1 == 1 && type2 == 1) {//C-C
    if(n == 1) return 12388.79197798;
    if(n == 2) return 17.56740646509;
    if(n == 3) return 30.71493208065;
    //otherwise
    std::cerr << "Error in B, invalid n.\n";
    return 0.0;
  }
  if((type1 == 1 && type2 == 2)||(type1 == 2 && type2 == 1)) { //C-H or H-C
    if(n == 1) return 32.3551866587;
    if(n == 2 || n == 3) return 0.0;
    //otherwise:
    std::cerr << "Error in B, invalid n.\n";
    return 0.0;
  }
  if(type1 == 2 && type2 == 2) { //H-H
    if(n == 1) return 29.632593;
    if(n == 2 || n == 3) return  0.0;
    //otherwise
    std::cerr << "Error in B, invalid n.\n";
    return 0.0;
  }
  //otherwise:
  std::cerr << "Error in B, invalid type.\n";
  return 0.0;
}


//this is beta as used in equation 6 of Brenner's new paper  
double betaGen(int type1, int type2, int n) {
  if(type1 == 1 && type2 == 1) {//C-C
    if(n == 1) return 4.7204523127;
    if(n == 2) return 1.4332132499;
    if(n == 3) return 1.3826912506;
    //otherwise
    std::cerr << "Error in beta, invalid n.\n";
    return 0.0;
  }
  if((type1 == 1 && type2 == 2)||(type1 == 2 && type2 == 1)) { //C-H or H-C
    if(n == 1) return 1.43445805925;
    if(n == 2 || n == 3) return 0.0;
    //otherwise:
    std::cerr << "Error in beta, invalid n.\n";
    return 0.0;
  }
  if(type1 == 2 && type2 == 2) { //H-H
    if(n == 1) return 1.71589217;
    if(n == 2 || n == 3) return  0.0;
    //otherwise
    std::cerr << "Error in beta, invalid n.\n";
    return 0.0;
  }
  //otherwise:
  std::cerr << "Error in beta, invalid type.\n";
  return 0.0;
}


void constantsGen() {
  if(initialized) return;
  int type1, type2, type3, n;
  for(type1 = 1; type1 <= 2; type1++) 
    for(type2 = 1; type2 <= 2; type2++)
      for(type3 = 1; type3 <= 2; type3++) {
	XDBVal[type1][type2][type3] = XDBGen(type1, type2, type3);
	REGVal[type1][type2][type3] = REGGen(type1, type2, type3);
      }
  for(type1 = 1; type1 <= 2; type1++) 
    for(type2 = 1; type2 <= 2; type2++)
      for(n = 1; n <= 3; n++) {
	BVal[type1][type2][n] = BGen(type1, type2, n);
	betaVal[type1][type2][n] = betaGen(type1, type2, n);	
      }
  initialized = 1;
}


double Dmax(int type1, int type2) { return DmaxVal[type1][type2]; }


double Dmin(int type1, int type2) { return DminVal[type1][type2]; }


double XDB(int type1, int type2, int type3) {
  return XDBVal[type1][type2][type3];
}


double REG(int type1, int type2, int type3) {
  return REGVal[type1][type2][type3];
}


double Q(int type1, int type2) { return QVal[type1][type2]; }


double A(int type1, int type2) { return AVal[type1][type2]; }


double alpha(int type1, int type2) { return alphaVal[type1][type2]; }


double B(int type1, int type2, int n) { return BVal[type1][type2][n]; }


double beta(int type1, int type2, int n) { return betaVal[type1][type2][n]; }


double fc(double r, int type1, int type2) {
  double dmin = DminVal[type1][type2];
  double dmax = DmaxVal[type1][type2];
  if(r < dmin) return 1.0;
  else if(r > dmax) return 0.0;
  return 0.5*(1.0+cos(PI*(r-dmin)/(dmax-dmin)));
}


double dfc_dr(double r, int type1, int type2) {
  double dmin = DminVal[type1][type2];
  double dmax = DmaxVal[type1][type2];
    if(r < dmin || r > dmax) return 0.0;
  return -0.5*(PI/(dmax-dmin))*sin(PI*(r-dmin)/(dmax-dmin));
}


double F(double N) {
  if(N > 3.0) return 0.0;
  else if(N < 2.0) return 1.0;
  return 0.5*(1.0+cos(PI*(N-2.0)));
}


double dF_dN(double N) {
  if(N < 2.0 || N > 3.0) return 0.0;
  return -0.5*PI*sin(PI*(N-2.0));
}

}/* tbtools */
