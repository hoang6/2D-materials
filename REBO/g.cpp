/*--------------------------------g.c------------------------------------
  This file contains routines to calculate g as in equations 8 and 11,
  as defined in DW Brenner J. Phys Cond. Mat 14 (2002) 783-802.  This is
  accomplished throught the use of 1 dimensional 5th degree polynomial
  interpolation.  This software was written by Stephen Jordan and modeled
  on code from Zyvex as well as DW Brenner's original Fortran.
  -----------------------------------------------------------------------*/

#include "g.h"
#include "globalConstants.h"
#include <cmath>  //for cos and floor

namespace tbtools {

/* The next table implements the special cases for G sub C (cos (theta))
   described on pages 11 through 13.  The index to igc is cos (theta) suitably
   quantized, and the output is an index to SPGC, so it says which sixth order
   polynomial spline we will use. */
/* 16*4,2*3,2*2,5*1 */
static const int igc[26] = {0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                            4, 4, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1 };

/* The next table implements the special cases for G sub H (cos (theta))
   described in the middle of page 17.  The index to igh is cos (theta)
   suitably quantized, and the output is an index to SPGH, so it says which
   sixth order polynomial spline we will use. */
/* 18*3,4*2,3*1 */
static const int igh[26] = {0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                            3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1 };


//The first five groups contain the coefficients for G sub c
//The last group contains the coefficients for gamma sub c
static const double SPGC[5+1][6+1] = {
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0.2817216000000E+00, 0.1062912000000E+01,  0.2136736000000E+01 ,
      0.2533952000000E+01, 0.1554736000000E+01,  0.3863296000000E+00},
  {0, 0.2817216000000E+00, 0.1062912000000E+01,  0.2136736000000E+01,
      0.2533952000000E+01, 0.1554736000000E+01,  0.3863296000000E+00},
  {0, 0.6900668660000E+00, 0.5460691360000E+01,  0.2301345680000E+02,
      0.5491519344000E+02, 0.6862037040000E+02,  0.3470897779200E+02},
  {0, 0.3754490870000E+00, 0.1407252749388E+01,  0.2255103926323E+01,
      0.2028902219952E+01, 0.1426981217906E+01,  0.5063107994308E+00}, 
  {0, 0.2718558000000E+00, 0.4892727456293E+00, -0.4328199017473E+00,
     -0.5616795197048E+00, 0.1270874966906E+01, -0.3750409108350E-01}
};

//This contains the coefficients for G sub h    
static const double SPGH[3+1][6+1] = {
  {0, 0, 0, 0, 0, 0, 0},
  {0, 270.467795364007301, 1549.701314596994564, 3781.927258631323866,
     4582.337619544424228, 2721.538161662818368,  630.658598136730774},
  {0,  16.956325544514659,  -21.059084522755980, -102.394184748124742,
     -210.527926707779059, -229.759473570467513,  -94.968528666251945},
  {0,  19.065031149937783,    2.017732531534021,    -2.566444502991983,
        3.291353893907436,   -2.653536801884563,     0.837650930130006}
};


//this implements equation 12
double Qi(double Nit) {
  if(Nit < 3.2) return 1;
  if(Nit > 3.7) return 0;
  //otherwise:
  return (1.0 + std::cos(2*PI*(Nit-3.2)))/2;
}

double G(int type, double costh, double Nit) {
  int ig, ig1;
  double qi, gangle = 0.0, gangle1 = 0.0;
  double costh2 = costh*costh;
  double costh3 = costh2*costh;
  double costh4 = costh3*costh;
  double costh5 = costh4*costh;

  //  printf("G: %i\t%e\t%e\t", type, costh, Nit);
  if(type == 1) { //carbon
    ig = igc[(int)std::floor(-costh*12.0)+13];
    gangle = SPGC[ig][1]
           + SPGC[ig][2]*costh
           + SPGC[ig][3]*costh2
           + SPGC[ig][4]*costh3
           + SPGC[ig][5]*costh4
           + SPGC[ig][6]*costh5;
    qi = Qi(Nit);
    if(ig == 4 && qi != 0.0) {
      // Now we're revising g sub C as in equation 11
      ig1 = ig + 1; 
      // ig was 4, so ig1 is 5.  This is the row of SPGC that defines gamma sub C
      //gangle1 is gamma sub c, and gangle is currently capital G
      gangle1 = SPGC[ig1][1]
	      + SPGC[ig1][2]*costh
	      + SPGC[ig1][3]*costh2
              + SPGC[ig1][4]*costh3
              + SPGC[ig1][5]*costh4
	      + SPGC[ig1][6]*costh5;
      // Next line is equation 11
      gangle += qi*(gangle1-gangle); //now gangle is lowercase g
    }//end if(ig == 4)
  }//end if(type == 1)
  else if(type == 2) { //hydrogen
    ig = igh[(int)std::floor(-costh*12.0)+13];
    gangle = SPGH[ig][1]
           + SPGH[ig][2]*costh
           + SPGH[ig][3]*costh2
           + SPGH[ig][4]*costh3
           + SPGH[ig][5]*costh4
           + SPGH[ig][6]*costh5;
  }//end else if(type == 2)
  //  printf("%e\n", gangle);
  return gangle;
}

double dQi_dNit(double Nit) {
  if(Nit < 3.2 || Nit > 3.7) return 0.0;
  return -0.5*2*PI*sin(2*PI*(Nit-3.2));
}

double dG_dcosth(int type, double costh, double Nit) {
  int ig, ig1;
  double gangle, gangle1;
  double costh2 = costh*costh;
  double costh3 = costh2*costh;
  double costh4 = costh3*costh;

  gangle = 0.0;
  //  printf("G: %i\t%e\t%e\t", type, costh, Nit);
  if(type == 1) { //carbon
    ig = igc[(int)std::floor(-costh*12.0)+13];
    gangle = 0.0
           + SPGC[ig][2]
           + SPGC[ig][3]*2*costh
           + SPGC[ig][4]*3*costh2
           + SPGC[ig][5]*4*costh3
           + SPGC[ig][6]*5*costh4;
    if(ig == 4) {
      // Now we're revising g sub C as in equation 11
      ig1 = ig + 1; 
      // ig was 4, so ig1 is 5.  This is the row of SPGC that defines gamma sub C
      //gangle1 is gamma sub c, and gangle is currently capital G
      gangle1 = 0.0
              + SPGC[ig1][2]
              + SPGC[ig1][3]*2*costh
              + SPGC[ig1][4]*3*costh2
              + SPGC[ig1][5]*4*costh3
	      + SPGC[ig1][6]*5*costh4;
      // Next line is equation 11
      gangle += Qi(Nit)*(gangle1-gangle); //now gangle is lowercase g
    }
  }
  if(type == 2) { //hydrogen
    ig = igh[(int)std::floor(-costh*12.0)+13];
    gangle = 0.0
           + SPGH[ig][2]*
           + SPGH[ig][3]*2*costh
           + SPGH[ig][4]*3*costh2
           + SPGH[ig][5]*4*costh3
           + SPGH[ig][6]*5*costh4;
  }
  //  printf("%e\n", gangle);
  return gangle;
}

double dG_dNit(int type, double costh, double Nit) {
  int ig, ig1;
  double gangle, gangle1;
  double costh2 = costh*costh;
  double costh3 = costh2*costh;
  double costh4 = costh3*costh;
  double costh5 = costh4*costh;

  gangle = 0.0;
  //  printf("G: %i\t%e\t%e\t", type, costh, Nit);
  if(type == 1) { //carbon
    ig = igc[(int)std::floor(-costh*12.0)+13];
    gangle = 0.0;
    if(ig == 4) {
      // Now we're revising g sub C as in equation 11
      ig1 = ig + 1; 
      // ig was 4, so ig1 is 5.  This is the row of SPGC that defines gamma sub C
      //gangle1 is gamma sub c, and gangle is currently capital G
      gangle = SPGC[ig][1]
	     + SPGC[ig][2]*costh
             + SPGC[ig][3]*costh2
             + SPGC[ig][4]*costh3
             + SPGC[ig][5]*costh4
             + SPGC[ig][6]*costh5;
      gangle1 = SPGC[ig1][1]
              + SPGC[ig1][2]*costh
              + SPGC[ig1][3]*costh2
              + SPGC[ig1][4]*costh3
              + SPGC[ig1][5]*costh4
	      + SPGC[ig1][6]*costh5;
      // Next line is equation 11
      gangle = dQi_dNit(Nit)*(gangle1-gangle); //now gangle is lowercase g
    }
  }
  if(type == 2) { //hydrogen
    gangle = 0.0;
  }
  //  printf("%e\n", gangle);
  return gangle;  
}

}/* tbtools */
