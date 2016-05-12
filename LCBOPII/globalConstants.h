#ifndef _LCBOPII_GLOBAL_CONSTANTS_H
#define _LCBOPII_GLOBAL_CONSTANTS_H

#include <limits>

namespace LCBOPII {

namespace Const {

/********** constants **********/
const double PI         = 3.14159265358979323846;
const double EPS_DOUBLE = std::numeric_limits<double>::epsilon();
const double SQRT3      = 1.7320508075688772935;

/********** VR & VA **********/
const double Asr   = 53026.92613715;
const double alpha = 6.74750993;
const double Bsr1  = 27618.35705738;
const double beta1 = 6.34503890;
const double Bsr2  = 34.07142502;
const double beta2 = 1.19712839;

/********** G **********/
const double gmin  = 0.0020588719;
const double ggr   = 0.0831047003;
const double gmax  = 16.0;
const double g1[3] = {0.7233666272, 1.7334665088, 1.8701997632};
const double g2[5] = {0.739945277952, -1.999211817114, -17.43251544992, 
		      -33.961271099802, -21.385299384112002};
const double g3[5] = {-15.19, -25.616855239779, -21.517283970381, 
		      0.98990809933, 13.664161601096};
const double Ay0   = -0.4;
const double By0   = 0.01875;
const double Ag    = 5.6304664723;
const double Bg    = 0.1516943990;
const double Cg    = 0.009832975891;
const double Dg    = -0.189175977654;
const double Eg    = 0.050977653631;

/********** H **********/
const double d     = 0.14;
const double d2    = d*d;
const double d3    = d*d2;
const double d4    = d*d3;
const double d5    = d*d4;
const double d6    = d*d5;
const double C1    = 3.335;
const double C12   = C1*C1; 
const double C4    = 220.0;
const double C6    = (-C12-12.0*C4*d2)/(30.0*d4);
const double L     = 1.0-C1*d+C12*d2/2.0+C4*d4+C6*d6;
const double R0    = 1.0+C1*d+C12*d2/2.0+C4*d4+C6*d6;
const double kappa = (C1-C12*d-4.0*C4*d3-6.0*C6*d5)/L;
const double R1    = C1+C12*d+4.0*C4*d3+6.0*C6*d5;

/********** Fconj **********/
const double Fconj0[4][4] = {
  { 0.00000000,  0.02068400, -0.00464600, -0.12779750},
  { 0.02068400,  0.00000000,  0.01858496, -0.10428500},
  {-0.00464600,  0.01858496,  0.00000000, -0.02733950},
  {-0.12779750, -0.10428500, -0.02733950,  0.00000000},
};

const double Fconj1[4][4] = {
  { 0.00000000,  0.05844900,  0.04163700, -0.12779750},
  { 0.05844900,  0.13789300,  0.09259214, -0.12428500},
  { 0.04163700,  0.09259214,  0.09355750, -0.03933950},
  {-0.12779750, -0.12428500, -0.03933950,  0.00000000},
};

/********** A **********/
const double alpha0 = -0.95;

/********** T **********/
const double At  = -52.61163955;
const double Bt1 = -0.1947358464;
const double Bt2 = 3.8;
const double Bt3 = 0.62;
const double Bt4 = 0.005;

/********** Vlr **********/
const double r0      = 3.715735;
const double eps2    = 0.002827918;
const double lambda1 = 1.338162;
const double lambda2 = 2.260479;
const double eps1    = eps2*lambda2*lambda2/(lambda1*lambda1);
const double v1      = eps1-eps2;
const double v2      = 0.0;

/********** Vmr **********/
const double rmr1 = 4.0;
const double rmr2 = 2.9;
const double Amr0 = -0.2345;
const double Amr1 = -0.67;
const double Amr2 = -4.94;
const double B    = 2.1;

}/* Const */

}/* LCBOPII */

#endif
