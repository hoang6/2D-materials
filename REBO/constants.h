/*------------------------------constants.h------------------------------------
  This contains the function prototypes for constant.c.  Furthere documentation
  can be found there.  This software was written by Stephen Jordan.
  ---------------------------------------------------------------------------*/

#ifndef _TBTOOLS_CONSTANTS_H
#define _TBTOOLS_CONSTANTS_H

namespace tbtools {

void constantsGen();

double Dmax(int type1, int type2);

double Dmin(int type1, int type2);

double F(double N);

double dF_dN(double N);

double XDB(int type1, int type2, int type3);

double REG(int type1, int type2, int type3);

double Q(int type1, int type2);

double A(int type1, int type2);

double alpha(int type1, int type2);

double B(int type1, int type2, int n);

double beta(int type1, int type2, int n);

double fc(double r, int type1, int type2);

double dfc_dr(double r, int type1, int type2);

} /* tbtools */

#endif

