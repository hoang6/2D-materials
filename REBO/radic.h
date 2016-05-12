/*----------------------------radic.h------------------------------
  This file contains the function prototypes for Fij and Tij, which
  are defined in radic.c.  This software was written by Stephen
  Jordan.
  -----------------------------------------------------------------*/

#ifndef _TBTOOLS_RADIC_H
#define _TBTOOLS_RADIC_H

namespace tbtools {

double Fij(int KI, int KJ, double XNT1, double XNT2, double CONJUG);

double Tij(double XNT1, double XNT2, double CONJUG);

double dFij_dXNT1(int KI, int KJ, double XNT1, double XNT2, double CONJUG);

double dFij_dXNT2(int KI, int KJ, double XNT1, double XNT2, double CONJUG);
  
double dFij_dCONJUG(int KI, int KJ, double XNT1, double XNT2, double CONJUG);

double dTij_dXNT1(double XNT1, double XNT2, double CONJUG);

double dTij_dXNT2(double XNT1, double XNT2, double CONJUG);

double dTij_dCONJUG(double XNT1, double XNT2, double CONJUG);

}/* tbtools */

#endif
