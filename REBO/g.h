/*-----------------------------g.h---------------------------- 
  This file contains the function prototypes for g.c.
  This software was written by Stephen Jordan and modeled on
  code from Zyvex as well as DW Brenner's original Fortran.
  ------------------------------------------------------------*/

#ifndef _TBTOOLS_G_H
#define _TBTOOLS_G_H

namespace tbtools {

double G(int type, double costh, double Nit);

double dG_dcosth(int type, double costh, double Nit);

double dG_dNit(int type, double costh, double Nit);

}/* tbtools */

#endif

