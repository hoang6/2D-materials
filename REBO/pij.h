/*---------------------------Pij.h---------------------------- 
  This file contains the function prototypes for Pij.c.
  This software was written by Stephen Jordan and modeled
  on code from Zyvex as well as DW Brenner's original Fortran.
  ------------------------------------------------------------*/

#ifndef _TBTOOLS_PIJ_H
#define _TBTOOLS_PIJ_H

namespace tbtools {

double Pij(int KJ, double Nih, double Nic);

double dPij_dNih(int KJ, double Nih, double Nic);
  
double dPij_dNic(int KJ, double Nih, double Nic);

}/* tbtools */

#endif
