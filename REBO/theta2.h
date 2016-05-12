#ifndef _TBTOOLS_THETA2_H
#define _TBTOOLS_THETA2_H

namespace tbtools {
  
double theta2(const double* a, const double* b, const double* r);
  
void dtheta2_dr(const double* a, const double* b, const double* r, 
		double* result);

void dtheta2_da(const double* a, const double* b, const double* r,
		double* result);

void dtheta2_db(const double* a, const double* b, const double* r,
		double* result);

}/* tbtools */

#endif
