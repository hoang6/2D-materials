#include "PTE.h"

namespace AMod {

const double PTE::AtomicMassArray[] = {
  0.0,
  1.00794,
  4.002602,
  6.941,
  9.012182,
  10.811,
  12.0107,
  14.0067,
  15.9994,
  18.9984032,
  20.1797
  //...
};
  
double PTE::atomicMass(int _num) { return AtomicMassArray[_num]; }

void PTE::atomicMass(const int* _num, double* _mass, int _n) {
  for(int k = 0; k < _n; k++) _mass[k] = atomicMass(_num[k]);
}

} /* AMod */
