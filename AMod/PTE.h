#ifndef _AMOD_PTE_H
#define _AMOD_PTE_H

namespace AMod {

class PTE {
public:
  enum AtomName {
    Unknown = 0, 
    Hydrogen, 
    Helium,
    Lithium,
    Beryllium,
    Boron,
    Carbon,
    Nitrogen,
    Oxygen,
    Fluorine,
    Neon,
    nAtomTypes};

  static double atomicMass(int _num);
  static void atomicMass(const int* _num, double* _mass, int _n);

protected:
  static const double AtomicMassArray[];
};
  
} /* AMod */

#endif
