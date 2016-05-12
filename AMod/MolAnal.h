#ifndef _AMOD_MOL_ANAL_H
#define _AMOD_MOL_ANAL_H

#include "Molecule.h"
#include "../Util/Accessors.h"
#include <utility>
#include <vector>

namespace AMod {

class Molecule;
class MolAnal {
public:
  typedef std::vector<IDpos> Chain;
  typedef std::vector<Chain> Chains;
  static const int DEFAULT_MAX_CHAIN_LENGTH;

  MolAnal(int _maxChainLen = DEFAULT_MAX_CHAIN_LENGTH);
  void reset();
  static bool arrowNorm(const Arrow& ar, double* n, const double* ref = NULL);
  static bool bondNorm(const Bond& bond, double* n, const double* ref = NULL);
  static bool arrowNeighbs(const Arrow& ar, const double* axis, 
			   int& left, int& right);
  /****** find ring(s) ******/
  bool ring(const Arrow& ar, Chain& chain, const double* ref = NULL) const;
  void rings(const Molecule& mol, Chains& chains);
  template<class InputIterator>
  void rings(const Molecule& mol, 
	     InputIterator bond_id_first,
	     InputIterator bond_id_last,
	     Chains& chains);
  void rings(const Molecule& mol, int bondID, Chains& chainsA, Chains& chainsB);
  /****** find other pattern(s) ******/
  bool stair(const Arrow& ar, Chain& chain, const double* ref = NULL) const;
  /****** set & get ******/
  Util::Accessors<int> maxChainLength;

protected:
  static void arrowOrder(const Arrow& ar, const double* axis,
			 std::vector< std::pair<double,int> >& order);
  static double angle(const double* x, const double* z, const double* t);
  void regChain(Chain& chain);
  bool checkDup(const Chain& chain) const;

  struct BondMark {
    Chain* pc[2];
    BondMark();
  };

private:
  std::vector<BondMark> marks;
};

std::ostream& operator<< (std::ostream& fout, const MolAnal::Chain& chain);

std::ostream& operator<< (std::ostream& fout, const MolAnal::Chains& chains);

/************************************************************/
template<class InputIterator>
void MolAnal::rings(const Molecule& mol, 
		    InputIterator bond_id_first,
		    InputIterator bond_id_last,
		    Chains& chains) {
  double ref[3];
  int i, j;
  marks.assign(mol.nbonds(),BondMark());
  chains.assign(1,Chain());
  for(; bond_id_first != bond_id_last; bond_id_first++) {
    const Arrow& ar = mol.bond(*bond_id_first).arrow(0);
    if(!bondNorm(ar.bond(), ref, NULL)) continue;
    for(j = 0; j < 2; j++) {
      if(j == 1) { for(i = 0; i < 3; i++) ref[i] = -ref[i]; }
      if(ring(ar, chains.back(), ref) && 
	 !checkDup(chains.back())) {
	regChain(chains.back());
	chains.push_back(Chain());
      }
    }//end for(i...
  }//end for(k...
  chains.pop_back();
}

} /* AMod */

#endif
