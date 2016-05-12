#ifndef _SW_SW_STIMULUS_H
#define _SW_SW_STIMULUS_H

#include "../AMod/Molecule.h"
#include <vector>
#include <iostream>

namespace SW {

class SWStimulus;
typedef std::vector<SWStimulus> SWStimuli;

class SWStimulus {
public:
  int bondID;
  bool rotCCW;
  double rotAxis[3];

  SWStimulus(int _bondID = -1, bool _rotCCW = false);
  SWStimulus(int _bondID, bool _rotCCW, const double _rotAxis[3]);
  SWStimulus(const SWStimulus& swStimulus);
  SWStimulus& operator= (const SWStimulus& swStimulus);
  void reset(); //reset to bad
  void set(int _bondID = -1, bool _rotCCW = false);
  void set(int _bondID, bool _rotCCW, const double _rotAxis[3]);
  bool commit(const AMod::Molecule& _mol, 
	      int minRingSize = 5, int maxRingSize = 8);
  int sourceTxt(std::istream& fin);
  void dumpTxt(std::ostream& fout) const;

protected:
  bool check(const AMod::Molecule& _mol, 
	     int minRingSize, int maxRingSize) const;
  bool invalidRotAxis() const;
};

std::istream& operator>> (std::istream& fin, SWStimulus& swStimulus);

std::ostream& operator<< (std::ostream& fout, const SWStimulus& swStimulus);

} /* SW */

#endif
