#include "SWStimulus.h"
#include "../AMod/MolAnal.h"
#include "../Util/util.h"
#include "../Util/constants.h"
#include <set>
#include <string>

namespace SW {

SWStimulus::SWStimulus(int _bondID, bool _rotCCW) {
  set(_bondID, _rotCCW);
}

SWStimulus::SWStimulus(int _bondID, bool _rotCCW, const double _rotAxis[3]) {
  set(_bondID, _rotCCW, _rotAxis);
}

SWStimulus::SWStimulus(const SWStimulus& swStimulus) {
  *this = swStimulus;
}

SWStimulus& SWStimulus::operator= (const SWStimulus& swStimulus) {
  if(this != &swStimulus) {
    bondID = swStimulus.bondID;
    rotCCW = swStimulus.rotCCW;
    for(int k = 0; k < 3; k++) 
      rotAxis[k] = swStimulus.rotAxis[k];
  }
  return *this;
}

void SWStimulus::reset() {//reset to bad
  set();
}

void SWStimulus::set(int _bondID, bool _rotCCW) {
  bondID = _bondID;
  rotCCW = _rotCCW;
  for(int k = 0; k < 3; k++) 
    rotAxis[k] = 0.0;
}

void SWStimulus::set(int _bondID, bool _rotCCW, const double _rotAxis[3]) {
  bondID = _bondID;
  rotCCW = _rotCCW;
  for(int k = 0; k < 3; k++) 
    rotAxis[k] = _rotAxis[k];
}

bool SWStimulus::commit(const AMod::Molecule& _mol,
			int minRingSize, int maxRingSize) {
  int j;
  double tmpAxis[3];
  const double* ref;

  if(bondID < 0 || bondID >= _mol.nbonds()) return false;  
  ref = (invalidRotAxis() ? NULL : rotAxis);
  if(!AMod::MolAnal::bondNorm(_mol.bond(bondID), tmpAxis, ref)) return false;
  for(j = 0; j < 3; j++) rotAxis[j] = tmpAxis[j];
  
  if(minRingSize >= 0 && maxRingSize >= 0) 
    return check(_mol, minRingSize, maxRingSize);

  return true;
}

bool SWStimulus::check(const AMod::Molecule& _mol, 
		       int minRingSize, int maxRingSize) const {
  AMod::MolAnal::Chains chainsA, chainsB;

  AMod::MolAnal().rings(_mol, bondID, chainsA, chainsB);  
  for(int k = 0, kmax = chainsA.size(); k < kmax; k++) {
    int chain_sz = chainsA[k].size();
    chain_sz--;
    if(chain_sz < minRingSize || chain_sz > maxRingSize) return false;
  }
  for(int k = 0, kmax = chainsB.size(); k < kmax; k++) {
    int chain_sz = chainsB[k].size();
    chain_sz++;
    if(chain_sz < minRingSize || chain_sz > maxRingSize) return false;
  }  

  return true;
}

bool SWStimulus::invalidRotAxis() const {
  return (Util::vdotv(rotAxis, rotAxis, 3) < Util::EPS_DOUBLE);
}

int SWStimulus::sourceTxt(std::istream& fin) {
  std::string tword;

  if(!fin.good()) return 1;
  
  fin >> tword >> bondID;
  fin >> tword >> rotCCW;
  fin >> tword >> rotAxis[0] >> rotAxis[1] >> rotAxis[2];
  
  return (fin.good() ? 0 : 1);
}

void SWStimulus::dumpTxt(std::ostream& fout) const {
  fout << "bondID " << bondID << std::endl;
  fout << "rotCCW " << int(rotCCW) << std::endl;
  fout << "rotAxis " 
       << rotAxis[0] << " " << rotAxis[1] << " " << rotAxis[2] << std::endl;
}

std::istream& operator>> (std::istream& fin, SWStimulus& swStimulus) {
  swStimulus.sourceTxt(fin);
  return fin;
}

std::ostream& operator<< (std::ostream& fout, const SWStimulus& swStimulus) {
  swStimulus.dumpTxt(fout);
  return fout;
}

} /* SW */
