#include "IO.h"
#include "../Util/util.h"
#include "../Util/RNG.h"
#include "../Util/constants.h"
#include <sstream>

namespace MD {

IO::IO() { reset(); }

void IO::reset() {
  ensemble = "";
  ensemblePara.clear();
  dtime = 0.0;
  maxStep = 0;
  dumpStep = 0;
  inMol = "";
  inVel = "";
  initTemp = 0.0;
  outTag = "";
  foutMol.clear();
  foutVel.clear();
  foutEne.clear();
  keys.clear();
}

int IO::readInput(std::istream& _fin) {
  std::string line, keyword, tmp;
  std::istringstream iss;

  reset();

  if(!_fin.good()) return 1;
  while(std::getline(_fin,line).good()) {
    iss.clear(); iss.str(line);
    keyword = ""; iss >> keyword;
    if(keyword == "POTENTIAL") {
      keys.insert(keyword);
      iss >> potential; if(iss.fail()) return 2;
    }
    else if(keyword == "ENSEMBLE") {
      keys.insert(keyword);
      iss >> ensemble;
      while(iss.good()) { 
	iss >> tmp;
	if(!iss.fail())
	  ensemblePara.push_back(tmp);
      }
    }
    else if(keyword == "DTIME") {
      keys.insert(keyword);
      iss >> dtime; if(iss.fail()) return 3;
    }
    else if(keyword == "MAX_STEP") {
      keys.insert(keyword);
      iss >> maxStep; if(iss.fail()) return 4;
    }
    else if(keyword == "DUMP_STEP") {
      keys.insert(keyword);
      iss >> dumpStep; if(iss.fail()) return 5;
    }
    else if(keyword == "IN_MOL") {
      keys.insert(keyword);
      iss >> inMol; if(iss.fail()) return 6;
    }
    else if(keyword == "IN_VEL") {
      keys.insert(keyword);
      iss >> inVel; if(iss.fail()) return 7;
    }
    else if(keyword == "INIT_TEMP") {
      keys.insert(keyword);
      iss >> initTemp; if(iss.fail()) return 8;
    }
    else if(keyword == "OUT_TAG") {
      keys.insert(keyword);
      iss >> outTag; if(iss.fail()) return 9;
    }
    else if(keyword == "END") {
      keys.insert(keyword);
      break;
    }
  }

  return (findKey("END")?0:10);
}

void IO::echoInput(std::ostream& _fout) const {
  if(findKey("POTENTIAL")) {
    _fout << "POTENTIAL       " << potential << std::endl;
  }
  if(findKey("ENSEMBLE")) {
    _fout << "ENSEMBLE        " << ensemble << " ";
    for(int k = 0, kmax = ensemblePara.size(); k < kmax; k++)
      _fout << ensemblePara[k] << " ";
    _fout << std::endl;
  }
  if(findKey("DTIME"))
    _fout << "DTIME           " << dtime << std::endl;
  if(findKey("MAX_STEP"))
    _fout << "MAX_STEP        " << maxStep << std::endl;
  if(findKey("DUMP_STEP"))
    _fout << "DUMP_STEP       " << dumpStep << std::endl;
  if(findKey("IN_MOL"))
    _fout << "IN_MOL          " << inMol << std::endl;
  if(findKey("IN_VEL"))
    _fout << "IN_VEL          " << inVel << std::endl;
  if(findKey("INIT_TEMP"))
    _fout << "INIT_TEMP       " << initTemp << std::endl;
  if(findKey("OUT_TAG"))
    _fout << "OUT_TAG         " << outTag << std::endl;
  if(findKey("END"))
    _fout << "END" << std::endl;
}

int IO::getData(Data& _data) const {
  std::ifstream finMol;
  std::ifstream finVel;
  int i, natoms, ndata;

  _data.reset();
  
  if(!findKey("POTENTIAL")) return 1;
  if(potential == "REBO") _data.pot = AMod::REBO;
  else if(potential == "LCBOPIIN") _data.pot = AMod::LCBOPIIN;
  else return 2;

  if(!findKey("IN_MOL")) return 3;
  if(!openFile(inMol,finMol)) return 4;
  finMol >> _data.mol;
  natoms = _data.mol.natoms();
  ndata = 3*natoms;
  _data.vel.resize(ndata);
  if(natoms == 0) return 5;
  
  if(findKey("IN_VEL")) {
    if(!openFile(inVel,finVel)) return 6;
    for(i = 0; i < ndata; i++) finVel >> _data.vel[i];
    if(finVel.fail()) return 7;
  }
  else if(findKey("INIT_TEMP")) {
    for(i = 0; i < ndata; i++)
      _data.vel[i] = 2.0*Util::RNG::uniform_SFMT()-1.0;
    _data.resetTotalMomentum();
    _data.rescaleVel(initTemp*Util::BOLTZMANN_CONSTANT);
  }
  else return 8;

  return 0;
}

int IO::getPara(NVT::Para& _para) const {
  _para.reset();
  if(ensemble != "NVT") return 1;
  int n = ensemblePara.size();
  if(n >= 1) {
    if(!Util::readWord(ensemblePara[0],_para.kT)) return 2;
    _para.kT *= Util::BOLTZMANN_CONSTANT;
  }
  if(n >= 2)
    if(!Util::readWord(ensemblePara[1],_para.tNH)) return 3;
  return 0;
}

int IO::getPara(NVE1He::Para& _para) const {
  _para.reset();
  if(ensemble != "NVE1He") return 1;
  int n = ensemblePara.size();
  if(n >= 3)
    for(int i = 0; i < 3; i++)
      if(!Util::readWord(ensemblePara[i],_para.posHe[i])) return 2;
  if(n >= 6)
    for(int i = 0; i < 3; i++)
      if(!Util::readWord(ensemblePara[i+3],_para.velHe[i])) return 3;
  return 0;
}

int IO::openFiles() {
  if(findKey("OUT_TAG")) {
    if(!openFile(outTag+".mol",foutMol)) return 1;
    if(!openFile(outTag+".vel",foutVel)) return 2;
    if(!openFile(outTag+".ene",foutEne)) return 3;
    setFormat(foutMol,12);
    setFormat(foutVel,12);
    setFormat(foutEne,12);
  }

  return 0;
}

void IO::closeFiles() {
  if(findKey("OUT_TAG")) {
    foutMol.close();
    foutVel.close();
    foutEne.close();
  }
}

void IO::dumpBlank() {
  foutMol << std::endl;
  foutVel << std::endl;
  foutEne << std::endl;
}

bool IO::findKey(const std::string& _key) const {
  return (keys.find(_key) != keys.end());
}

void IO::setFormat(std::ostream& _fout, int _precision) {
  _fout.setf(std::ios::scientific, std::ios::floatfield);
  _fout.precision(_precision);
  _fout << std::left;
}

}/* MD */
