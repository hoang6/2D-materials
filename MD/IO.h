#ifndef _MD_IO_H
#define _MD_IO_H

#include "NVE.h"
#include "NVT.h"
#include "NVE1He.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

namespace MD {

class IO {
public:
  IO();
  void reset();
  int readInput(std::istream& _fin);
  void echoInput(std::ostream& _fout) const;
  int getData(Data& _data) const;
  int getPara(NVT::Para& _para) const;
  int getPara(NVE1He::Para& _para) const;
  int openFiles();
  void closeFiles();
  void dumpBlank();
  bool findKey(const std::string& _key) const;

  std::string potential;
  std::string ensemble;
  std::vector<std::string> ensemblePara;
  double dtime;
  int maxStep;
  int dumpStep;
  std::string inMol;
  std::string inVel;
  double initTemp;
  std::string outTag;
  std::ofstream foutMol;
  std::ofstream foutVel;
  std::ofstream foutEne;

private:
  template<class FSTREAM>
  static bool openFile(const std::string& _fname, FSTREAM& _file);
  void setFormat(std::ostream& _fout, int _precision);
  
  std::set<std::string> keys;
};

/************************************************************/
template<class FSTREAM>
bool IO::openFile(const std::string& _fname, FSTREAM& _file) {
  _file.open(_fname.c_str());
  return _file.good();
}

}/* MD */

#endif

