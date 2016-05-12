#include "../MD/Sim.h"
#include "../Util/RNG.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <ctime>
using namespace std;

class MDSimNVE1He: public MD::Sim {
public:
  virtual bool stop(MD::NVE& _ensem, int _step) {
    static bool flag = false;
    MD::NVE1He& _ensem1 = static_cast<MD::NVE1He&>(_ensem);
    if(!flag && _ensem1.eHeC() > 1.0e-3) flag = true;
    if(flag && _step>=1000 && _ensem1.eHeC()==0.0) return true;
    return false;
  }
};

class MDSimNVE: public MD::Sim {
public:
  virtual bool stop(MD::NVE& _ensem, int _step) {
    MD::NVT& _ensem1 = static_cast<MD::NVT&>(_ensem);
    double eKMax = _ensem1.eKinMax();
    double eKMean = _ensem1.eKin()/_ensem1.natoms();
    if(_step>=1000 && Util::abs(eKMax-eKMean)/eKMean < 0.05) return true;
    return false;
  }
};

int main(int argc, char* argv[]) {
  if(argc != 3 && argc != 4) {
    cout << "Usage: testMDRadHe <input> <ncycle> [seed]" << endl;
    return 1;
  }

  /********** Parse Arguments **********/
  MD::IO mdIOA, mdIOB, mdIOC;
  int ncycle;
  ifstream finput;
  long seed;

  finput.open(argv[1]);
  if(mdIOA.readInput(finput) != 0 || mdIOA.ensemble != "NVE1He") {
    cout << "Error: fail to read input file (NVE1He): " << argv[1] << endl;
    return 2;
  }
  if(!mdIOA.findKey("IN_VEL")) {
    cout << "Error: missing IN_VEL" << endl;
    return 3;
  }
  mdIOA.echoInput(cout);

  if(mdIOB.readInput(finput) != 0 || mdIOB.ensemble != "NVE") {
    cout << "Error: fail to read input file (NVE): " << argv[1] << endl;
    return 4;
  }
  mdIOB.echoInput(cout);

  if(mdIOC.readInput(finput) != 0 || mdIOC.ensemble != "NVT") {
    cout << "Error: fail to read input file (NVT): " << argv[1] << endl;
    return 5;
  }
  mdIOC.echoInput(cout);
  finput.close();
  cout << "--------------------" << endl;

  if(!Util::readWord(argv[2],ncycle)) {
    cout << "Error: fail to read <ncycle>" << endl;
    return 6;
  }

  seed = time(NULL);
  if(argc == 4) Util::readWord(argv[3], seed);
  cout <<  "SEED            " << seed << endl;
  Util::RNG::seed(seed);

  /********** Get MD Data **********/
  MD::Data mdData;

  if(mdIOA.getData(mdData) != 0) {
    cout << "Error: fail to get MD::Data" << endl;
    return 7;
  }

  /********** Run MD **********/
  MDSimNVE1He simA; MD::NVE1He mdEnsemA; MD::NVE1He::Para mdParaA;
  MDSimNVE1He simB; MD::NVE    mdEnsemB;
  MD::Sim     simC; MD::NVT    mdEnsemC; MD::NVT::Para    mdParaC;
  std::string tagA, tagB, tagC;
  double mLx, mLy;
  char buff[256];
  int cycle;
  int exitCode;
  
  if(mdIOA.getPara(mdParaA) != 0) {
    cout << "Error: fail to get MDPara (A)" << endl;
    return 8;
  }
  if(mdIOC.getPara(mdParaC) != 0) {
    cout << "Error: fail to get MDPara (C)" << endl;
    return 9;
  }

  mLx = mdData.mol.axis(0).x[0];
  mLy = mdData.mol.axis(1).x[1];  
  tagA = mdIOA.outTag;
  tagB = mdIOB.outTag;
  tagC = mdIOC.outTag;
  exitCode = 0;

  for(cycle = 1; cycle <= ncycle; cycle++) {
    std::sprintf(buff,"%03d",cycle);
    cout << "cycle = " << cycle << endl; cout.flush();
    cout << "Phase: Radiation(He)" << endl;
    mdParaA.posHe[0] = Util::RNG::uniform_SFMT()*mLx;
    mdParaA.posHe[1] = Util::RNG::uniform_SFMT()*mLy;
    mdIOA.outTag = tagA+buff;
    mdEnsemA.init(mdData,mdParaA);
    exitCode = simA.run(mdEnsemA,mdIOA);
    mdEnsemA.final();
    mdData.rmFreeC(simA.freeAtomIDs());
    cout << "Phase: Relaxation" << endl;
    mdIOB.outTag = tagB+buff;
    mdData.resetTotalMomentum();
    mdEnsemB.init(mdData);
    exitCode = simB.run(mdEnsemB,mdIOB);
    mdEnsemB.final();
    mdData.rmFreeC(simB.freeAtomIDs());
    cout << "Phase: Anneal" << endl;
    mdIOC.outTag = tagC+buff;
    mdData.resetTotalMomentum();
    mdEnsemC.init(mdData,mdParaC);
    exitCode = simC.run(mdEnsemC,mdIOC);
    mdEnsemC.final();
    mdData.rmFreeC(simC.freeAtomIDs());
  }

  if(exitCode != 0) {
    cout << "Error: fail to run MD, exitCode = " << exitCode << endl;
    return 10;
  }

  /********** Final MD **********/
  return 0;
}

