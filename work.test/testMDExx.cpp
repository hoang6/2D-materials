#include "../MD/Sim.h"
#include "../Util/util.h"
#include "../Util/RNG.h"
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 3 && argc != 4) {
    cout << "Usage: testMDExx <input> <Exx> [seed]" << endl;
    return 1;
  }

  /********** Parse Argument **********/
  MD::IO mdIO;
  ifstream finput;
  double exx;
  long seed;

  finput.open(argv[1]);
  if(mdIO.readInput(finput) != 0) {
    cout << "Error: fail to read input file: " << argv[1] << endl;
    return 2;
  }
  finput.close();
  mdIO.echoInput(cout);
  cout << "--------------------" << endl;
  if(mdIO.ensemble != "NVT") {
    cout << "Error: must use NVT ensemble" << endl;
    return 3;
  }

  if(!Util::readWord(argv[2], exx)) {
    cout << "Error: fail to read Exx: " << argv[2] << endl;
    return 4;
  }

  seed = time(NULL);
  if(argc == 4) Util::readWord(argv[3], seed);
  cout <<  "SEED            " << seed << endl;
  Util::RNG::seed(seed);

  /********** Init MD **********/
  MD::Data mdData;
  
  if(mdIO.getData(mdData) != 0) {
    cout << "Error: fail to get MDData" << endl;
    return 5;
  }

  //Apply exx
  mdData.mol.axis(0).x[0] *= (1+exx);
  for(int k = 0; k < mdData.mol.natoms(); k++)
    mdData.mol.atom(k)->x[0] *= (1+exx);  
  mdData.mol.topo().update(AMod::CASE_BRUTAL);

  /********** Run MD **********/
  MD::Sim sim;
  MD::NVT mdEnsem;
  MD::NVT::Para mdPara;
  
  int exitCode;

  if(mdIO.getPara(mdPara) != 0) {
    cout << "Error: fail to read MDNVT::Para" << endl;
    return 6;
  }

  mdData.resetTotalMomentum();
  mdEnsem.init(mdData,mdPara);
  exitCode = sim.run(mdEnsem,mdIO);
  mdEnsem.final();
  if(exitCode != 0) {
    cout << "Error: fail to run MD, exitCode = " << exitCode << endl;
    return 7;
  }  

  return 0;
}
