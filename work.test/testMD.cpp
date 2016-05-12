#include "../MD/Sim.h"
#include "../Util/RNG.h"
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 2 && argc != 3) {
    cout << "Usage: testMD <input> [seed]" << endl;
    return 1;
  }

  /********** Parse Arguments **********/
  MD::IO mdIO;
  ifstream finput;
  long seed;

  finput.open(argv[1]);
  if(mdIO.readInput(finput) != 0) {
    cout << "Error: fail to read input file: " << argv[1] << endl;
    return 2;
  }
  finput.close();
  mdIO.echoInput(cout);
  cout << "--------------------" << endl;

  seed = time(NULL);
  if(argc == 3) Util::readWord(argv[2], seed);
  cout <<  "SEED            " << seed << endl;
  Util::RNG::seed(seed);

  /********** Init MD **********/
  MD::Data mdData;

  if(mdIO.getData(mdData) != 0) {
    cout << "Error: fail to get MD::Data" << endl;
    return 3;
  }

  /********** Run MD **********/
  MD::Sim sim;
  int exitCode;

  if(mdIO.ensemble == "NVE") {
    MD::NVE mdEnsem;
    mdEnsem.init(mdData);
    cout << "Run MD::NVE" << endl;
    exitCode = sim.run(mdEnsem,mdIO);
    mdEnsem.final();
  }
  else if(mdIO.ensemble == "NVT") {
    MD::NVT mdEnsem;
    MD::NVT::Para mdPara;
    if(mdIO.getPara(mdPara) != 0) {
      cout << "Error: fail to read MD::NVT::Para" << endl;
      return 4;
    };
    mdEnsem.init(mdData,mdPara);
    cout << "Run MD::NVT: "
	 << "kT = " << mdPara.kT << "; tNH = " << mdPara.tNH << endl;
    exitCode = sim.run(mdEnsem,mdIO);
    mdEnsem.final();
  }
  else if(mdIO.ensemble == "NVE1He") {
    MD::NVE1He mdEnsem;
    MD::NVE1He::Para mdPara;
    if(mdIO.getPara(mdPara) != 0) {
      cout << "Error: fail to read MD::NVE1He::Para" << endl;
      return 5;
    }
    mdEnsem.init(mdData,mdPara);
    cout << "Run MD::NVE1He: ";
    cout << "posHe = (";
    for(int i = 0; i < 3; i++) cout << mdPara.posHe[i] << (i<2?",":"); ");
    cout << "velHe = ";							   
    for(int i = 0; i < 3; i++) cout << mdPara.velHe[i] << (i<2?",":")");
    exitCode = sim.run(mdEnsem,mdIO);
    mdEnsem.final();
  }
  else {
    cout << "Error: unrecognized ensemble: " << mdIO.ensemble << endl;
    return 6;
  }

  if(exitCode != 0) {
    cout << "Error: fail to run MD" << endl;
    return 7;
  }

  return 0;
}
