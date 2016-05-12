#include "../MC/MC.h"
#include <cmath>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {
  /********** PRE PARSE **********/
  long seed;
  string potTag;
  string agentTag;
  double kT;
  double dxmax;
  long rStep;
  long eStep;
  string workdir;
  AMod::Molecule mol;
  MC::MCBasic* pmc;
  MC::MCAgent* pmcAgent;
  MC::Stream stream;
  int digits;
  string keyword;

  cin >> keyword >> seed
      >> keyword >> potTag
      >> keyword >> agentTag
      >> keyword >> kT;

  pmc = MC::allocMC(potTag);
  if(pmc == NULL) {
    cout << "runMCNVT: unrecognized potential: " << potTag << endl;
    return 1;
  }  

  pmcAgent = MC::allocMCAgent(agentTag, cin);
  if(pmcAgent == NULL) {
    cout << "runMCNVT: unrecognized MCAgent: " << agentTag << endl;
    return 2;
  }

  cin >> keyword >> dxmax
      >> keyword >> rStep
      >> keyword >> eStep
      >> keyword >> stream.showHeader()
      >> keyword >> stream.dAction() 
      >> keyword >> stream.dStateStep() 
      >> keyword >> stream.dStatStep() 
      >> keyword >> stream.dInfoStep()
      >> keyword >> workdir
      >> keyword >> stream.actionFile()
      >> keyword >> stream.stateFile()
      >> keyword >> stream.statFile()
      >> keyword >> stream.infoFile()
      >> keyword >> digits;

  if(cin.fail()) {
    cout << "runMCNVT: illegal parameters" << endl;
    return 3;
  }
  cin >> mol;
  if(cin.fail() || mol.natoms() == 0) {
    cout << "runMCNVT: no (initial) molecule is read" << endl;
    return 4;
  }
  if(mol.naxes() != 3) {
    cout << "runMCNVT: molecule must be in a 3D box" << endl;
    return 5;
  }

  /********** RUN **********/
  stream.setWorkDir(workdir);
  stream.init(digits,ios::app,false);
  MC::useNVT(seed,mol,*pmc,*pmcAgent,kT,dxmax,rStep,eStep,stream);
  stream.final();

  return 0;
}
