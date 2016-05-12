#include "../MC/MC.h"
#include <cmath>
#include <sstream>
using namespace std;

int main() {
  /********** PRE PARSE **********/
  long seed;
  string potTag;
  string agentTag;
  double pressure;
  double tension[9];
  double enthalpy;
  double volume0;
  double h0[9];
  double dxmax;
  double dhmax;
  long rStep;
  long eStep;
  long hStep;
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
      >> keyword >> pressure;
  cin >> keyword; for(int i = 0; i < 9; i++) cin >> tension[i];
  cin >> keyword >> enthalpy  
      >> keyword >> volume0;
  cin >> keyword; for(int i = 0; i < 9; i++) cin >> h0[i];

  pmc = MC::allocMC(potTag);
  if(pmc == NULL) {
    cout << "runMCNtH: unrecognized potential: " << potTag << endl;
    return 1;
  }  

  pmcAgent = MC::allocMCAgent(agentTag, cin);
  if(pmcAgent == NULL) {
    cout << "runMCNtH: unrecognized MCAgent: " << agentTag << endl;
    return 2;
  }

  cin >> keyword >> dxmax
      >> keyword >> dhmax
      >> keyword >> rStep
      >> keyword >> eStep
      >> keyword >> hStep
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
    cout << "runMCNtH: illegal parameters" << endl;
    return 3;
  }
  cin >> mol;
  if(cin.fail() || mol.natoms() == 0) {
    cout << "runMCNtH: no (initial) molecule is read" << endl;
    return 4;
  }
  if(mol.naxes() != 3) {
    cout << "runMCNtH: molecule must be in a 3D box" << endl;
    return 5;
  }

  /********** RUN **********/
  enthalpy *= mol.natoms();
  stream.setWorkDir(workdir);
  stream.init(digits,ios::app,false);
  MC::useNtH(seed,mol,*pmc,*pmcAgent,pressure,tension,enthalpy,
	     volume0,h0,dxmax,dhmax,rStep,eStep,hStep,stream);
  stream.final();

  return 0;
}
