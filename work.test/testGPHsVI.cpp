#include "../AMod/AMod.h"
#include "../Util/util.h"
#include "../Util/RNG.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
  if(argc != 7) {
    cout << "Usage: testGPHsVI <in_mol> <out_mol> "
	 << "<layer_#1> <n_rm_atoms> <layer_#2> <n_ad_atoms>" << endl;
    cout << "       Remove <n_rm_atoms> in <layer_#1> and add <n_ad_atoms> UNDER <layer_#2>" << endl;
    return 1;
  }

  /********** PARSE ARGUMENT **********/
  AMod::Molecule mol;
  std::ifstream fin;
  std::ofstream fout;
  int layerNO1; int nRmAtoms;
  int layerNO2; int nAdAtoms;
  
  fin.open(argv[1]);
  if(!fin.good()) {
    cout << "Error: cannot open file: " << argv[1] << endl;
    return 2;
  }
  fin >> mol;
  fin.close();
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 3;
  }
  mol.setStdAxes();

  fout.open(argv[2]);
  if(!fout.good()) {
    cout << "Error: cannot open file: " << argv[2] << endl;
    return 4;
  }
  AMod::MolO::setFormat(fout,12);

  if(!Util::readWord(argv[3], layerNO1)) {
    cout << "Error: cannot read <layer_#1>" << endl;
    return 5;
  }

  if(!Util::readWord(argv[4], nRmAtoms)) {
    cout << "Error: cannot read <n_rm_atoms>" << endl;
    return 6;
  }

  if(!Util::readWord(argv[5], layerNO2)) {
    cout << "Error: cannot read <layer_#2>" << endl;
    return 7;
  }

  if(!Util::readWord(argv[6], nAdAtoms)) {
    cout << "Error: cannot read <n_ad_atoms>" << endl;
    return 8;
  }

  /********** ANALYZE LAYERS **********/
  double zmin, zmax, tmp;
  int natoms = mol.natoms();
  const double distGPH = 3.35;

  zmin = 0.0; zmax = 0.0;
  for(int i = 0; i < natoms; i++) {
    tmp = mol.atom(i)->x[2];
    if(i == 0 || zmin > tmp) zmin = tmp;
    if(i == 0 || zmax < tmp) zmax = tmp;
  }
  const int nlayers = int((zmax-zmin)/distGPH)+1;
  if(layerNO1 < 0 || layerNO1 >= nlayers ||
     layerNO2 < 0 || layerNO2 >= nlayers) {
    cout << "Error: <layer_#> should be in [0," << nlayers << ")" << endl;
    cout << "       <layer_#1> = " << layerNO1 << "; <layer_#2> = " << layerNO2 << endl;
    return 9;
  }

  std::vector< std::vector<int> > layerAtomIDs;
  layerAtomIDs.resize(nlayers);
  for(int i = 0; i < natoms; i++) {
    int layerNO = int(std::floor((zmax-mol.atom(i)->x[2])/distGPH+0.5));
    layerAtomIDs[layerNO].push_back(i);
  }

  /********** REMOVE ATOMS **********/
  Util::RNG::seed(time(NULL));
  int randIndex, atomID, fixed;
  std::vector<AMod::Atom*> rmAtomPtrs;
  rmAtomPtrs.resize(nRmAtoms);
  std::vector<int> oneLayer = layerAtomIDs[layerNO1];
  for(int i = 0; i < nRmAtoms; i++) {
    while(true) {
      randIndex = int(Util::RNG::uniform_SFMT()*oneLayer.size());
      atomID = oneLayer[randIndex];
      fixed = 0; for(int j = 0; j < 3; j++) fixed += mol.atom(atomID)->fixed[j];
      if(fixed == 0) { 
	oneLayer.erase(oneLayer.begin()+randIndex);
	break;
      }
    }
    rmAtomPtrs[i] = &mol.atom(atomID);
  }
  for(int i = 0; i < nRmAtoms; i++)
    mol.deleteAtom(rmAtomPtrs[i]->id());

  /********** ADD NEW ATOMS **********/
  const double* lA = mol.axis(0).x;
  const double* lB = mol.axis(1).x;
  double new_z = zmax-(layerNO2+0.5)*distGPH;
  double new_x, new_y;
  std::vector<int> newAtomIDs;
  for(int i = 0; i < nAdAtoms; i++) {
    double randA = Util::RNG::uniform_SFMT();
    double randB = Util::RNG::uniform_SFMT();
    new_x = randA*lA[0]+randB*lB[0];
    new_y = randA*lA[1]+randB*lB[1];
    AMod::Atom& a = mol.newAtom(-1,AMod::Molecule::AData(6,new_x,new_y,new_z));
    newAtomIDs.push_back(a.id());
  }

  /********** FIX ATOMS **********/
  AMod::Molecule molBak = mol;
  for(int i = 0; i < mol.naxes(); i++)
    for(int j = 0; j < 3; j++)
      mol.axis(i).fixed[j] = 1;
  for(int i = 0; i < natoms; i++)
    for(int j = 0; j < 3; j++)
      mol.atom(i)->fixed[j] = 1;
  for(int i = 0; i < int(newAtomIDs.size()); i++) {
    mol.atom(newAtomIDs[i])->fixed[0] = 0;
    mol.atom(newAtomIDs[i])->fixed[1] = 0;
  }

  /********** RELAX **********/
  AMod::MolOptim mopt(mol);
  mopt.init(AMod::LCBOPIIN);
  mopt.optimize(AMod::MolOptim::DonlpPara(argv[2],false,0.01));
  mopt.final();

  /********** RESTORE **********/
  for(int i = 0; i < mol.naxes(); i++)
    for(int j = 0; j < 3; j++)
      mol.axis(i).fixed[j] = molBak.axis(i).fixed[j];
  for(int i = 0; i < natoms; i++)
    for(int j = 0; j < 3; j++)
      mol.atom(i)->fixed[j] = molBak.atom(i)->fixed[j];

  /********** MOVE ATOMS INTO BOX **********/
  double cx[3];
  cx[0] = (lA[0]+lB[0])/2.0;
  cx[1] = (lA[1]+lB[1])/2.0;
  cx[2] = new_z;
  for(int i = 0; i < int(newAtomIDs.size()); i++)
    mol.topo().around(cx,newAtomIDs[i]);

  /********** OUTPUT MOL **********/
  mol.io().dumpTxt(fout);

  return 0;
}
