#include "../AMod/AMod.h"
#include "../Util/util.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

int main(int argc, char* argv[]) {
  using std::cout;
  using std::endl;

  if(argc != 5) {
    cout << "Usage: testGPHonHole <in_mol> <out_mol> <radius> <depth>" << endl;
    cout << "       input mol must be graphene layers" << endl;
    cout << "       <depth> indicates how deep the probe tip goes" << endl;
    return 1;
  }

  /********** PARSE ARGUMENT **********/
  AMod::Molecule mol;
  std::ifstream fin;
  std::ofstream fout;
  double radius;
  double depth;
  
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

  if(!Util::readWord(argv[3], radius)) {
    cout << "Error: cannot read <radius>" << endl;
    return 5;
  }

  if(!Util::readWord(argv[4], depth)) {
    cout << "Error: cannot read <depth>" << endl;
    return 6;
  }

  /********** FIX ATOMS OUTSIDE R **********/
  int i, j;
  int naxes = mol.naxes();
  int natoms = mol.natoms();
  //fix all axes
  for(i = 0; i < naxes; i++)
    for(j = 0; j < 3; j++) mol.axis(i).fixed[j] = 1;
  //unfix all atoms
  for(i = 0; i < natoms; i++)
    for(j = 0; j < 3; j++) mol.atom(i)->fixed[j] = 0;

  double zmin = 0.0, zmax = 0.0, tmp;
  const double distGPH = 3.35;
  //find zmin & zmax;
  for(i = 0; i < natoms; i++) {
    tmp = mol.atom(i)->x[2];
    if(i == 0 || zmin > tmp) zmin = tmp;
    if(i == 0 || zmax < tmp) zmax = tmp;
  }
  const int nlayers = int((zmax-zmin)/distGPH)+1;

  int topAtomID = -1; double topR = Util::MAX_DOUBLE;
  int botAtomID = -1; double botR = Util::MAX_DOUBLE;
  AMod::MolAdjuster::Point mc = AMod::MolAdjuster::massCenter(mol);
  std::vector<int> atomLayerNO(natoms,0);
  //find top atom & bottom atom that are closest to mc
  for(i = 0; i < natoms; i++) {
    tmp = std::sqrt(Util::vvdist2(mol.atom(i)->x,mc.data,2));
    if(tmp < topR && std::abs(mol.atom(i)->x[2]-zmax) < 1.0) {
      topAtomID = i; topR = tmp;
    }
    if(tmp < botR && std::abs(mol.atom(i)->x[2]-zmin) < 1.0) {
      botAtomID = i; botR = tmp;
    }
    atomLayerNO[i] = int(std::floor((zmax-mol.atom(i)->x[2])/distGPH+0.5));
  }

  const double minDist = 0.5*distGPH; //may be adjusted
  const double* xc = mol.atom(topAtomID)->x;
  std::vector<double> depths(nlayers,depth);
  double adisp;
  int count = 0;
  //distance to go for each layer
  for(i = 1; i < nlayers; i++)
    depths[i] = Util::max(depths[i-1]-(distGPH-minDist),0.0);
  //fix the top atom
  mol.atom(topAtomID)->fixed[2] = 1;
  //fix atoms outside topAtom+radius
  for(i = 0; i < natoms; i++) {
    tmp = std::sqrt(Util::vvdist2(mol.atom(i)->x,xc,2));
    if(tmp > radius) 
      for(j = 0; j < 3; j++) mol.atom(i)->fixed[j] = 1;
    else {//press down
      adisp = (1.0-tmp/radius)*depths[atomLayerNO[i]];
      mol.atom(i)->x[2] -= adisp;
      count++;
    }
  }
  //update topology
  mol.topo().update(AMod::CASE_GENERAL);

  //output info
  cout << "# of layers: " << nlayers << endl;
  cout << "# of atoms : " << count 
       << " inside the hole (r = " << radius << " A)" << endl;
  cout << "Center atom: atomID = " << topAtomID << ", "
       << "x = (" << xc[0] << "," << xc[1] << "," << xc[2] << ")" << endl;

  /********** MAKE OUTPUT GPH **********/
  mol.io().dumpTxt(fout);

  return 0;
}
