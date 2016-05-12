#include "MolIO.h"
#include "Molecule.h"
#include "MolAnal.h"
#include <fstream>
#include <string>

namespace AMod {

/********** class MolI**********/
MolI::MolI() {}

MolI::MolI(Molecule& _mol): MolServ(_mol) {}

void MolI::reset() {}

int MolI::sourceTxt(std::istream& fin, Molecule& mol) {
  int k, tid;
  std::string tword;

  if(!fin.good()) return ISTREAM_ERROR;

  //read header   
  double ew[2];
  int ns[3], nt;
  fin >> tword >> tword; //extract "Energy xxx" and ignore
  for(k = 0; k < 2; k++) fin >> tword >> ew[k];
  for(k = 0; k < 3; k++) fin >> tword >> ns[k];
  if(!fin.good()) return ISTREAM_ERROR;
  if(ns[0] < 0 || ns[0] > 3) return BOX_ERROR;
  if(ns[1] < 0) return ATOM_ERROR;
  if(ns[2] < 0) return BOND_ERROR;
  mol.axes().resize(ns[0]);
  mol.atomsData().resize(ns[1]);
  mol.potential(ew[0]);
  mol.work(ew[1]);
  nt = ns[0]+ns[1];

  //read box and atoms
  Atom::Data* pdata;
  for(k = 0; k < nt; k++) {
    if(k < ns[0]) 
      pdata = &mol.axis(k);
    else 
      pdata = &mol.atomsData()[k-ns[0]];
    fin >> tid >> pdata->type 
	>> pdata->x[0] >> pdata->x[1] >> pdata->x[2] 
	>> pdata->fpot[0] >> pdata->fpot[1] >> pdata->fpot[2]
	>> pdata->fext[0] >> pdata->fext[1] >> pdata->fext[2]
	>> pdata->fixed[0] >> pdata->fixed[1] >> pdata->fixed[2];
    if(!fin.good()) return ISTREAM_ERROR;
  }
  mol.topo().sync();
  mol.topo().update(CASE_AXES);
  mol.topo().update(CASE_LONELY);

  //read bonds
  int id[3], n[3];
  for(k = 0; k < ns[2]; k++) {
    fin >> id[0] >> id[1] >> id[2] >> n[0] >> n[1] >> n[2];
    if(!fin.good()) return ISTREAM_ERROR;
    if(id[1] >= ns[1] || id[2] >= ns[1]) return BOND_ERROR;
    mol.topo().newBond(id[1],id[2],n);
  }

  //read(& skip) chains
  int nchains;
  std::streampos fpos = fin.tellg(); 
  fin >> tword;
  if(!fin.good()) {
    fin.clear();
    fin.seekg(fpos);
  }
  else if(tword != "Chains") {
    fin.seekg(fpos);
  }
  else {
    fin >> nchains;
    std::getline(fin, tword); //eat end of line
    for(k = 0; k < nchains; k++) std::getline(fin, tword);
  }

  return (fin.good() ? SUCCESS : ISTREAM_ERROR);
}

int MolI::sourceBin(std::istream& fin, Molecule& mol) {
  int k, tid;
  
  if(!fin.good()) return ISTREAM_ERROR;

  //read header
  double ew[2];
  int ns[3], nt;
  fin.read((char*)ew, sizeof(ew));
  fin.read((char*)ns, sizeof(ns));
  if(!fin.good()) return ISTREAM_ERROR;
  if(ns[0] < 0 || ns[0] > 3) return BOX_ERROR;
  if(ns[1] < 0) return ATOM_ERROR;
  if(ns[2] < 0) return BOND_ERROR;
  mol.axes().resize(ns[0]);
  mol.atomsData().resize(ns[1]);
  mol.potential(ew[0]);
  mol.work(ew[1]);
  nt = ns[0]+ns[1];

  //read lv and atoms
  Atom::Data* pdata;
  for(k = 0; k < nt; k++) {
    if(k < ns[0]) 
      pdata = &mol.axis(k);
    else 
      pdata = &mol.atomsData()[k-ns[0]];
    fin.read((char*)&tid, sizeof(tid));
    fin.read((char*)&pdata->type, sizeof(pdata->type));
    fin.read((char*)pdata->x, sizeof(pdata->x));
    fin.read((char*)pdata->fpot, sizeof(pdata->fpot));
    fin.read((char*)pdata->fext, sizeof(pdata->fext));
    fin.read((char*)pdata->fixed, sizeof(pdata->fixed));
  }
  mol.topo().sync();
  mol.topo().update(CASE_AXES);
  mol.topo().update(CASE_LONELY);

  //read bonds
  int id[3], n[3];
  for(k = 0; k < ns[2]; k++) {
    fin.read((char*)id, sizeof(id));
    fin.read((char*)n, sizeof(n));
    if(!fin.good()) return ISTREAM_ERROR;
    if(id[1] >= ns[1] || id[2] >= ns[1]) return BOND_ERROR;
    mol.topo().newBond(id[1],id[2],n);
  }

  return (fin.good() ? SUCCESS : ISTREAM_ERROR);
}

/********** class MolO **********/
MolO::MolO() { reset(); }

MolO::MolO(const Molecule& _mol): ConstMolServ(_mol) { reset(); }

void MolO::reset() { dumpChains(true); }

void MolO::setFormat(std::ostream& fout, int _precision) {
  fout.setf(std::ios::scientific, std::ios::floatfield);
  fout.precision(_precision);
  fout << std::left;
}

void MolO::dumpTxt(std::ostream& fout, const Molecule&mol, const MolO& molO) {
  dumpTxt(fout, mol, molO.dumpChains());
}

void MolO::dumpTxt(std::ostream& fout, const Molecule& mol, bool _dumpCh) {
  using std::endl;

  int k;
  int ns[3] = { mol.naxes(), mol.natoms(), mol.nbonds() };
  MolAnal::Chains chains;
  fout << "Energy    " << mol.totalEnergy() << endl;
  fout << "Potential " << mol.potential() << endl;
  fout << "Work      " << mol.work() << endl;  
  fout << "Box       " << ns[0] << endl;
  fout << "Atoms     " << ns[1] << endl;
  fout << "Bonds     " << ns[2] << endl;
  for(k = 0; k < ns[0]; k++) fout << mol.axis(k) << endl;
  for(k = 0; k < ns[1]; k++) fout << mol.atom(k) << endl;
  for(k = 0; k < ns[2]; k++) fout << mol.bond(k) << endl;
  if(_dumpCh) MolAnal().rings(mol, chains);
  fout << "Chains " << chains.size() << std::endl;
  fout << chains;
}

void MolO::dumpBin(std::ostream& fout, const Molecule& mol) {
  int k;
  double ew[2] = { mol.potential(), mol.work() };
  int ns[3] = { mol.naxes(), mol.natoms(), mol.nbonds() };
  fout.write((char*)ew, sizeof(ew));
  fout.write((char*)ns, sizeof(ns));
  for(k = 0; k < ns[0]; k++) mol.axis(k).dumpBin(fout);
  for(k = 0; k < ns[1]; k++) mol.atom(k).dumpBin(fout);
  for(k = 0; k < ns[2]; k++) mol.bond(k).dumpBin(fout);
}

/********** class MolIO **********/
MolIO::MolIO(): MolI(), MolO() {} 

MolIO::MolIO(Molecule& _mol): MolI(_mol), MolO(_mol) {}

void MolIO::attach(Molecule& _mol) {
  MolI::attach(_mol);
  MolO::attach(_mol);
}

void MolIO::reset() { MolI::reset(); MolO::reset(); }

} /* AMod */
