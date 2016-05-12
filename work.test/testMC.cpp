#include "../MC/MC.h"
#include "../MC/MC_REBO.h"
#include "../MC/MC_LCBOPII.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
using namespace std;

const double THRESHOLD = 1.0e-5;

template<class T>
ostream& dumpV(const T& t, int n, ostream& fout);

void dumpArrow(const AMod::Arrow& ar, ostream& fout);

void dumpBondDiff_REBO(const REBO::BondDatum& bd,
		       const REBO::BondDatum& bdREF,
		       ostream& fout);

void dumpDiff_REBO(const MC::MC_REBO& mc, 
		   const MC::MC_REBO& mcREF, 
		   ostream& fout);

void dumpDiff_LCBOPII(const MC::MC_LCBOPII& mc, 
		      const MC::MC_LCBOPII& mcREF, 
		      ostream& fout);

int main(int argc, char* argv[]) {
  if(argc != 5 && argc != 6) {
    cout << "Usage: testMC <file_mol> <potential> <nStep> <file_debug> <seed>" 
	 << endl;
    return 1;
  }

  
  /********** READ ARGUMENTS **********/
  string fmol;
  string pot;
  string fdebug;
  long nStep;
  long seed;
  istringstream iss;

  fmol = argv[1];
  pot = argv[2];
  iss.clear();
  iss.str(argv[3]);
  iss >> nStep;
  if(iss.fail()) {
    cout << "Error: illegal nStep" << endl;
    return 2;
  }
  fdebug = argv[4];
  if(argc < 6) {
    seed = time(NULL);
  }
  else {
    iss.clear();
    iss.str(argv[5]);
    iss >> seed;
    if(iss.fail()) {
      cout << "Error: illegal seed" << endl;
      return 3;
    }
  }
  
  /********** PROCESS ARGUMENTS **********/
  AMod::Molecule mol, molREF;
  MC::MCBasic *pmc, *pmcREF;
  MC::MC_REBO *pmc_REBO = NULL, *pmcREF_REBO = NULL;
  MC::MC_LCBOPII *pmc_LCBOPII = NULL, *pmcREF_LCBOPII = NULL;
  double dxmax, dhmax;
  long natoms;
  ifstream fin;

  if(pot == "REBO") {
    pmc_REBO = new MC::MC_REBO;
    pmcREF_REBO = new MC::MC_REBO;
    pmc = pmc_REBO;
    pmcREF = pmcREF_REBO;
  }
  else if(pot == "LCBOPII") {
    pmc_LCBOPII = new MC::MC_LCBOPII;
    pmcREF_LCBOPII = new MC::MC_LCBOPII;
    pmc = pmc_LCBOPII;
    pmcREF = pmcREF_LCBOPII;
  }
  else {
    cout << "Error: illegal potential" << endl;;
    return 4;
  } 
  fin.open(fmol.c_str());
  if(!fin.good()) {
    cout << "Error: cannot open file: " << fmol << endl;
    return 5;
  }
  fin >> mol;
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 6;
  }
  fin.close();
  dxmax = 2.0;
  dhmax = 1.0;
  natoms = mol.natoms();

  /********** SETUP (PSEUDO MC) **********/
  Util::RNG::seed(seed);
  cout << "seed = " << seed << endl;
  cout << "******************************" << endl;
  mol.topo().mode(AMod::MODE_CELL);
  mol.topo().update(AMod::CASE_CHANGE_MODE);

  /********** CHECK **********/
  if(!AMod::MolChecker(mol).check()) return 7;
  cout << "******************************" << endl;

  pmc->attach(mol);
  pmc->init();
  
  /********** RUN **********/
  long kStep;
  long kkStep;
  long atomID;
  long nT[3];
  int moveType;
  double old_x[9];
  double new_x[9];
  double old_h[9];
  double new_h[9];
  double inv_old_h[9];
  double inv_new_h[9];
  double d_h[9];
  double s_x[3];
  double dE;
  double energyi;
  double energyf;
  double dEREF;
  double energyREFi;
  double energyREFf;
  double volume;
  bool accept;
  double xErrorMax;
  double hErrorMax;
  ofstream dfout;
  bool debugStep;
  
  xErrorMax = hErrorMax = 0.0;
  for(kStep = 0; kStep < nStep; kStep++) {
    for(kkStep = 0; kkStep < natoms; kkStep++) {
      debugStep = false; //(kStep == 798 && kkStep == 11);
      molREF = mol;
      moveType = int(Util::RNG::uniform_SFMT()*2.0);
      atomID = Util::RNG::uniform_SFMT()*natoms;
      for(int k = 0; k < 3; k++) {
	old_x[k] = mol.atom(atomID)->x[k];
	new_x[k] = old_x[k];
      }
      for(int k = 0; k < 3; k++) nT[k] = 0;
      if(moveType == 0) {
	for(int k = 0; k < 3; k++) {
	  if(!mol.atom(atomID)->fixed[k]) {
	    new_x[k] += dxmax*(2.0*Util::RNG::uniform_SFMT()-1);
	  }
	}
      }//end if(moveType == 0)
      else if(moveType == 1) {
	for(int j = 0; j < mol.naxes(); j++) 
	  nT[j] = int(Util::RNG::uniform_SFMT()*10.0);
	for(int k = 0; k < 3; k++)
	  for(int j = 0; j < mol.naxes(); j++)
	    new_x[k] += nT[j]*mol.axis(j).x[k];
      }//end else if(moveType == 1)
      //--------------------------------------------------
      molREF.topo().update(AMod::CASE_AXES);
      molREF.topo().update(AMod::CASE_GENERAL);
      pmcREF->attach(molREF);
      pmcREF->init();
      energyREFi = pmcREF->computeEnergy();
      pmcREF->final();
      if(debugStep) {
	dfout.open("debug.MC.I");
	int bID = mol.topo().bondID(25,38);
	int bIDREF = molREF.topo().bondID(25,38);
	dfout << "bID = " << bID << endl;
	dfout << "bIDREF = " << bIDREF << endl;
	REBO::BondDatum& bd = pmc_REBO->energy().bondData[bID];
	REBO::BondDatum& bdREF = pmcREF_REBO->energy().bondData[bIDREF];
	dumpBondDiff_REBO(bd,bdREF,dfout);
	dfout.close();
      }
      //**------------------------------------------------
      if(debugStep) {
	dfout.open("mol.I.log");
	AMod::MolO::setFormat(dfout);
	mol.io().dumpTxt(dfout);
	dfout.close();
      }
      if(debugStep) {
      }
      dE = pmc->dEnergyMoveAtom(atomID,new_x);
      if(debugStep) {
      }
      if(debugStep) {
	dfout.open("mol.F.log");
	AMod::MolO::setFormat(dfout);
	mol.io().dumpTxt(dfout);
	dfout.close();
      }
      //--------------------------------------------------
      if(debugStep) {
      }
      molREF = mol;
      molREF.topo().update(AMod::CASE_AXES);
      molREF.topo().update(AMod::CASE_GENERAL);
      pmcREF->attach(molREF);
      pmcREF->init();
      energyREFf = pmcREF->computeEnergy();
      pmcREF->final();
      dEREF = energyREFf-energyREFi;
      if(debugStep) {
      }
      if(debugStep) {
	dfout.open("debug.MC.F");
	int bID = mol.topo().bondID(25,38);
	int bIDREF = molREF.topo().bondID(25,38);
	dfout << "bID = " << bID << endl;
	dfout << "bIDREF = " << bIDREF << endl;
	REBO::BondDatum& bd = pmc_REBO->energy().bondData[bID];
	REBO::BondDatum& bdREF = pmcREF_REBO->energy().bondData[bIDREF];
	dumpBondDiff_REBO(bd,bdREF,dfout);
	dfout.close();
      }
      //--------------------------------------------------
      //dump
      cout << kStep << " " << kkStep << " ";
      cout << "[MoveAtom--" << atomID << "] "
	   << "[MoveType--" << moveType << "~("
	   << nT[0] << "," << nT[1] << "," << nT[2] << ")] "
	   << "|dE-dEREF| = " << fabs(dE-dEREF)
	   << " = |" << dE << " - " << dEREF << "|" << endl;
      if(fabs(dE-dEREF) > xErrorMax) xErrorMax = fabs(dE-dEREF);
      //compare
      if(fabs(dE-dEREF) > THRESHOLD) {//something may wrong
	cout << "WARNING!" << endl;
	cout << "atomID = " << atomID << endl;
	cout << "atom.old.x = (" 
	     << old_x[0] << "," << old_x[1] << "," << old_x[2] << ")" << endl;
	cout << "atom.new.x = ("
	     << new_x[0] << "," << new_x[1] << "," << new_x[2] << ")" << endl;
	cout << "natoms = " << mol.natoms() << " ~ " << molREF.natoms() << endl;
	cout << "nbonds = " << mol.nbonds() << " ~ " << molREF.nbonds() << endl;
	cout << "dE = " << dE << endl;
	cout << "dEREF = " << energyREFf << " - " << energyREFi << endl;
	dfout.open(fdebug.c_str());
	if(pot == "REBO")
	  dumpDiff_REBO(*pmc_REBO, *pmcREF_REBO, dfout);
	else if(pot == "LCBOPII")
	  dumpDiff_LCBOPII(*pmc_LCBOPII, *pmcREF_LCBOPII, dfout);
	dfout.close();
	return 8;
      }//end if(fabs(dE-dEREF)...
      //compare
      if(moveType == 1 && fabs(dE) > THRESHOLD) {
	cout << "WARNING!" << endl;
	return 9;
      }
      //pseudo MC
      accept = (Util::RNG::uniform_SFMT() < 0.5);
      if(accept) {//accept
	cout << "accept " << endl;
	cout << "---------------------------------------------------------"  
	     << endl;
      }
      else {//reject
	pmc->unMoveAtom();
	cout << "reject" << endl;
	cout << "---------------------------------------------------------"  
	     << endl;
      }
    }//end for(kkStep...
    //change volume
    molREF = mol;
    pmc->save();
    for(int k = 0; k < 3; k++) Util::vcopy(mol.axis(k).x, old_h+3*k, 3);
    Util::minverse(old_h, inv_old_h, 3);
    //dump
    cout << "[ChangeVolume]" << endl;
    cout << "old.h = ";
    dumpV(old_h, 9, cout) << endl;
    do {
      Util::vassign(0.0,d_h,9);
      for(int i = 0; i < 3; i++)
	for(int j = i; j < 3; j++) {
	  d_h[3*i+j] = dhmax*(2.0*Util::RNG::uniform_SFMT()-1);
	  if(i != j) d_h[3*j+i] = d_h[3*i+j];
	}
      Util::madd(old_h, d_h, new_h, 3, 3);
      volume = fabs(Util::cross_dot3(new_h,new_h+3,new_h+6));
      cout << "new.h = ";
      dumpV(new_h, 9, cout) << endl;
      cout << "volume = " << volume << endl;
    } while(volume < 10.0);
    Util::minverse(new_h, inv_new_h, 3);
    for(int k = 0; k < 3; k++) Util::vcopy(new_h+3*k, mol.axis(k).x, 3);
    for(int i = 0; i < natoms; i++) {
      Util::mdotv(inv_old_h, mol.atom(i)->x, s_x, 3, 3);
      Util::mdotv(new_h, s_x, mol.atom(i)->x, 3, 3);
    }
    mol.topo().update(AMod::CASE_AXES);
    mol.topo().update(AMod::CASE_GENERAL);
    energyi = mol.potential();
    energyf = pmc->computeEnergy();
    dE = energyf-energyi;
    //--------------------------------------------------
    molREF.topo().update(AMod::CASE_AXES);
    molREF.topo().update(AMod::CASE_GENERAL);
    pmcREF->attach(molREF);
    pmcREF->init();
    energyREFi = pmcREF->computeEnergy();
    pmcREF->final();
    //--------------------------------------------------
    molREF = mol;
    molREF.topo().update(AMod::CASE_AXES);
    molREF.topo().update(AMod::CASE_GENERAL);
    pmcREF->attach(molREF);
    pmcREF->init();
    energyREFf = pmcREF->computeEnergy();
    pmcREF->final();
    dEREF = energyREFf-energyREFi;
    //--------------------------------------------------
    //dump
    cout << "|dE-dEREF| = " << fabs(dE-dEREF)
	 << " = |" << dE << " - " << dEREF << "|" << endl;
    cout << "Ei = " << energyi << " ~ " << energyREFi << endl;
    cout << "Ef = " << energyf << " ~ " << energyREFf << endl;
    if(fabs(dE-dEREF) > xErrorMax) hErrorMax = fabs(dE-dEREF);
    //compare
    if(fabs(dE-dEREF) > THRESHOLD) {//something may wrong
      cout << "WARNING!" << endl;
      cout << "old.h = (" 
	   << old_h[0] << "," << old_h[1] << "," << old_h[2] << ")" << endl;
      cout << "new.h = ("
	   << new_h[0] << "," << new_h[1] << "," << new_h[2] << ")" << endl;
      cout << "natoms = " << mol.natoms() << " ~ " << molREF.natoms() << endl;
      cout << "nbonds = " << mol.nbonds() << " ~ " << molREF.nbonds() << endl;
      return 10;
    }//end if(fabs(dE-dEREF)... 
    //pseudo MC
    accept = (Util::RNG::uniform_SFMT() < 0.5);
    if(accept) {//accept
    }
    else {//reject
      pmc->restore();
    }
  }//end for(kStep...
  cout << "[Summary]" << endl;
  cout << "dEErrorMax (MoveAtom)     = " << xErrorMax << endl;
  cout << "dEErrorMax (ChangeVolume) = " << hErrorMax << endl;

  /********** FINALIZE **********/
  pmc->final();

  return 0;
}

template<class T>
ostream& dumpV(const T& t, int n, ostream& fout) {
  fout << "(";
  for(int k = 0; k < n; k++) {
    fout << t[k];
    if(k < n-1) fout << ",";
    else fout << ")";
  }
  return fout;
}

void dumpArrow(const AMod::Arrow& ar, ostream& fout) {
  fout << "atomID(host) = " << ar.host().id() << endl;
  fout << "atomID(moon) = " << ar.moon().id() << endl;
  fout << "atom.x(host) = "; dumpV(ar.host()->x,3,fout) << endl;
  fout << "atom.x(moon) = "; dumpV(ar.moon()->x,3,fout) << endl;
  fout << "bond.r = " << ar.bond().r << endl;
  fout << "arrow.dx = "; dumpV(ar.dx,3,fout) << endl;
  fout << "arrow.n = "; dumpV(ar.n,3,fout) << endl;
}

void dumpBondDiff_REBO(const REBO::BondDatum& bd,
		       const REBO::BondDatum& bdREF,
		       ostream& fout) {
  fout << "type = ";
  dumpV(bd.type,2,fout) << " ~ ";
  dumpV(bdREF.type,2,fout) << endl;
  fout << "bsigmapi = ";
  dumpV(bd.bsigmapi,2,fout) << " ~ ";
  dumpV(bdREF.bsigmapi,2,fout) << endl;
  fout << "N[0] = ";
  dumpV(bd.N[0],3,fout) << " ~ ";
  dumpV(bdREF.N[0],3,fout) << endl;
  fout << "N[1] = ";
  dumpV(bd.N[1],3,fout) << " ~ ";
  dumpV(bdREF.N[1],3,fout) << endl;
  fout << "F = ";
  dumpV(bd.F,3,fout) << " ~ ";
  dumpV(bdREF.F,3,fout) << endl;
  fout << "Q = " << bd.Q << " ~ " << bdREF.Q << endl;
  fout << "A = " << bd.A << " ~ " << bdREF.A << endl;
  fout << "alpha = " << bd.alpha << " ~ " << bdREF.alpha << endl;
  fout << "B = ";
  dumpV(bd.B,3,fout) << " ~ ";
  dumpV(bdREF.B,3,fout) << endl;
  fout << "beta = ";
  dumpV(bd.beta,3,fout) << " ~ ";
  dumpV(bdREF.beta,3,fout) << endl;
  fout << "sum_theta2 = " 
       << bd.sum_theta2 << " ~ " << bdREF.sum_theta2 << endl;
  fout << "bDH = " << bd.bDH << " ~ " << bdREF.bDH << endl;
  fout << "bbar = " << bd.bbar << " ~ " << bdREF.bbar << endl;
  fout << "fc = " << bd.fc << " ~ " << bdREF.fc << endl;
  fout << "sumFfc = ";
  dumpV(bd.sumFfc,2,fout) << " ~ ";
  dumpV(bdREF.sumFfc,2,fout) << endl;
  fout << "Nconj = " << bd.Nconj << " ~ " << bdREF.Nconj << endl;
  fout << "Fij = " << bd.Fij << " ~ " << bdREF.Fij << endl;
  fout << "Tij = " << bd.Tij << " ~ " << bdREF.Tij << endl;
  fout << "VApre = " << bd.VApre << " ~ " << bdREF.VApre << endl;
  fout << "VR = " << bd.VR << " ~ " << bdREF.VR << endl;
  fout << "VA = " << bd.VA << " ~ " << bdREF.VA << endl;
  fout << "E = " << bd.E << " ~ " << bdREF.E << endl;
}

void dumpDiff_REBO(const MC::MC_REBO& mc, 
		   const MC::MC_REBO& mcREF, 
		   ostream& fout) {
  double maxdE;
  const AMod::Molecule& mol = mc.molecule();
  const AMod::Molecule& molREF = mcREF.molecule();
  fout << "############### NAB ###############" << endl;
  if(mol.natoms() != molREF.natoms()) {//error
    fout << "Error[mol.natoms]:" << endl 
	 << mol.natoms() << " != " 
	 << molREF.natoms() << endl;
    fout << "--------------------" << endl;
    return;
  }
  if(mol.nbonds() != molREF.nbonds()) {//error
    fout << "Error[mol.nbonds]:" << endl 
	 << mol.nbonds() << " != " 
	 << molREF.nbonds() << endl;
    fout << "--------------------" << endl;
    return;
  }
  fout << "############### BOND ###############" << endl;
  maxdE = 0.0;
  for(int k = 0, kmax = mol.natoms(); k < kmax; k++) {
    const AMod::Atom& a = mol.atom(k);
    const AMod::Atom& aREF = molREF.atom(k);
    int jmax = a.arrows.size();
    int jmaxREF = aREF.arrows.size();
    if(jmax != jmaxREF) {//error
      fout << "Error[arrows.size]:" << endl
	   << jmax << " != " << jmaxREF << endl;
      fout << "atomID = " << k << endl
	   << "atom.x = "; dumpV(a->x,3,fout) << endl;
      fout << "--------------------" << endl;
      return; 
    }
    for(int j = 0; j < jmax; j++) {
      const AMod::Atom& ap = a.arrows[j].moon();
      int kp = ap.id();
      if(k > kp) continue;
      int bID = a.arrows[j].bond().id();
      int bIDREF = molREF.topo().bondID(k,kp);
      if(bIDREF < 0) {//error
	fout << "Error[neighbor]:" << endl
	     << "atomID(host) = " << k << " is not bonded to "
	     << "atomID(moon) = " << kp << endl;
	fout << "atom.x(host) = "; dumpV(a->x,3,fout) << endl;
	fout << "atom.x(moon) = "; dumpV(ap->x,3,fout) << endl;
	fout << "--------------------" << endl;
	return;
      }
      int arID = j;
      int arIDREF = molREF.atom(k).arrowID(kp);
      const AMod::Arrow& ar = mol.atom(k).arrows[arID];
      const AMod::Arrow& arREF = molREF.atom(k).arrows[arIDREF];
      double ddx[3];
      int dn[3];
      Util::vsub(ar.dx, arREF.dx, ddx, 3);
      Util::vsub(ar.n, arREF.n, dn, 3);
      if(Util::norm2(ddx,3) > THRESHOLD ||
	 std::sqrt(double(Util::vdotv(dn,dn,3))) > THRESHOLD) {//error
	fout << "Error[ar]:" << endl;
	fout << "[Arrow: ME]" << endl;
	dumpArrow(ar, fout);
	for(int c = 0; c < mol.naxes(); c++) {
	  fout << "T_" << c << " = "; dumpV(mol.axis(c),3,fout) << endl; 
	}
	fout << "[Arrow: REF]" << endl;
	dumpArrow(arREF, fout);
	for(int c = 0; c < molREF.naxes(); c++) {
	  fout << "T_" << c << " = "; dumpV(molREF.axis(c),3,fout) << endl; 
	}
	fout << "--------------------" << endl;
      }
      REBO::BondDatum& bd = mc.energy().bondData[bID];
      REBO::BondDatum& bdREF = mcREF.energy().bondData[bIDREF];
      if(fabs(bd.E-bdREF.E) > THRESHOLD/mc.energy().bondData.size()) {//error
	fout << "Error[bond.E]:" << endl
	     << "|" << bd.E << "-" << bdREF.E << "| = " 
	     << fabs(bd.E-bdREF.E) << endl;
	fout << "[Arrow]" << endl;
	dumpArrow(ar, fout);
	for(int c = 0; c < mol.naxes(); c++) {
	  fout << "T_" << c << " = "; dumpV(mol.axis(c),3,fout) << endl; 
	}
	dumpBondDiff_REBO(bd,bdREF,fout);
	fout << "--------------------" << endl;
	return;
      }
      else if(fabs(bd.E-bdREF.E) > maxdE) maxdE = fabs(bd.E-bdREF.E);
    }//end for(int j...
  }//end for(int k...
}

void dumpDiff_LCBOPII(const MC::MC_LCBOPII& mc, 
		      const MC::MC_LCBOPII& mcREF, 
		      ostream& fout) {
  double maxdE;
  const AMod::Molecule& mol = mc.molecule();
  const AMod::Molecule& molREF = mcREF.molecule();
  fout << "############### NAB ###############" << endl;
  if(mol.natoms() != molREF.natoms()) {//error
    fout << "Error[mol.natoms]:" << endl 
	 << mol.natoms() << " != " 
	 << molREF.natoms() << endl;
    fout << "--------------------" << endl;
    return;
  }
  if(mol.nbonds() != molREF.nbonds()) {//error
    fout << "Error[mol.nbonds]:" << endl 
	 << mol.nbonds() << " != " 
	 << molREF.nbonds() << endl;
    fout << "--------------------" << endl;
    return;
  }
  fout << "############### SR ###############" << endl;
  maxdE = 0.0;
  for(int k = 0, kmax = mol.natoms(); k < kmax; k++) {
    const AMod::Atom& a = mol.atom(k);
    const AMod::Atom& aREF = molREF.atom(k);
    int jmax= a.arrows.size();
    int jmaxREF = aREF.arrows.size();
    if(jmax != jmaxREF) {//error
      fout << "Error[arrows.size]:" << endl
	   << jmax << " != " << jmaxREF << endl;
      fout << "atomSRID = " << k << endl
	   << "atomSR.x = "; dumpV(a->x,3,fout) << endl;
      fout << "--------------------" << endl;
      return;
    }
    for(int j = 0; j < jmax; j++) {
      const AMod::Atom& ap = a.arrows[j].moon();
      int kp = ap.id();
      if(k > kp) continue;
      int bID = a.arrows[j].bond().id();
      int bIDREF = molREF.topo().bondID(k,kp);
      if(bIDREF < 0) {//error
	fout << "Error[neighbor]:" << endl
	     << "atomSRID(host) = " << k << " is not bonded to "
	     << "atomSRID(moon) = " << kp << endl;
	fout << "atomSR.x(host) = "; dumpV(a->x,3,fout) << endl;
	fout << "atomSR.x(moon) = "; dumpV(ap->x,3,fout) << endl;
	fout << "--------------------" << endl;
	return;
      }
      LCBOPII::BondDatumSR& bd = mc.energy().bondDataSR[bID];
      LCBOPII::BondDatumSR& bdREF = mcREF.energy().bondDataSR[bIDREF];
      if(fabs(bd.E-bdREF.E) > THRESHOLD/mc.energy().bondDataSR.size()) {//error
	fout << "Error[bondSR.E]:" << endl
	     << "|" << bd.E << "-" << bdREF.E << "| = " 
	     << fabs(bd.E-bdREF.E) << endl;
	fout << "atomSRID(host) = " << k << endl
	     << "atomSRID(moon) = " << kp << endl;
	fout << "atomSR.x(host) = "; dumpV(a->x,3,fout) << endl;
	fout << "atomSR.x(moon) = "; dumpV(ap->x,3,fout) << endl;
	fout << "SdN = " << bd.SdN << " ~ " << bdREF.SdN << endl;
	fout << "N_ = "; 
	dumpV(bd.N_,2,fout) << " ~ "; 
	dumpV(bdREF.N_,2,fout) << endl;
	fout << "M_ = "; 
	dumpV(bd.M_,2,fout) << " ~ "; 
	dumpV(bdREF.M_,2,fout) << endl;
	fout << "Suel_ = "; 
	dumpV(bd.Suel_,2,fout) << " ~ "; 
	dumpV(bdREF.Suel_,2,fout) << endl;
	fout << "E = " << bd.E << " ~ " << bdREF.E << endl;
	fout << "--------------------" << endl;
	return;
      }
      else if(fabs(bd.E-bdREF.E) > maxdE) maxdE = fabs(bd.E-bdREF.E);
    }//end for(int j...
  }//end for(int k...
  fout << "maxdE = " << maxdE << endl;
  fout << "############### LR ###############" << endl;  
  maxdE = 0.0;
  for(int k = 0, kmax = mol.natoms(); k < kmax; k++) {
    const AMod::Atom& a = mc.molTopoLR().atom(k);
    const AMod::Atom& aREF = mcREF.molTopoLR().atom(k);
    int jmax= a.arrows.size();
    int jmaxREF = aREF.arrows.size();
    if(jmax != jmaxREF) {//error
      fout << "Error[arrows.size]:" << endl
	   << jmax << " != " << jmaxREF << endl;
      fout << "atomLRID = " << k << endl
	   << "atomLR.x = "; dumpV(a->x,3,fout) << endl;
      fout << "--------------------" << endl;
      return;
    }
    for(int j = 0; j < jmax; j++) {
      const AMod::Atom& ap = a.arrows[j].moon();
      int kp = ap.id();
      if(k > kp) continue;
      int bID = a.arrows[j].bond().id();
      int bIDREF = mcREF.molTopoLR().bondID(k,kp);
      if(bIDREF < 0) {//error
	fout << "Error[neighbor]:" << endl
	     << "atomLRID(host) = " << k << " is not bonded to "
	     << "atomLRID(moon) = " << kp << endl;
	fout << "atomLR.x(host) = "; dumpV(a->x,3,fout) << endl;
	fout << "atomLR.x(moon) = "; dumpV(ap->x,3,fout) << endl;
	fout << "--------------------" << endl;
	return;
      }
      LCBOPII::BondDatumLR& bd = mc.energy().bondDataLR[bID];
      LCBOPII::BondDatumLR& bdREF = mcREF.energy().bondDataLR[bIDREF];
      if(fabs(bd.E-bdREF.E) > THRESHOLD/mc.energy().bondDataLR.size()) {//error
	fout << "Error[bondLR.E]:" << endl
	     << "|" << bd.E << "-" << bdREF.E << "| = " 
	     << fabs(bd.E-bdREF.E) << endl;
	fout << "atomLRID(host) = " << k << endl
	     << "atomLRID(moon) = " << kp << endl;
	fout << "atomLR.x(host) = "; dumpV(a->x,3,fout) << endl;
	fout << "atomLR.x(moon) = "; dumpV(ap->x,3,fout) << endl;
	fout << "E = " << bd.E << " ~ " << bdREF.E << endl;
	fout << "--------------------" << endl;
	return;
      }
      else if(fabs(bd.E-bdREF.E) > maxdE) maxdE = fabs(bd.E-bdREF.E);
    }//end for(int j...
  }//end for(int k...
  fout << "maxdE = " << maxdE << endl;
  fout << "############### MR ###############" << endl;
  maxdE = 0.0;
  for(int k = 0, kmax = mol.natoms(); k < kmax; k++) {
    const AMod::Atom& a = mc.molTopoMR().atom(k);
    const AMod::Atom& aREF = mcREF.molTopoMR().atom(k);
    int jmax= a.arrows.size();
    int jmaxREF = aREF.arrows.size();
    if(jmax != jmaxREF) {//error
      fout << "Error[arrows.size]:" << endl
	   << jmax << " != " << jmaxREF << endl;
      fout << "atomMRID = " << k << endl
	   << "atomMR.x = "; dumpV(a->x,3,fout) << endl;
      fout << "--------------------" << endl;
      return;
    }
    for(int j = 0; j < jmax; j++) {
      const AMod::Atom& ap = a.arrows[j].moon();
      int kp = ap.id();
      if(k > kp) continue;
      int bIDREF = mcREF.molTopoMR().bondID(k,kp);
      if(bIDREF < 0) {//error
	fout << "Error[neighbor]:" << endl
	     << "atomMRID(host) = " << k << " is not bonded to "
	     << "atomMRID(moon) = " << kp << endl;
	fout << "atomMR.x(host) = "; dumpV(a->x,3,fout) << endl;
	fout << "atomMR.x(moon) = "; dumpV(ap->x,3,fout) << endl;
	fout << "--------------------" << endl;
	return;
      }
      LCBOPII::AtomDatumMR& ad = mc.energy().atomDataMR[k];
      LCBOPII::AtomDatumMR& adREF = mcREF.energy().atomDataMR[k];
      if(fabs(ad.E-adREF.E) > THRESHOLD/mc.energy().atomDataMR.size()) {//error
	fout << "Error[atomMR.E]:" << endl
	     << "|" << ad.E << "-" << adREF.E << "| = " 
	     << fabs(ad.E-adREF.E) << endl;
	fout << "atomMRID(host) = " << k << endl
	     << "atomMRID(moon) = " << kp << endl;
	fout << "atomMR.x(host) = "; dumpV(a->x,3,fout) << endl;
	fout << "atomMR.x(moon) = "; dumpV(ap->x,3,fout) << endl;
	fout << "E = " << ad.E << " ~ " << adREF.E << endl;
	fout << "--------------------" << endl;
	return;
      }
      else if(fabs(ad.E-adREF.E) > maxdE) maxdE = fabs(ad.E-adREF.E);
    }//end for(int j...
  }//end for(int k...
  fout << "maxdE = " << maxdE << endl;
}
