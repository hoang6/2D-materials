#include "../AMod/AMod.h"
#include "../Util/util.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void molInfo(const AMod::Molecule& _mol, 
	     double _h[9], double& _bondr, double& _r);

int main(int argc, char* argv[]) {
  if(argc != 4 && argc != 5) {
    cout << "Usage: analGPH <file_mols> <skip> <drec> [digits = 9]" << endl;
    return 1;
  }

  /********** READ ARGUMENTS **********/
  ifstream mfin;
  int skip;
  int drec;
  int digits;

  mfin.open(argv[1]);
  if(!mfin.good()) {
    cout << "Error: cannot open file " << argv[1] << endl;
    return 2;
  }
  if(!Util::readWord(argv[2],skip)) {
    cout << "Error: illegal argument: skip" << endl;
    return 3;
  }
  if(!Util::readWord(argv[3],drec)) {
    cout << "Error: illegal argument: drec" << endl;
    return 4;
  }
  if(argc == 5) {
    if(!Util::readWord(argv[4], digits, 9)) {
      cout << "Error: illegal argument: digits" << endl;
      return 5;
    }
  }
  AMod::MolO::setFormat(cout,digits);

  /********** DO ANALYSIS **********/
  AMod::Molecule mol;
  double data[11];
  double data_mean[11], data_var[11];
  int k, trigger, count, c;
  double tmp, v;

  trigger = 0;
  count = 0;
  while(!mol.io().sourceTxt(mfin)) {
    if(trigger < skip) { trigger++; continue; }

    //get mol information
    molInfo(mol,data,data[9],data[10]);

    if(count%drec == 0) {
      c = count/drec;
      if(c == 0) {
	Util::vcopy(data,data_mean,11);
	Util::vassign(0.0,data_var,11);
      }
      else {
	for(k = 0; k < 11; k++) {
	  tmp = data_mean[k];
	  v = data[k];
	  data_mean[k] = 
	    v/(c+1)+(c*(data_mean[k]))/(c+1);
	  data_var[k] = 
	    (v-(data_mean[k]))*(v-tmp)/c+((c-1)*(data_var[k]))/c;
	}
      }
    }
    //counter
    count++;    
  }//end while(...
  mfin.close();
  
  /********** FINALIZE **********/
  cout << "skip       " << skip << endl;
  cout << "nrecords   " << count << endl;
  cout << "drecords   " << drec << endl;
  cout << "natoms     " << mol.natoms() << endl;
  cout << "axes_mean  ";
  for(k = 0; k < 9; k++) cout << data_mean[k] << " "; cout << endl;
  cout << "axes_var   ";
  for(k = 0; k < 9; k++) cout << data_var[k] << " "; cout << endl;
  cout << "bondr_mean " << data_mean[9] << endl;
  cout << "bondr_var  " << data_var[9] << endl;
  cout << "r_mean     " << data_mean[10] << endl;
  cout << "r_var      " << data_var[10] << endl;

  return 0;
}

void molInfo(const AMod::Molecule& _mol, 
	     double _h[9], double& _bondr, double& _r) {
  int k, natoms = _mol.natoms(), nbonds = _mol.nbonds();
  double tmp;

  for(k = 0; k < 3; k++)
    Util::vcopy(_mol.axis(k).x, _h+3*k, 3);

  _bondr = 0.0;
  for(k = 0; k < nbonds; k++) 
    _bondr += _mol.bond(k).r;
  _bondr /= nbonds;

  _r = 0.0;
  AMod::MolAdjuster::Point mc = AMod::MolAdjuster::massCenter(_mol);
  for(k = 0; k < natoms; k++) {
    tmp = 
      Util::square(_mol.atom(k)->x[0]-mc[0])+
      Util::square(_mol.atom(k)->x[1]-mc[1]);
    _r += std::sqrt(tmp);
  }
  _r /= natoms;
}
