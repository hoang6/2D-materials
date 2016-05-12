#include "../Util/util.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 4 && argc != 5) {
    cout << "Usage: analMCNPH <stat.log> <skip> <drec> [digits = 9]" << endl;
    return 1;
  }

  /********** READ ARGUMENTS **********/
  ifstream sfin;
  int skip;
  int drec;
  int digits;

  sfin.open(argv[1]);
  if(!sfin.good()) {
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
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(digits);
  cout << left;
  
  /********** DO ANALYSIS **********/
  int sN;
  double sP;
  double sH;
  double sPot, sEk, sh[9];
  double h0[9], inv_h0[9];
  double strain[9];
  double dataEkStrain[10];
  double meanEkStrain[10];
  double covEkStrain[10][10];
  double mtmp[10], mtmp_tr[10];
  int j, k, trigger, count, c;
  string keyword;
  streampos spos;

  /* Read header */
  keyword = "";
  while(keyword != "N" && sfin.good()) sfin >> keyword;
  sfin >> sN;
  sfin >> keyword >> sP;
  sfin >> keyword >> sH;
  for(k = 0; k < 11; k++) sfin >> keyword; //skip table header
  if(sfin.fail()) {
    cout << "Error: fail to read header" << endl;
    return 6;
  }

  /* First pass, get h0 */ 
  trigger = 0;
  count = 0;
  spos = sfin.tellg();
  while(true) {
    sfin >> sPot >> sEk;
    for(k = 0; k < 9; k++) sfin >> sh[k];
    if(sfin.fail()) break;
    if(trigger++ < skip) continue;

    if(count%drec == 0) {
      c = count/drec;
      if(c == 0) {
	Util::vcopy(sh,h0,9);
      }
      else {
	for(k = 0; k < 9; k++) 
	  h0[k] = sh[k]/(c+1)+(c*h0[k])/(c+1);
      }
    }
    count++;
  }
  Util::minverse(h0, inv_h0, 3);

  /* Second pass, get cov(Ek, strain) */
  trigger = 0;
  count = 0;
  sfin.clear();
  sfin.seekg(spos);
  while(true) {
    sfin >> sPot >> sEk;
    for(k = 0; k < 9; k++) sfin >> sh[k];
    if(sfin.fail()) break;
    if(trigger++ < skip) continue;
    
    if(count%drec == 0) {
      c = count/drec;
      //compute finite strain
      Util::mdotm(sh, inv_h0, mtmp, 3, 3, 3);
      Util::mtrans(mtmp, mtmp_tr, 3, 3);
      Util::mdotm(mtmp_tr, mtmp, strain, 3, 3, 3);
      for(k = 0; k < 3; k++) strain[4*k] -= 1.0;
      Util::kdotm(0.5, strain, strain, 3, 3);
      //compute cov(Ek, strain)
      dataEkStrain[0] = sEk;
      for(k = 1; k < 10; k++) dataEkStrain[k] = strain[k-1];
      if(c == 0) {
	for(k = 0; k < 10; k++) meanEkStrain[k] = dataEkStrain[k];
	for(j = 0; j < 10; j++)
	  for(k = 0; k < 10; k++) covEkStrain[j][k] = 0.0;
      }
      else {
	for(k = 0; k < 10; k++) {
	  mtmp[k] = meanEkStrain[k];
	  meanEkStrain[k] = 
	    dataEkStrain[k]/(c+1)+(c*meanEkStrain[k])/(c+1);
	}//end for(k...
	for(j = 0; j < 10; j++) {
	  for(k = j; k < 10; k++) {
	    covEkStrain[j][k] = 
	      (dataEkStrain[j]-meanEkStrain[j])*(dataEkStrain[k]-mtmp[k])/c+
	      ((c-1)*covEkStrain[j][k])/c;
	  }//end for(k...
	}//end for(j...
      }
    }
    count++;
  }//end while(sfin...
  for(j = 0; j < 10; j++)
    for(k = 0; k < j; k++) covEkStrain[j][k] = covEkStrain[k][j];

  /********** FINALIZE **********/
  cout << "skip       " << skip << endl;
  cout << "nrecords   " << count << endl;
  cout << "drecords   " << drec << endl;
  cout << "N          " << sN << endl;
  cout << "P          " << showpos << sP << endl;
  cout << "H          " << showpos << sH << endl;
  cout << "Ek_mean    " << showpos << meanEkStrain[0] << endl;
  cout << "h0         ";
  for(k = 0; k < 9; k++) {
    cout.width(digits+9);
    cout << showpos << h0[k] << " ";
  }
  cout << endl;
  cout << "cov(Ek,e)" << endl;
  for(j = 0; j < 10; j++) {
    for(k = 0; k < 10; k++) {
      cout.width(digits+9);
      cout << showpos << covEkStrain[j][k] << " ";
    }
    cout << endl;
  }
  cout << endl;

  return 0;
}
