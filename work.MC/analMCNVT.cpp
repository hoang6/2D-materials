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
    cout << "Usage: analMCNVT <stat.log> <skip> <drec> [digits = 9]" << endl;
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
  double sV;
  double skT;
  double data[10];
  double data_mean[10], data_var[10];
  int k, trigger, count, c;
  string keyword;
  double tmp, v;

  /* Read header */
  keyword = "";
  while(keyword != "N" && sfin.good()) sfin >> keyword;
  sfin >> sN;
  sfin >> keyword >> sV;
  sfin >> keyword >> skT;
  for(k = 0; k < 10; k++) sfin >> keyword; //skip table header
  if(sfin.fail()) {
    cout << "Error: fail to read header" << endl;
    return 6;
  }

  /* Get stat */
  trigger = 0;
  count = 0;
  sfin.clear();
  while(true) {
    for(k = 0; k < 10; k++) sfin >> data[k];
    if(sfin.fail()) break;
    if(trigger < skip) { trigger++; continue; }

    if(count%drec == 0) {
      c = count/drec;
      if(c == 0) {
	Util::vcopy(data,data_mean,10);
	Util::vassign(0.0,data_var,10);
      }
      else {
	for(k = 0; k < 10; k++) {
	  tmp = data_mean[k];
	  v = data[k];
	  data_mean[k] = 
	    v/(c+1)+(c*(data_mean[k]))/(c+1);
	  data_var[k] = 
	    (v-(data_mean[k]))*(v-tmp)/c+((c-1)*(data_var[k]))/c;
	}
      }
    }
    count++;
  }//end while(sfin...
  sfin.close();

  /********** FINALIZE **********/
  cout << "skip       " << skip << endl;
  cout << "nrecords   " << count << endl;
  cout << "drecords   " << drec << endl;
  cout << "N          " << sN << endl;
  cout << "V          " << sV << endl;
  cout << "kT         " << skT << endl;
  cout << "Epot_mean  " << showpos << data_mean[0] << endl;
  cout << "Epot_var   " << showpos << data_var[0] << endl;
  cout << "F_h_mean   ";
  for(k = 1; k < 10; k++) cout << data_mean[k] << " "; cout << endl;
  cout << "F_h_var    ";
  for(k = 1; k < 10; k++) cout << data_var[k] << " "; cout << endl;

  return 0;
}
