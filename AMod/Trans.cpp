#include <iostream>
#include "Trans.h"
#include "../Util/util.h"
#include <cmath>

namespace AMod {

const double Trans::SAFE_COEF = 1.01;

Trans::Trans() { reset(); }

Trans::Trans(const Trans& trans) { *this = trans; }

Trans& Trans::operator= (const Trans& trans) { 
  if(this != &trans) {
    mymode = trans.mymode;
    Util::vcopy(trans.T1, T1, 3);
    Util::vcopy(trans.T2, T2, 3);
    Util::vcopy(trans.T3, T3, 3);
    Util::vcopy(trans.TS1, TS1, 3);
    Util::vcopy(trans.TS2, TS2, 3);
    Util::vcopy(trans.TS3, TS3, 3);
    T11 = trans.T11;
    T22 = trans.T22;
    T33 = trans.T33;
    TS11 = trans.TS11;
    TS22 = trans.TS22;
    TS33 = trans.TS33;
    cos11 = trans.cos11;
    cos22 = trans.cos22;
    cos33 = trans.cos33;
  }
  return *this; 
}

void Trans::reset() { 
  mymode = 0; 
  Util::vassign(0.0, T1, 3);
  Util::vassign(0.0, T2, 3);
  Util::vassign(0.0, T3, 3);
  Util::vassign(0.0, TS1, 3);
  Util::vassign(0.0, TS2, 3);
  Util::vassign(0.0, TS3, 3);
  T1[0] = T2[1] = T3[2] = 1.0;
  TS1[0] = TS2[1] = TS3[2] = 1.0;
  T11 = T22 = T33 = 1.0;
  TS11 = TS22 = TS33 = 1.0;
  cos11 = cos22 = cos33 = 1.0;
  Util::vassign(0.0, disp, 3);
  distSqr = 0.0;
  Util::vassign(0, n, 3);
}

void Trans::proj(const double* _dr, double* _x) const {
  if(mode() == 1) {
    _x[0] = Util::vdotv(_dr, TS1, 3);
  }
  else if(mode() == 2 || mode() == 3) {
    _x[0] = Util::vdotv(_dr, TS1, 3);
    _x[1] = Util::vdotv(_dr, TS2, 3);
    _x[2] = Util::vdotv(_dr, TS3, 3);
  }
}

void Trans::findImage(double* _myPos, const double* _taPos) const {
  findImage(_myPos, _taPos, -1.0);
}

bool Trans::findImage(double* _myPos, const double* _taPos, double _cutoff) const {
  if(mode() == 0) return findImage0(_myPos, _taPos, _cutoff);
  else if(mode() == 1) return findImage1(_myPos, _taPos, _cutoff);
  else if(mode() == 2) return findImage2(_myPos, _taPos, _cutoff);
  else if(mode() == 3) return findImage3(_myPos, _taPos, _cutoff);
  return false;
}

bool Trans::findImage0(double* _myPos, const double* _taPos, double _cutoff) const {
  if(_cutoff < 0) {
    Util::vassign(0, n, 3);
    Util::vsub(_myPos, _taPos, disp, 3);
    distSqr = Util::vdotv(disp, disp, 3);
  }
  else {
    Util::vassign(0, n, 3);
    double rc2 = _cutoff*_cutoff;
    double dr[3];
    //for efficiency
    distSqr = 0.0;
    for(int k = 0; k < 3; k++) {
      dr[k] = _myPos[k]-_taPos[k];
      distSqr += dr[k]*dr[k];
      if(distSqr-rc2 > Util::EPS_DOUBLE) return false;
    }
    Util::vsub(_myPos, _taPos, disp, 3);
  }
  return true;
}

bool Trans::findImage1(double* _myPos, const double* _taPos, double _cutoff) const {
  n[1] = n[2] = 0;
  double rc2;
  double dr[3], x[3];
  Util::vsub(_myPos, _taPos, dr, 3);
  //project dr on T1
  proj(dr, x);
  double k1 = x[0];
  double res2 = Util::vdotv(dr, dr, 3)-k1*k1*T11;
  if(_cutoff < 0) {
    double _dr[3];
    Util::kdotvadd(1.0, dr, -std::floor(k1+0.5), T1, _dr, 3);
    rc2 = Util::vdotv(_dr, _dr, 3)*SAFE_COEF;
  }
  else {
    rc2 = _cutoff*_cutoff;
  }
  if(res2-rc2 > Util::EPS_DOUBLE) return false;
  //draw cutoff on T1
  double k1max = sqrt((rc2-res2)/T11);
  //find all images
  int _n1;
  int n1min = int(ceil(-k1max-k1));
  int n1max = int(floor(k1max-k1));
  if(n1min > n1max) return false;
  double _distSqr;
  bool flag = 1;
  for(_n1 = n1min; _n1 <= n1max; _n1++) {
    _distSqr = res2+(k1+_n1)*(k1+_n1)*T11;
    if(flag || _distSqr < distSqr) {
      n[0] = _n1;
      distSqr = _distSqr;
      flag = 0;
    } 
  }
  if(distSqr-rc2 > Util::EPS_DOUBLE) return false;
  Util::kdotvadd(1.0, _myPos, 
		 double(n[0]), T1, 
		 _myPos, 3);
  Util::vsub(_myPos, _taPos, disp, 3);
  return true;
}

bool Trans::findImage2(double* _myPos, const double* _taPos, double _cutoff) const {
  n[2] = 0;
  double rc2;
  double dr[3], x[3];
  Util::vsub(_myPos, _taPos, dr, 3);
  proj(dr, x);
  double k1 = x[0];
  double k2 = x[1];
  double res2 = x[2]*x[2];
  if(_cutoff < 0) {
    double _dr[3];
    Util::kdotvadd3(1.0, dr, 
		    -std::floor(k1+0.5), T1, 
		    -std::floor(k2+0.5), T2, _dr, 3);
    rc2 = Util::vdotv(_dr, _dr, 3)*SAFE_COEF;
  }
  else {
    rc2 = _cutoff*_cutoff;
  }
  if(res2-rc2 > Util::EPS_DOUBLE) return false;
  //draw cutoff on (T1, T2)
  double k1max = std::sqrt((rc2-res2)/T11)/cos11;
  double k2max = std::sqrt((rc2-res2)/T22)/cos22;
  int _n1, _n2;
  int n1min = int(ceil(-k1max-k1));
  int n1max = int(floor(k1max-k1));
  int n2min = int(ceil(-k2max-k2));
  int n2max = int(floor(k2max-k2));
  if(n1min > n1max || n2min > n2max) return false;
  double _distSqr, tmp[3];
  bool flag = 1;
  for(_n1 = n1min; _n1 <= n1max; _n1++) 
    for(_n2 = n2min; _n2 <= n2max; _n2++) {
      Util::kdotvadd(double(k1+_n1), T1, double(k2+_n2), T2, tmp, 3);
      _distSqr = Util::vdotv(tmp, tmp, 3)+res2;
      if(flag || _distSqr < distSqr) {
	n[0] = _n1;
	n[1] = _n2;
	distSqr = _distSqr;
	flag = 0;
      }
    }
  if(distSqr-rc2 > Util::EPS_DOUBLE) return false;
  Util::kdotvadd3(1.0, _myPos, 
		  double(n[0]), T1, 
		  double(n[1]), T2, 
		  _myPos, 3);
  Util::vsub(_myPos, _taPos, disp, 3);
  return true;
}

bool Trans::findImage3(double* _myPos, const double* _taPos, double _cutoff) const {
  double rc2;
  double dr[3], x[3];
  Util::vsub(_myPos, _taPos, dr, 3);
  proj(dr, x);
  double k1 = x[0];
  double k2 = x[1];
  double k3 = x[2];
  if(_cutoff < 0) {
    double _dr[3];
    Util::kdotvadd4(1.0, dr, 
		    -std::floor(k1+0.5), T1, 
		    -std::floor(k2+0.5), T2, 
		    -std::floor(k3+0.5), T3, _dr, 3);
    rc2 = Util::vdotv(_dr, _dr, 3)*SAFE_COEF;
  }
  else {
    rc2 = _cutoff*_cutoff;
  }
  //draw cutoff on (T1, T2)
  double k1max = std::sqrt(rc2/T11)/cos11;
  double k2max = std::sqrt(rc2/T22)/cos22;
  double k3max = std::sqrt(rc2/T33)/cos33;
  int _n1, _n2, _n3;
  int n1min = int(ceil(-k1max-k1));
  int n1max = int(floor(k1max-k1));
  int n2min = int(ceil(-k2max-k2));
  int n2max = int(floor(k2max-k2));
  int n3min = int(ceil(-k3max-k3));
  int n3max = int(floor(k3max-k3));
  if(n1min > n1max || n2min > n2max || n3min > n3max) return false;
  double _distSqr, tmp[3];
  bool flag = 1;
  for(_n1 = n1min; _n1 <= n1max; _n1++) 
    for(_n2 = n2min; _n2 <= n2max; _n2++)
      for(_n3 = n3min; _n3 <= n3max; _n3++) {
	Util::kdotvadd3(double(k1+_n1), T1, 
			double(k2+_n2), T2, 
			double(k3+_n3), T3, 
			tmp, 3);
	_distSqr = Util::vdotv(tmp, tmp, 3);
	if(flag || _distSqr < distSqr) {
	  n[0] = _n1;
	  n[1] = _n2;
	  n[2] = _n3;
	  distSqr = _distSqr;
	  flag = 0;
	}
      }
  if(distSqr-rc2 > Util::EPS_DOUBLE) return false;
  Util::kdotvadd4(1.0, _myPos, 
		  double(n[0]), T1, 
		  double(n[1]), T2, 
		  double(n[2]), T3, 
		  _myPos, 3);
  Util::vsub(_myPos, _taPos, disp, 3); 
  return true;
}

void Trans::init1(const double* _T1) {
  Util::vcopy(_T1, T1, 3);
  T11 = Util::vdotv(T1, T1, 3);
  Util::kdotv(1.0/T11, T1, TS1, 3);
}

void Trans::init2(const double* _T1, const double* _T2) {
  double _T3[3];
  Util::cross3(_T1, _T2, _T3);
  Util::normalize(_T3, 3);
  init3(_T1, _T2, _T3);
}

void Trans::init3(const double* _T1, const double* _T2, const double* _T3) {
  Util::vcopy(_T1, T1, 3);
  Util::vcopy(_T2, T2, 3);
  Util::vcopy(_T3, T3, 3);
  Util::dual_axes(T1, T2, T3, TS1, TS2, TS3);
  T11 = Util::vdotv(T1, T1, 3);
  T22 = Util::vdotv(T2, T2, 3);
  T33 = Util::vdotv(T3, T3, 3);
  TS11 = Util::vdotv(TS1, TS1, 3);
  TS22 = Util::vdotv(TS2, TS2, 3);
  TS33 = Util::vdotv(TS3, TS3, 3);
  cos11 = 1.0/std::sqrt(T11*TS11);
  cos22 = 1.0/std::sqrt(T22*TS22);
  cos33 = 1.0/std::sqrt(T33*TS33);
}

} /* AMod */
