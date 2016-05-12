#ifndef _UTIL_NEB_H
#define _UTIL_NEB_H

#include "Accessors.h"
#include "util.h"
#include <cstddef>
#include <vector>

namespace Util {

template<class F, class T = double> class NEB;  

/******************* The interface of F ********************
void operator()(T* x, T* g, T* energy);
where:
length(x) = (imageNumber-2)*imageSize
length(g) = (imageNumber-2)*imageSize
length(energy) = imageNumber-2
************************************************************/

/************************************************************/
template<class F, class T>
class NEB {
public:
  NEB(T _sprConstMin, T _sprConstMax,
      int _imageSZ, int _imageNO, F& _funcd);
  virtual ~NEB();
  void load_xa(const T* _xa, const T& _energy);
  void load_xz(const T* _xz, const T& _energy);
  void load_xi(const T* _xi, int i, const T& _energy);
  void load_mass(const T* _mass);
  void load_fixed(const int* _fixed);
  void load_xi();
  void load_mass(T _mass = 1.0);
  void load_fixed(int _fixed = 0);
  double operator()(T* _x);
  void operator()(T* _x, T* _g);
  void operator()(T* _x, T* _g, T _f);
  virtual void initialize();
  virtual void finalize();
  //set and get
  T* work_x() const;
  const T* work_mass() const;
  int work_size() const;
  int imageSize() const;
  int imageNumber() const;
  int climbImageIndex() const;
  T climbImageEnergy() const;
  int peakImageIndex() const;
  T peakImageEnergy() const;
  template<class U>
  U* seek_i(int _kimage, U* _p);
  T* seek_xi(int _kimage, T* _x);
  T* seek_xi(int _kimage);

  Accessors<bool> climb;
  Accessors<std::ostream*> callBackStreamPtr;

protected:
  void intrp();
  void intrp_ij(int i, int j);
  void compdR(int i, int j, T* _x, T* dR);
  void compEnergyForce(T* _x, T* _g);
  void compSprConsts();
  void compTau(T* _x);
  virtual void dfunc(T* _x, T* _g);
  virtual void callBack(T* _x, T* _g, T _f);

  T sprConstMin, sprConstMax;
  int climbImageInd;
  F* pfuncd;                      //provide func, dfunc for one image
  int imageSZ, imageNO;           //imageNO includes start and end images
  int nDOF;                       //nDOF = imageSize*(imageNO-2)
  T* xa, *xz;                     //start image and end image, size = imageSize
  int* fixed_i;                   //size = imageSize
  T* tauTemp;                     //size = imageSize
  T* x;                           //size = nDOF
  T* mass;                        //size = nDOF
  T* tau;                         //size = nDOF
  T* energy;                      //size = imageNO
  T* sprConsts;                   //size = imageNO-1
  int iter;

private:
  int maxEInd;
  T maxE;
  std::vector<int> imageIndices;

  NEB(const NEB& neb);
  NEB& operator= (const NEB& neb);
};

/************************************************************/
template<class F, class T>
NEB<F, T>::NEB(T _sprConstMin, T _sprConstMax,
	       int _imageSZ, int _imageNO, F& _funcd) {
  sprConstMin = _sprConstMin;
  sprConstMax = _sprConstMax;
  pfuncd = &_funcd;
  imageSZ = _imageSZ;
  imageNO = _imageNO;
  nDOF = imageSZ*(imageNO-2);

  tauTemp = new T [imageSZ];
  xa = new T [imageSZ];
  xz = new T [imageSZ];
  fixed_i = new int [imageSZ];

  x = new T [nDOF];
  mass = new T[nDOF];
  tau = new T [nDOF];

  energy = new T [imageNO];
  sprConsts = new T [imageNO-1];

  climb(false);
  callBackStreamPtr((std::ostream*)NULL);

  load_xi();
  load_mass();
  load_fixed();
}

template<class F, class T>
NEB<F, T>::~NEB() {
  delete [] sprConsts;
  delete [] energy;
  delete [] tau;
  delete [] mass;
  delete [] x;
  delete [] fixed_i;
  delete [] xz;
  delete [] xa;
  delete [] tauTemp;
}

template<class F, class T>
inline void NEB<F, T>::load_xa(const T* _xa, const T& _energy) {
  vcopy(_xa, xa, imageSZ);
  energy[0] = _energy;
  
}

template<class F, class T>
inline void NEB<F, T>::load_xz(const T* _xz, const T& _energy) {
  vcopy(_xz, xz, imageSZ);
  energy[imageNO-1] = _energy;
}

template<class F, class T>
void NEB<F, T>::load_xi(const T* _xi, int i, const T& _energy) {
  T* xi = seek_xi(i,x);
  if(xi == (T*)NULL) return;
  vcopy(_xi, xi, imageSZ);
  imageIndices.push_back(i);
  energy[i] = _energy;
}

template<class F, class T>
inline void NEB<F, T>::load_mass(const T* _mass) {
  for(int kimage = 1; kimage <= imageNO-2; kimage++)
    vcopy(_mass, seek_i(kimage, mass), imageSZ);
}

template<class F, class T>
inline void NEB<F, T>::load_fixed(const int* _fixed) {
  vcopy(_fixed, fixed_i, imageSZ);
}

template<class F, class T>
inline void NEB<F, T>::load_xi() { imageIndices.clear(); }

template<class F, class T>
void NEB<F, T>::load_mass(T _mass) {
  vassign(_mass, mass, nDOF);
}

template<class F, class T>
inline void NEB<F, T>::load_fixed(int _fixed) {
  vassign(_fixed, fixed_i, imageSZ);
}

template<class F, class T>
inline double NEB<F, T>::operator()(T* _x) {
  return 0.0;
}

template<class F, class T>
inline void NEB<F, T>::operator()(T* _x, T* _g) { 
  dfunc(_x, _g); 
}

template<class F, class T>
inline void NEB<F, T>::operator()(T* _x, T* _g, T _f) {
  callBack(_x, _g, _f);
}

template<class F, class T>
void NEB<F, T>::initialize() {
  intrp();
  climbImageInd = -1;
  iter = 0;
}

template<class F, class T>
inline void NEB<F, T>::finalize() {}

template<class F, class T>
inline T* NEB<F, T>::work_x() const { return x; }

template<class F, class T>
inline const T* NEB<F, T>::work_mass() const { return mass; }

template<class F, class T>
inline int NEB<F, T>::work_size() const { return nDOF; }

template<class F, class T>
inline int NEB<F, T>::imageSize() const { return imageSZ; }

template<class F, class T>
inline int NEB<F, T>::imageNumber() const { return imageNO; }

template<class F, class T>
inline int NEB<F, T>::climbImageIndex() const { return climbImageInd; }

template<class F, class T>
inline T NEB<F, T>::climbImageEnergy() const { 
  return (climbImageInd<0 ? peakImageEnergy() : energy[climbImageInd]); 
}

template<class F, class T>
inline int NEB<F, T>::peakImageIndex() const { return maxEInd; }

template<class F, class T>
inline T NEB<F, T>::peakImageEnergy() const { return maxE; }

template<class F, class T>
template<class U>
U* NEB<F, T>::seek_i(int _kimage, U* _p) { 
  if(_kimage <= 0) return (U*)NULL;
  else if(_kimage >= imageNO-1) return (U*)NULL;
  return _p+(_kimage-1)*imageSZ; 
}

template<class F, class T>
T* NEB<F, T>::seek_xi(int _kimage, T* _x) {
  if(_kimage == 0) return xa;
  else if(_kimage == imageNO-1) return xz;
  return seek_i(_kimage, _x);
}

template<class F, class T>
inline T* NEB<F, T>::seek_xi(int _kimage) { return seek_xi(_kimage, x); }

template<class F, class T>
void NEB<F, T>::intrp() {
  int k, sz = imageIndices.size();
  if(sz == 0) {
    intrp_ij(0, imageNO-1);
  }
  else {
    intrp_ij(0,imageIndices[0]);
    for(k = 0; k < sz-1; k++)
      intrp_ij(imageIndices[k],imageIndices[k+1]);
    intrp_ij(imageIndices[sz-1],imageNO-1);
  }
}

template<class F, class T>
void NEB<F, T>::intrp_ij(int i, int j) {
  T* xStart = seek_xi(i,x);
  T* xEnd = seek_xi(j,x);
  T* xi;
  T frac;
  for(int kimage = i+1; kimage <= j-1; kimage++) {
    xi = seek_xi(kimage,x);
    frac = T(kimage-i)/(j-i);
    for(int k = 0; k < imageSZ; k++)
      xi[k] = (1-frac)*xStart[k]+frac*xEnd[k];
  }
}

template<class F, class T>
inline void NEB<F, T>::compdR(int i, int j, T* _x, T* dR) {
  //dR = Ri-Rj
  kdotvadd(1.0, seek_xi(i, _x), -1.0, seek_xi(j, _x), dR, imageSZ);
}

template<class F, class T>
void NEB<F, T>::compEnergyForce(T* _x, T* _g) {
  int kimage;
  
  (*pfuncd)(_x, _g, energy+1); //kimage = 1~imageNO-2

  for(kimage = 0; kimage <= imageNO-1; kimage++) {
    if(kimage == 0 || maxE < energy[kimage]) {
      maxEInd = kimage;
      maxE = energy[kimage];
    }
  } 

  if(climb() && climbImageInd < 0) climbImageInd = maxEInd;

  /* debug
  if(climb())
    std::cout << "climb is triggered: " <<  
	      << "image " << climbImageInd 
	      << std::endl;
  */
}

template<class F, class T>
void NEB<F, T>::compSprConsts() {
  if(sprConstMax == sprConstMin) 
    vassign(sprConstMax, sprConsts, imageNO-1);    

  T kmax = sprConstMax, dk = sprConstMax-sprConstMin;
  T Ei, Eref = MAX(energy[0], energy[imageNO-1]);
  for(int kimage = 1; kimage <= imageNO-1; kimage++) {
    Ei = MAX(energy[kimage], energy[kimage-1]);
    if(Ei > Eref) 
      sprConsts[kimage-1] = kmax-dk*((maxE-Ei)/(maxE-Eref));
    else
      sprConsts[kimage-1] = sprConstMin;
  }
}

template<class F, class T>
void NEB<F, T>::compTau(T* _x) {
  T Ei, Epre, Epost;
  T dEmax, dEmin;
  T* tau_i;
  for(int kimage = 1; kimage <= imageNO-2; kimage++) {
    Ei = energy[kimage];
    Epre = energy[kimage-1];
    Epost = energy[kimage+1];
    tau_i = seek_i(kimage, tau);
    if(Epost > Ei && Ei > Epre) 
      compdR(kimage+1, kimage, _x, tau_i);
    else if(Epost < Ei && Ei < Epre) 
      compdR(kimage, kimage-1, _x, tau_i);
    else if(Epost > Epre) {
      dEmax = MAX(fabs(Epost-Ei), fabs(Ei-Epre));
      dEmin = MIN(fabs(Epost-Ei), fabs(Ei-Epre));
      compdR(kimage+1, kimage, _x, tau_i);
      compdR(kimage, kimage-1, _x, tauTemp);
      kdotvadd(dEmax, tau_i, dEmin, tauTemp, 
		   tau_i, imageSZ);
    }
    else {
      dEmax = MAX(fabs(Epost-Ei), fabs(Ei-Epre));
      dEmin = MIN(fabs(Epost-Ei), fabs(Ei-Epre));
      compdR(kimage+1, kimage, _x, tau_i);
      compdR(kimage, kimage-1, _x, tauTemp);
      kdotvadd(dEmin, tau_i, dEmax, tauTemp, 
		   tau_i, imageSZ);
    }
    kdotv(1.0/sqrt(norm2(tau_i, imageSZ)), tau_i, imageSZ);
  }//end for(kimage = 1...
}

template<class F, class T>
void NEB<F, T>::dfunc(T* _x, T* _g) {
  //must be done in correct order
  compEnergyForce(_x, _g);
  compSprConsts();
  compTau(_x);

  int k, kimage;
  T temp;
  T* tau_i, *gi;
  for(kimage = 1; kimage <= imageNO-2; kimage++) {
    tau_i = seek_i(kimage, tau);
    gi = seek_i(kimage, _g);
    temp = vdotv(gi, tau_i, imageSZ);
    if(climb() && kimage == climbImageInd) {	  
      kdotvadd(+1.0, gi, -2*temp, tau_i, gi, imageSZ);	  
    }
    else {
      kdotvadd(+1.0, gi, -temp, tau_i, 
		   gi, imageSZ);
      compdR(kimage+1, kimage, _x, tauTemp);
      temp = sprConsts[kimage]*sqrt(vdotv(tauTemp, tauTemp, imageSZ));
      compdR(kimage, kimage-1, _x, tauTemp);
      temp -= sprConsts[kimage-1]*sqrt(vdotv(tauTemp, tauTemp, imageSZ));
      for(k = 0; k < imageSZ; k++) 
	tauTemp[k] = (fixed_i[k] ? 0.0 : temp*tau_i[k]);
      kdotvadd(+1.0, gi, -1.0, tauTemp, gi, imageSZ);
    }
  }//end for(kimage = 1...
}

template<class F, class T>
void NEB<F, T>::callBack(T* _x, T* _g, T _f) {
  int kimage, k;
  double* xout;

  iter++;

  if(callBackStreamPtr() != (std::ostream*)NULL) {
    std::ostream& fout = *callBackStreamPtr();
    fout << "iter " << iter << std::endl;
    for(kimage = 0; kimage < imageNO; kimage++) {
      fout << kimage << " " << energy[kimage] << " ";
      xout = seek_xi(kimage, _x);
      for(k = 0; k < imageSZ; k++) fout << xout[k] << " ";
      fout << std::endl;
    }
    fout << std::endl;
  }/* end if(callBackStram()... */

  /* debug
  std::cout << "<NEB::callBack> " 
	    << iter << " " 
	    << norm1(_g, work_size()) 
	    << std::endl;
  */
}

}/* Util */

#endif
