#ifndef _AMOD_TRANS_H
#define _AMOD_TRANS_H

#include "../Util/Accessors.h"

namespace AMod {

/* Follow the minimum image convention */

class Trans {
public:
  Trans();
  template<class M>
  Trans(int _mode, const M& _mat);
  Trans(const Trans& trans);
  Trans& operator= (const Trans& trans);
  void reset();
  template<class M>
  void init(int _mode, const M& _mat);
  void proj(const double* _dr, double* _x) const; 
  /*
    find the image of myPos
  */
  void findImage(double* _myPos, const double* _taPos) const;
  /*
    find the image of myPos nearest to taPos and within cutoff
  */
  bool findImage(double* _myPos, const double* _taPos, double _cutoff) const;
  /* get */
  int mode() const;
  const double* displacement() const;
  double displacement(int i) const;
  double distSquare() const;
  int nT(int i) const;
  const int* nT() const;

protected:
  void init1(const double* _T1);
  void init2(const double* _T1, const double* _T2);
  void init3(const double* _T1, const double* _T2, const double* _T3);
  bool findImage0(double* _myPos, const double* _taPos, double _cutoff) const;
  bool findImage1(double* _myPos, const double* _taPos, double _cutoff) const;
  bool findImage2(double* _myPos, const double* _taPos, double _cutoff) const;
  bool findImage3(double* _myPos, const double* _taPos, double _cutoff) const;

  int mymode;
  double T1[3];
  double T2[3];
  double T3[3];
  double TS1[3];
  double TS2[3];
  double TS3[3];
  double T11, T22, T33;
  double TS11, TS22, TS33;
  double cos11, cos22, cos33;

private:
  static const double SAFE_COEF;
  mutable double disp[3]; //set by findImage*; read only
  mutable double distSqr; //set by findImage*; read only
  mutable int n[3];       //set by findImage*; read only
};

/************************************************************/
template<class M>
Trans::Trans(int _mode, const M& _mat) { init(_mode, _mat); }

template<class M>
void Trans::init(int _mode, const M& _mat) {
  reset();
  mymode = _mode;
  if(mode() == 1) init1(&_mat[0][0]);
  else if(mode() == 2) init2(&_mat[0][0], &_mat[1][0]);
  else if(mode() == 3) init3(&_mat[0][0], &_mat[1][0], &_mat[2][0]);
}

inline int Trans::mode() const { return mymode; }

inline const double* Trans::displacement() const { return disp; }

inline double Trans::displacement(int i) const { return disp[i]; }

inline double Trans::distSquare() const { return distSqr; }

inline int Trans::nT(int i) const { return n[i]; }

inline const int* Trans::nT() const { return n; }

} /* AMod */

#endif
