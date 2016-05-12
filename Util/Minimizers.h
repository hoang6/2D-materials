#ifndef _UTIL_MINIMIZERS_H
#define _UTIL_MINIMIZERS_H

#include "FuncDFunc.h"
#include "Verlet.h"
#include "QuenchVerlet.h"
#include "Fire.h"
#include "ConjGrad.h"
#include <string>
#include <iostream>

namespace Util {

template<class U, class T = double> class Minimizers;

template<class U, class T>
class Minimizers {
public:
  /* Minimizer Types */
  typedef U                     TypeFuncDFunc;
  typedef MinimizerBasic<U, T>  TypeBasic;
  typedef Verlet<U, T>          TypeVerlet;
  typedef QuenchVerlet<U, T>    TypeQuenchVerlet;
  typedef Fire<U, T>            TypeFire;
  typedef ConjGrad<U, T>        TypeConjGrad;
  /* Minimizer Labels */
  typedef int Selection;
  static const int VERLET;
  static const int QUENCH_VERLET;
  static const int FIRE;
  static const int CONJ_GRAD;
 
  Minimizers();
  Minimizers(const Minimizers& minimizers);
  virtual ~Minimizers();
  Minimizers& operator= (const Minimizers& minimizers);
  void reset();
  void set(Selection _sel, 
	   T _tol, int _maxIterNumber, T _deltaTime);
  TypeBasic& selected();
  const TypeBasic& selected() const;
  /********** source & dump **********/
  int sourceTxt(std::istream& fin);
  void dumpTxt(std::ostream& fout) const;

  Accessors<Selection>        selection;
  Accessors<TypeVerlet>       verlet;
  Accessors<TypeQuenchVerlet> quenchVerlet;
  Accessors<TypeFire>         fire;
  Accessors<TypeConjGrad>     conjGrad;

private:
  void setBasicPtrArray();
  
  TypeBasic* pbasic[4];
};

/************************************************************/
template<class U, class T>
const int Minimizers<U, T>::VERLET = 0;

template<class U, class T>
const int Minimizers<U, T>::QUENCH_VERLET = 1;

template<class U, class T>
const int Minimizers<U, T>::FIRE = 2;

template<class U, class T>
const int Minimizers<U, T>::CONJ_GRAD = 3;

template<class U, class T>
inline Minimizers<U, T>::Minimizers() { setBasicPtrArray(); reset(); }

template<class U, class T>
inline Minimizers<U, T>::Minimizers(const Minimizers<U, T>& minimizers) {
  setBasicPtrArray();
  *this = minimizers;
}

template<class U, class T>
inline Minimizers<U, T>::~Minimizers() {}

template<class U, class T>
Minimizers<U, T>& Minimizers<U, T>::operator= (const Minimizers<U, T>& minimizers) {
  if(this != &minimizers) {
    selection = minimizers.selection;
    verlet = minimizers.verlet;
    quenchVerlet = minimizers.quenchVerlet;
    fire = minimizers.fire;
    conjGrad = minimizers.conjGrad;
  }
  return *this;
}

template<class U, class T>
void Minimizers<U, T>::reset() {
  selection(VERLET);
  verlet().reset();
  quenchVerlet().reset();
  fire().reset();
  conjGrad().reset();
}

template<class U, class T>
void Minimizers<U, T>::set(Selection _sel, 
			   T _tol, int _maxIterNumber, T _deltaTime) {
  selection(_sel);
  selected().tolerance(_tol);
  selected().maxIterNumber(_maxIterNumber);
  selected().deltaTime(_deltaTime);
}

template<class U, class T>
inline typename Minimizers<U, T>::TypeBasic& Minimizers<U, T>::selected() { return *pbasic[selection()]; }

template<class U, class T>
inline const typename Minimizers<U, T>::TypeBasic& Minimizers<U, T>::selected() const { return *pbasic[selection()]; }

template<class U, class T>
int Minimizers<U, T>::sourceTxt(std::istream& fin) {
  std::string tword;
  if(!fin.good()) return 1;
  fin >> tword >> selection();
  if(!fin.good()) return 2;
  fin >> tword;
  if(verlet().sourceTxt(fin)) return 3;
  fin >> tword;
  fin >> tword;
  if(quenchVerlet().sourceTxt(fin)) return 4;
  fin >> tword;
  fin >> tword;
  if(fire().sourceTxt(fin)) return 5;
  fin >> tword;
  fin >> tword;
  if(conjGrad().sourceTxt(fin)) return 6;
  fin >> tword;
  return (fin.good() ? 0 : 1);
}

template<class U, class T>
void Minimizers<U, T>::dumpTxt(std::ostream& fout) const {
  using namespace std;
  fout << "selection " << selection() << endl;
  fout << "[verlet" << endl;
  verlet().dumpTxt(fout);
  fout << "]" << endl;
  fout << "[quenchVerlet" << endl;
  quenchVerlet().dumpTxt(fout);
  fout << "]" << endl;
  fout << "[fire" << endl;
  fire().dumpTxt(fout);
  fout << "]" << endl;
  fout << "[conjGrad"<< endl;
  conjGrad().dumpTxt(fout);
  fout << "]" << endl;
}

template<class U, class T>
void Minimizers<U, T>::setBasicPtrArray() {
  pbasic[VERLET] = &verlet();
  pbasic[QUENCH_VERLET] = &quenchVerlet();
  pbasic[FIRE] = &fire();
  pbasic[CONJ_GRAD] = &conjGrad();
}

}/* Util */

#endif
