#ifndef _UTIL_PTRCONT_H
#define _UTIL_PTRCONT_H

#include <vector>
#include <algorithm>
#include "utilSTL.h"

namespace Util {

template<class T>
class PtrCont {
public:
  class ID;

  PtrCont();
  PtrCont(const PtrCont& ptrCont);
  virtual ~PtrCont();
  PtrCont& operator= (const PtrCont& ptrCont);
  void cleanup();
  int size() const;
  void resize(int sz);
  void resize(int sz, const T& c);
  int capacity() const;
  bool empty() const;
  void reserve(int n);
  void shrink(double p = 1.0);
  T& operator[] (int n) const;
  T& front() const;
  T& back() const;
  void push_back();
  void push_back(const T& x);
  void pop_back();
  void clear();
  template<class Compare>
  void sort(const Compare& compare);
  void exch(int i, int j);
  T& insert();
  T& insert(const T& x);
  T& insert(int pos);
  T& insert(int pos, const T& x);
  void erase(int pos);
  template<class InputIterator>
  void keep(InputIterator first, InputIterator last);
  template<class InputIterator>
  void erase(InputIterator first, InputIterator last);

  class ID {
  public:
    ID() { identity = 0; }
    ID(const ID& id) {}
    ID& operator= (const ID& id) { return *this; }
    int operator() () const { return identity; }
  private:
    int identity;
    friend class PtrCont<T>;
  };

protected:
  template<class U, class Compare> class PtrCompare;
  std::vector<T*> cont;
  std::vector<T*> junk;
  
  void auto_shrink();

  template<class U, class Compare>
  class PtrCompare {
  public:
    PtrCompare(const Compare& _compare): compare(_compare) {}
    bool operator() (U* i, U* j) const { return compare(*i, *j); }
  private:
    const Compare& compare;
  };
};

/************************************************************/
template<class T>
inline PtrCont<T>::PtrCont() {}

template<class T>
inline PtrCont<T>::PtrCont(const PtrCont<T>& ptrCont) { *this = ptrCont; }

template<class T>
inline PtrCont<T>::~PtrCont() { cleanup(); }

template<class T>
PtrCont<T>& PtrCont<T>::operator= (const PtrCont<T>& ptrCont) {
  if(this != &ptrCont) {
    resize(ptrCont.size());
    for(int k = 0; k < size(); k++) *cont[k] = *ptrCont.cont[k];
  }
  return *this;
}

template<class T>
void PtrCont<T>::cleanup() {
  purge(cont);
  cont.clear();
  purge(junk);
  junk.clear();
}

template<class T>
inline int PtrCont<T>::size() const { return cont.size(); }

template<class T>
void PtrCont<T>::resize(int sz) {
  int k, old_size = size(), old_capacity = capacity();
  if(sz <= old_size) {
    for(k = sz; k < old_size; k++) *cont[k] = T();
    junk.insert(junk.end(), cont.end()-(old_size-sz), cont.end());
    cont.erase(cont.end()-(old_size-sz), cont.end());
  }
  else if(sz <= old_capacity) {
    cont.insert(cont.end(), junk.end()-(sz-old_size), junk.end());
    junk.erase(junk.end()-(sz-old_size), junk.end());
  }
  else {
    if(junk.size()) {
      cont.insert(cont.end(), junk.end()-(old_capacity-old_size), junk.end());
      junk.clear();
    }
    cont.resize(sz);
    for(k = old_capacity; k < sz; k++) cont[k] = new T();
  }
  for(k = old_size; k < sz; k++) { cont[k]->id.identity = k; }
  auto_shrink();
}

template<class T>
void PtrCont<T>::resize(int sz, const T& c) {
  int k, old_size = size();
  resize(sz);
  for(k = old_size; k < sz; k++) *cont[k] = c;
}

template<class T>
inline int PtrCont<T>::capacity() const { return cont.size()+junk.size(); }

template<class T>
inline bool PtrCont<T>::empty() const { return cont.size() == 0; }

template<class T>
void PtrCont<T>::reserve(int n) {
  int k, old_capacity = capacity();
  if(old_capacity < n) {
    junk.resize(n);
    for(k = old_capacity; k < n; k++) junk[k] = new T();
  }
}

template<class T>
inline void PtrCont<T>::shrink(double p) {
  int n = junk.size();
  if(n) {
    int len = int(p*n+0.5);
    for(int k = n-len; k < n; k++) delete junk[k];
    junk.erase(junk.end()-len, junk.end());
  }
}

template<class T>
inline T& PtrCont<T>::operator[] (int n) const { return *cont[n];}

template<class T>
inline T& PtrCont<T>::front() const { return *cont.front(); }

template<class T>
inline T& PtrCont<T>::back() const { return *cont.back(); }

template<class T>
void PtrCont<T>::push_back() {
  if(junk.size() == 0) cont.push_back(new T());
  else {
    cont.push_back(junk.back());
    junk.pop_back();
  }
  cont.back()->id.identity = cont.size()-1;
}

template<class T>
void PtrCont<T>::push_back(const T& x) {
  push_back();
  *cont.back() = x;
}

template<class T>
void PtrCont<T>::pop_back() {
  if(!empty()) {
    junk.push_back(cont.back());
    *junk.back() = T();
    cont.pop_back();
  }
  auto_shrink();
}

template<class T>
inline void PtrCont<T>::clear() { resize(0); }

template<class T>
template<class Compare>
inline void PtrCont<T>::sort(const Compare& compare) {
  std::sort(cont.begin(), cont.begin()+size(), 
	    PtrCompare<T, Compare>(compare));
}

template<class T>
inline void PtrCont<T>::exch(int i, int j) { 
  std::swap(cont[i], cont[j]); 
  std::swap(cont[i]->id.identity, cont[j]->id.identity);
}

template<class T>
inline T& PtrCont<T>::insert() { push_back(); return *cont.back(); }

template<class T>
inline T& PtrCont<T>::insert(const T& x) { push_back(x); return *cont.back(); }

template<class T>
T& PtrCont<T>::insert(int pos) {
  push_back();
  exch(pos, back().id.identity);
  return *cont[pos];
}

template<class T>
T& PtrCont<T>::insert(int pos, const T& x) {
  push_back(x);
  exch(pos, back().id.identity);
  return *cont[pos];
}

template<class T>
void PtrCont<T>::erase(int pos) {
  exch(pos, back().id.identity);
  pop_back();
}

template<class T>
template<class InputIterator>
void PtrCont<T>::keep(InputIterator first, InputIterator last) {
  int k, sz = size();
  std::vector<T*> labels = cont;
  while(first != last) {
    if(int(*first) >= 0 && int(*first) < int(sz))
      labels[*first] = (T*)NULL;
    first++;
  }
  for(k = 0; k < sz; k++) 
    if(labels[k] != (T*)NULL) erase(labels[k]->id.identity);
}

template<class T>
template<class InputIterator>
void PtrCont<T>::erase(InputIterator first, InputIterator last) {
  int k, sz = size();
  std::vector<T*> labels(sz, (T*)NULL);
  while(first != last) {
    if(int(*first) >= 0 && int(*first) < int(sz))
      labels[*first] = cont[*first];
    first++;
  }
  for(k = 0; k < sz; k++) 
    if(labels[k] != (T*)NULL) erase(labels[k]->id.identity);
}

template<class T>
inline void PtrCont<T>::auto_shrink() {
  if(junk.size() > cont.size()) shrink(0.8);
}

}/* Util */

#endif
