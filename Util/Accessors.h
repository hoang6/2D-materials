/*
  The class Accessors simplifies setter and getter, copy from
  http:://www.kirit.com
  keywords: a simple meta accessor
 */

#ifndef _UTIL_ACCESSORS_H
#define _UTIL_ACCESSORS_H

namespace Util {

template<class T>
class Accessors {
public:
  Accessors() {}
  explicit Accessors(const T& _t): t(_t) {}
  T& operator()() { return t; }
  const T& operator()() const { return t; }
  void operator()(const T& _t) { t = _t; }
private:
  T t;
}; 

} /* Util */

#endif
