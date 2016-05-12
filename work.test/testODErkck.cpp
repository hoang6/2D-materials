#include "../Util/ODErkck.h"
#include "../Util/util.h"
#include <iostream>
#include <cmath>
using namespace std;

typedef Util::ODErkck<8> ODErkck;

class Orbit {
public:
  Orbit(double _m1, double _m2, double _m3): m1(_m1), m2(_m2), m3(_m3) {}
  void d(double _t, const ODErkck::Data& _u, ODErkck::Data& _du) {
    double R[2] = {_u[0],_u[1]};
    double r[2] = {_u[2],_u[3]};
    double Rr[2] = {R[0]-r[0],R[1]-r[1]};
    double nR3 = std::pow(Util::norm2(R,2),3);
    double nr3 = std::pow(Util::norm2(r,2),3);
    double nRr3 = std::pow(Util::norm2(Rr,2),3);
    for(int i = 0; i < 4; i++) _du[i] = _u[i+4];
    for(int i = 0; i < 2; i++)
      _du[i+4] = -m1*R[i]/nR3-m3*Rr[i]/nRr3;
    for(int i = 0; i < 2; i++)
      _du[i+6] = -m1*r[i]/nr3+m2*Rr[i]/nRr3;
  }
  bool stop(double _t, const ODErkck::Data& _u) { 
    return false;
  }

private:
  double m1, m2, m3;
};

int main() {
  //setup solar
  double m1 = 1.0;
  double m2 = 1.0/100;
  double m3 = 1.0/8000;
  double T = 30.0;
  double u0[8] = {5.0, 0.0, 5.1, 0.0, 0.0, 0.2, 0.0, 0.3};
  Orbit orbit(m1,m2,m3);
  ODErkck::Soln s;

  //solve ODE
  ODErkck().solve(orbit, 0.0, T, ODErkck::Data(u0), 1.0e-6, 1.0e-6, 1.0, s);

  //output result
  for(int k = 0, kmax = s.t.size(); k < kmax; k++) {
    cout << s.t[k] << "\t"; 
    for(int i = 0; i < 8; i++) cout << s.u[k][i] << "\t";
    cout << endl;
  }
}
