#include "../MD/LJ.h"
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
  double e = 0.0014;
  double s = 2.74;
  double c = 2.50;

  MD::LJ lj(e,s);
  MD::LJCutoff ljc(e,s,c);

  double r;
  int i, n = 1000;
  for(i = 0; i < n; i++) {
    r = (0.88+i*(1.8-0.88)/n)*s;
    cout << r << " " << lj.energy(r) << " " << ljc.energy(r) << endl;
  }
  cerr << ljc.energy(c*s) << " " << ljc.denergy(c*s) << endl;
  cerr << ljc.energy(c*s+1.0e-15) << " " << ljc.denergy(c*s+1.0e-15) << endl;

  return 0;
}
