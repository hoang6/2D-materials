#include "../Util/util.h"
#include "../Util/RNG.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>

int main(int argc, char* argv[]) {
  using std::cout;
  using std::endl;
  if(argc != 4) {
    cout << "Usage: testPreClockMC <n_atoms> <n_MC_moves> <dxmax>" << endl;
    return 1;
  }

  int natoms;
  if(!Util::readWord(argv[1],natoms)) {
    cout << "Error: cannot read <natoms>" << endl;
    return 2;
  }
  int nMCMoves;
  if(!Util::readWord(argv[2],nMCMoves)) {
    cout << "Error: cannot read <n_MC_moves>" << endl;
    return 3;
  }
  double dxmax;
  if(!Util::readWord(argv[3],dxmax)) {
    cout << "Error: cannot read <dxmax>" << endl;
    return 4;
  }

  Util::RNG::seed(std::time(NULL));
  int i, j, atomID;
  double dx[3];
  bool accept;
  const int digits = 12;
  cout.setf(std::ios::scientific, std::ios::floatfield);
  cout.precision(digits);
  cout << std::left;
  for(i = 0; i < nMCMoves; i++) {
    atomID = Util::RNG::uniform_SFMT()*natoms;
    for(j = 0; j < 3; j++)
      dx[j] = dxmax*(-1.0+2.0*Util::RNG::uniform_SFMT());
    accept = (Util::RNG::uniform_SFMT()<0.5?true:false);
    cout << std::setw(digits+2) << std::noshowpos << atomID << " "
         << std::setw(digits+9) << std::showpos << dx[0] << " "
         << std::setw(digits+9) << std::showpos << dx[1] << " "
         << std::setw(digits+9) << std::showpos << dx[2] << " "
         << std::noshowpos << accept << endl;
  }

  return 0;
}
