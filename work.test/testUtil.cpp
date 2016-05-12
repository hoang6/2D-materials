#include "../Util/util.h"
#include <iostream>
using namespace std;

void mdump(const double* m, int r, int s);

int main() {
  double A[][4] = {
    {1,10,2,8},
    {5,9,11,12},
    {8.0,4.6,1.5,3.8},
    {10.0,2.34,7.56,1.2}
  };
  double B[][4] = {
    {5,2,3,0},
    {4,9,3,8},
    {6.5,6.89,10.2,1.8},
    {5.6,7.8,9.12,12.2}
  };

  double mA[16];
  double mB[16];
  double mX[16];
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      mA[i*4+j] = A[i][j];
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      mB[i*4+j] = B[i][j];
  
  cout << "Matrix A:" << endl;
  mdump(mA, 4, 4);
  cout << endl;

  cout << "Matrix B:" << endl;
  mdump(mB, 4, 4);
  cout << endl;

  cout << "Matrix X: (X = A)" << endl;
  Util::mcopy(mA, mX, 4, 4);
  mdump(mX, 4, 4);
  cout << endl;

  cout << "Matrix X: (X = A')" << endl;
  Util::mtrans(mA, mX, 4, 4);
  mdump(mX, 4, 4);
  cout << endl; 

  cout << "Matrix X: (X = A+B)" << endl;
  Util::madd(mA, mB, mX, 4, 4);
  mdump(mX, 4, 4);
  cout << endl;

  cout << "Matrix X: (X = A-B)" << endl;
  Util::msub(mA, mB, mX, 4, 4);
  mdump(mX, 4, 4);
  cout << endl;

  cout << "Matrix X: (X = 0.2*A)" << endl;
  Util::kdotm(0.2, mA, mX, 4, 4);
  mdump(mX, 4, 4);
  cout << endl;

  cout << "Vector x: (x = A*[3 5 9 8]')" << endl;
  double x[4], b[4] = {3,5,9,8};
  Util::mdotv(mA, b, x, 4, 4);
  for(int k = 0; k < 4; k++) cout << x[k] << " ";
  cout << endl;
  cout << endl;

  cout << "Matrix X: (X = A*B)" << endl;
  Util::mdotm(mA, mB, mX, 4, 4, 4);
  mdump(mX, 4, 4);
  cout << endl;

  cout << "Trance(A): " << Util::mtrace(mA, 4) << endl << endl;

  cout << "det(A) = " << Util::mdet(mA, 4) << endl << endl;

  cout << "cross_dot3(A0,A1,A2) = "
       << Util::cross_dot3(mA,mA+4,mA+8) << endl << endl;

  cout << "Matrix X: (X = inv(A))" << endl;
  Util::minverse(mA, mX, 4);
  mdump(mX, 4, 4);
  cout << endl;

  return 0;
}

void mdump(const double* m, int r, int s) {
  for(int i = 0; i < r; i++) {
    for(int j = 0; j < s; j++) { 
      cout << m[i*s+j] << "\t";
    }
    cout << endl;
  }
}
