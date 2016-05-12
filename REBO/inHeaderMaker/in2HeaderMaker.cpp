#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

#define ERROR 1e-16
#define MAX_BC_NEIGHBORS 32

static int IN2[16+1][2+1];
static double CLM[2+1][MAX_BC_NEIGHBORS+1][MAX_BC_NEIGHBORS+1][16+1];
static double xh[3][5][5];

bool init_in2() {
  int IC=0;
  int I, J, K, L, M;
  double xhh;
  double I2D;
  FILE *inter2D;
  inter2D = fopen("inter2d_iv.d", "r");
  if(inter2D == NULL) {
    cerr << "Unable to open inter2d_iv.d\n" << endl;
    return 0;
  }
  fscanf(inter2D, "%lf\n",&I2D); //discard the first number
  //--------------initialize xh-----------
  do {
    fscanf(inter2D, "%d %d %d %lf\n",&I,&J,&K,&xhh);
    if (I > 0) xh[I][J][K] = xhh;
  }while(I > 0);
  //---------------initialize CLM---------
  for(I=1; I<5; I++) {
    for(J=1; J<5; J++) {
      IC=IC+1;
      IN2[IC][1]=I-1;
      IN2[IC][2]=J-1;
    }
  }
  //zero the coefficients
  for(I=1; I<3; I++)
    for(L=1; L <= 32; L++)
      for(M=1; M <= 32; M++)
	for(J=1; J<17; J++)
	  CLM[I][L][M][J]=0.0;
  for(K=1; K<73; K++) {
    fscanf(inter2D, "%d %d %d\n",&I,&L,&M);
    for(J=1; J<16; J+=4) {
      double j0, j1, j2, j3;
      fscanf(inter2D, "%lf %lf %lf %lf\n", &j0, &j1, &j2, &j3);
      CLM[I][L][M][J] = (double) j0;
      CLM[I][L][M][J+1] = (double) j1;
      CLM[I][L][M][J+2] = (double) j2;
      CLM[I][L][M][J+3] = (double) j3;
    }
  }
  return 1;
}

int main() {
  if(!init_in2()) exit(1);
  //--------------- make in2.h ---------------
  ofstream fout("../in2.h");
  fout.precision(16);
  fout << "#ifndef _TBTOOLS_IN2_H" << endl;
  fout << "#define _TBTOOLS_IN2_H" << endl;
  fout << endl;
  fout << "namespace tbtools {" << endl;
  fout << endl;
  fout << "#define MAX_BC_NEIGHBORS " << MAX_BC_NEIGHBORS << endl;
  fout << endl;
  //define global arrays: IN2, xh, CLM
  fout << "static int IN2[16+1][2+1];" << endl;
  fout << "static double xh[3][5][5];" << endl;
  fout << "static double CLM[2+1][MAX_BC_NEIGHBORS+1][MAX_BC_NEIGHBORS+1][16+1];" << endl;
  fout << "static int in2_initialized = 0;" << endl;
  fout << endl;
  //define function init_in2()
  int i, n;
  fout << "void init_in2() {" << endl;
  //IN2_VAL
  n = (16+1)*(2+1);
  fout << "double IN2_VAL[" << n << "] = {";
  for(i = 0; i < n; i++) {
    fout << *(&IN2[0][0]+i);
    if(i != n-1) fout << ", ";
  }
  fout << "};" << endl;
  //xh_VAL
  n = 3*5*5;
  fout << "double xh_VAL[" << n << "] = {";
  for(i = 0; i < n; i++) {
    fout << *(&xh[0][0][0]+i);
    if(i != n-1) fout << ", ";
  }
  fout << "};" << endl;
  //compress CLM
  n = (2+1)*(MAX_BC_NEIGHBORS+1)*(MAX_BC_NEIGHBORS+1)*(16+1);
  vector<int> CLM_IND;
  for(i = 0; i < n; i++) {
    if(fabs(*(&CLM[0][0][0][0]+i)) > ERROR) {
      CLM_IND.push_back(i);
    }
  }
  //CLM_IND
  n = CLM_IND.size();
  fout << "int CLM_IND[" << n << "] = {";
  for(i = 0; i < n; i++) {
    fout << CLM_IND[i];
    if(i != n-1) fout << ", ";
  }
  fout << "}; " << endl;
  //CLM_VAL
  fout << "double CLM_VAL[" << n << "] = {";
  for(i = 0; i < n; i++) {
    fout << *(&CLM[0][0][0][0]+CLM_IND[i]);
    if(i != n-1) fout << ", ";
  }
  fout << "}; " << endl;
  fout <<endl;
  //global variables assignment
  fout << "int* IN2_ptr = &IN2[0][0];" << endl;
  fout << "double* xh_ptr = &xh[0][0][0];" << endl;
  fout << "double* CLM_ptr = &CLM[0][0][0][0];" << endl;  
  n = (16+1)*(2+1);
  fout << "for(int i = 0; i < " << n << "; i++) {" << endl;
  fout << "*(IN2_ptr+i) = IN2_VAL[i];" << endl;
  fout << "}" << endl;
  n = 3*5*5;
  fout << "for(int i = 0; i < " << n << "; i++) {" << endl;
  fout << "*(xh_ptr+i) = xh_VAL[i];" << endl;
  fout << "}" << endl;
  n = CLM_IND.size();
  fout << "for(int i = 0; i < " << n << "; i++) {" << endl;
  fout << "*(CLM_ptr+CLM_IND[i]) = CLM_VAL[i];" << endl;
  fout << "}" << endl;
  fout << endl;
  fout << "in2_initialized = 1;" << endl;
  fout << "}" << endl;
  fout << endl;
  fout << "}" << endl;
  fout << endl;
  fout << "#endif" << endl;
  fout.close();
  return 0;
}

#undef MAX_BC_NEIGHBORS
#undef ERROR
