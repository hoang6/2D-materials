#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#define ERROR 1e-16

using namespace std;

static int IN3[64+1][3+1];
static double CLMN[3+1][10+1][10+1][10+1][64+1];
static double TLMN[10+1][10+1][10+1][64+1];

static bool init_in3() {
  int IC=0;
  int I3D,ITD;
  int i, I, J, K, L, M, N;
  double val;
  FILE *inter3Dtors;
  for (I=1; I<5; I++) {
    for (J=1; J<5; J++) {
      for (K=1; K<5; K++) {
        IC=IC+1;
        /* It's weird that IN3 and dIN3 are the same as each other. */
        IN3[IC][1]=I-1;
        IN3[IC][2]=J-1;
        IN3[IC][3]=K-1;
      }
    }
  }
  /* Read the tricubic spline coefficients for RADIC (that is, F). */
  for(i = 0; i < 3; ++i) {
    static const char *inter_file[]={"inter3d_iv_new.d",
				     "inter3d_ch.d",
				     "inter3d_h.d"};
    FILE *inter3Divnew;
    if((inter3Divnew = fopen(inter_file[i],"r"))==NULL) {
      cerr << "opening " << inter_file[i] << " failed" << endl;;
      return 0;
    }
    fscanf(inter3Divnew, "%d\n",&I3D);
    while(1) {
      if(fscanf(inter3Divnew, "%d %d %d\n",&L,&M,&N) != 3) break;
      for(I=1; I<=64; ++I) {
	double val;
	if(!fscanf(inter3Divnew,"%lf\n", &val)) {
	  cerr << "parse error in file " << inter_file[i] << ", "
	       << L << " " << M << " " << N << " " << I << endl;
	  return 0;
	}
        CLMN[i+1][L][M][N][I] = val;
      }
    }
    fclose(inter3Divnew);
  }   
  for (L=1; L<11; L++) {
    for (M=1; M<11; M++) {
      for (I=1; I<65; I++) {
	CLMN[1][L][M][10][I]=CLMN[1][L][M][9][I];
	CLMN[2][L][M][10][I]=CLMN[2][L][M][9][I];
	for (N=6; N<11; N++)
	CLMN[3][L][M][N][I]=CLMN[3][L][M][5][I];
      }
    }
  }

  // Read tricubic spline coefficients for torsional potential
  if ((inter3Dtors = fopen("inter3dtors.d","r"))==NULL)
    cerr << "opening inter3dtors failed" << endl;
  fscanf(inter3Dtors, "%d\n",&ITD);
  for (J=1; J<109; J++) {
    fscanf(inter3Dtors, "%d %d %d\n",&L,&M,&N);
    for (I=1; I<64; I+=3) {
      double i0, i1, i2;
      fscanf(inter3Dtors, "%lf %lf %lf\n", &i0, &i1, &i2);
      TLMN[L][M][N][I] = i0;
      TLMN[L][M][N][I+1] = i1;
      TLMN[L][M][N][I+2] = i2;
    }
    fscanf(inter3Dtors, "%lf\n",&val);
    TLMN[L][M][N][64] = val;
  }
  fclose(inter3Dtors);
  for (L=1; L<11; L++)
    for (M=1; M<11; M++)
      for (N=4; N<11; N++)
	for (I=1; I<65; I++) {
	  TLMN[L][M][N][I]=TLMN[L][M][3][I];
	}
  return 1;
}

int main() {
  if(!init_in3()) exit(1);
  //--------------- make in3.h ---------------
  ofstream fout("../in3.h");
  fout.precision(16);
  fout << "#ifndef _TBTOOLS_IN3_H" << endl;
  fout << "#define _TBTOOLS_IN3_H" << endl;
  fout << endl;
  fout << "namespace tbtools {" << endl;
  fout << endl;
  //define global arrays: IN3, CLMN, TLMN
  fout << "static int IN3[64+1][3+1];" << endl;
  fout << "static double CLMN[3+1][10+1][10+1][10+1][64+1];" << endl;
  fout << "static double TLMN[10+1][10+1][10+1][64+1];" << endl;
  fout << "static int in3_initialized = 0;" << endl;
  fout << endl;
  //define function init_in3()
  int i, n;
  fout << "void init_in3() {" << endl;
  //IN3_VAL
  n = (64+1)*(3+1);
  fout << "int IN3_VAL[" << n << "]  = {";
  for(i = 0; i < n; i++) {
    fout << *(&IN3[0][0]+i);
    if(i != n-1)
      fout << ", ";
  }
  fout << "};" << endl;
  //compress CLMN
  vector<int> CLMN_IND;
  n = (3+1)*(10+1)*(10+1)*(10+1)*(64+1);
  for(i = 0; i < n; i++) {
    if(fabs(*(&CLMN[0][0][0][0][0]+i)) > ERROR) {
      CLMN_IND.push_back(i);
    }
  }
  //CLMN_IND
  n = CLMN_IND.size();
  fout << "int CLMN_IND[" << n << "] = {";
  for(i = 0; i < n; i++) {
    fout << CLMN_IND[i];
    if(i != n-1) fout << ", ";
  }
  fout << "}; " << endl;
  //CLMN_VAL
  fout << "double CLMN_VAL[" << n << "] = {";
  for(i = 0; i < n; i++) {
    fout << *(&CLMN[0][0][0][0][0]+CLMN_IND[i]);
    if(i != n-1) fout << ", ";
  }
  fout << "}; " << endl;
  //compress TLMN
  vector<int> TLMN_IND;
  n = (10+1)*(10+1)*(10+1)*(64+1);
  for(i = 0; i < n; i++) {
    if(fabs(*(&TLMN[0][0][0][0]+i)) > ERROR) {
      TLMN_IND.push_back(i);
    }
  }
  //TLMN_IND
  n = TLMN_IND.size();
  fout << "int TLMN_IND[" << n << "] = {";
  for(i = 0; i < n; i++) {
    fout << TLMN_IND[i];
    if(i != n-1) fout << ", ";
  }
  fout << "}; " << endl;
  //TLMN_VAL
  fout << "double TLMN_VAL[" << n << "] = {";
  for(i = 0; i < n; i++) {
    fout << *(&TLMN[0][0][0][0]+TLMN_IND[i]);
    if(i != n-1) fout << ", ";
  }
  fout << "}; " << endl;
  fout << endl;
  //global variables assignment
  fout << "int* IN3_ptr = &IN3[0][0];" << endl;
  fout << "double* CLMN_ptr = &CLMN[0][0][0][0][0];" << endl;
  fout << "double* TLMN_ptr = &TLMN[0][0][0][0];" << endl;
  n = (64+1)*(3+1);
  fout << "for(int i = 0; i < " << n << "; i++) {" << endl;
  fout << "*(IN3_ptr+i) = IN3_VAL[i];" << endl;
  fout << "}" << endl;
  n = CLMN_IND.size();
  fout << "for(int i = 0; i < " << n << "; i++) {" << endl;
  fout << "*(CLMN_ptr+CLMN_IND[i]) = CLMN_VAL[i];" << endl;
  fout << "}" << endl;
  n = TLMN_IND.size();
  fout << "for(int i = 0; i < " << n << "; i++) {" << endl;
  fout << "*(TLMN_ptr+TLMN_IND[i]) = TLMN_VAL[i];" << endl;
  fout << "}" << endl;
  fout << "in3_initialized = 1;" << endl;
  fout << "}" << endl;
  fout << endl;
  fout << "}" << endl;
  fout << endl;
  fout << "#endif" << endl;
  fout.close();
  return 0;
}

#undef ERROR
