#include "MolOptim.h"
#include "../Util/RNG.h"
#include "../Util/util.h"
#include "KernelREBO.h"
#include "KernelLCBOPIIN.h"
#include <ctime>

static const AMod::MolOptim* G_MOL_OPTIM;
static const AMod::MolOptim::DonlpPara* G_DONLP_PARA;
static double* G_WS;
extern "C" void donlp2();

namespace AMod {

MolOptim::MolOptim() {}

MolOptim::MolOptim(Molecule& _mol): MolServ(_mol) {}

MolOptim::MolOptim(const MolOptim& molOptim) { *this = molOptim; }

MolOptim& MolOptim::operator= (const MolOptim& molOptim) {
  if(this != &molOptim) {
    MolServ::operator= (molOptim);
  }
  return *this; 
}

void MolOptim::init(PotType _potType) {
  if(_potType == REBO)
    _kernelPtr = new KernelREBO(molecule());
  else if(_potType == LCBOPIIN)
    _kernelPtr = new KernelLCBOPIIN(molecule());
  _kernelPtr->init();
}

void MolOptim::final() {
  _kernelPtr->final();
  delete _kernelPtr;
}

void MolOptim::optimize(const DonlpPara& _para) {
  G_MOL_OPTIM = &*this;
  G_DONLP_PARA = &_para;
  G_WS = new double[nDOF()+1];

  donlp2();

  delete [] G_WS;
}

int MolOptim::nDOF() const {
  return 3*(molecule().naxes()+molecule().natoms());
}

void MolOptim::set_x(const double* _x) const {
  Molecule& mol = molecule();
  int j, k, c = 0, naxes = mol.naxes(), natoms = mol.natoms();
  for(k = 0; k < naxes; k++)
    for(j = 0; j < 3; j++) mol.axis(k).x[j] = _x[c++];
  for(k = 0; k < natoms; k++)
    for(j = 0; j < 3; j++) mol.atom(k)->x[j] = _x[c++];
  _kernelPtr->update();
}

void MolOptim::get_x(double* _x) const {
  const Molecule& mol = molecule();
  int j, k, c = 0, naxes = mol.naxes(), natoms = mol.natoms();
  for(k = 0; k < naxes; k++)
    for(j = 0; j < 3; j++) _x[c++] = mol.axis(k).x[j];
  for(k = 0; k < natoms; k++) 
    for(j = 0; j < 3; j++) _x[c++] = mol.atom(k)->x[j];
}

void MolOptim::perturb_x(double* _x, double _amp) const {
  const Molecule& mol = molecule();
  int j, k, c = 0, naxes = mol.naxes(), natoms = mol.natoms();
  for(k = 0; k < naxes; k++) 
    for(j = 0; j < 3; j++) {
      if(!mol.axis(k).fixed[j])
	_x[c] += _amp*(2.0*Util::RNG::uniform_SFMT()-1.0);
      c++;
    }
  for(k = 0; k < natoms; k++) 
    for(j = 0; j < 3; j++) {
      if(!mol.atom(k)->fixed[j])
	_x[c] += _amp*(2.0*Util::RNG::uniform_SFMT()-1.0);
      c++;
    }
}

void MolOptim::get_g(double* _g) const {
  const Molecule& mol = molecule();
  int j, k, c = 0, naxes = mol.naxes(), natoms = mol.natoms();
  for(k = 0; k < naxes; k++)
    for(j = 0; j < 3; j++)
      _g[c++] = (mol.axis(k).fixed[j]?0.0:-mol.axis(k).ftotal(j));
  for(k = 0; k < natoms; k++)
    for(j = 0; j < 3; j++)
      _g[c++] = (mol.atom(k)->fixed[j]?0.0:-mol.atom(k)->ftotal(j));
}

double MolOptim::solver(double* _x, double* _g) const {
  bool fxOnly = (_g == NULL);
  double fx;

  set_x(_x);

  fx = _kernelPtr->energy();
  if(fxOnly) return fx;

  _kernelPtr->denergy();
  get_g(_g);

  return fx;
}

}/* AMod */

/****************************** donlp2 ******************************/
extern "C" {
  #include "../donlp2/o8para.h"

  void user_init_size() {
    #define  X extern
    #include "../donlp2/o8comm.h"
    #include "../donlp2/o8fint.h"
    #include "../donlp2/o8cons.h"
    #undef   X

    n = G_MOL_OPTIM->nDOF();
    nlin = 0;
    nonlin = 0;
    iterma = 10000;
    nstep = 20;
  }

  void user_init() {
    #define  X extern
    #include "../donlp2/o8comm.h"
    #include "../donlp2/o8fint.h"
    #include "../donlp2/o8cons.h"
    #undef   X
    
    Util::RNG::seed(std::time(NULL));

    G_MOL_OPTIM->get_x(x+1);
    G_MOL_OPTIM->perturb_x(x+1, G_DONLP_PARA->amp);

    strcpy(name, G_DONLP_PARA->tag.c_str());
    silent = G_DONLP_PARA->silent;
    bloc = TRUE;
    tau0 = 1.0;
    del0 = 0.05;    
    nreset = n;
    analyt = TRUE;
  }

  void eval_extern(INTEGER mode) {
    #define  X extern
    #include "../donlp2/o8comm.h"
    #include "../donlp2/o8fint.h"
    #undef   X
    #include "../donlp2/o8cons.h"

    int i;

    //compute the obj function & the nonlinear constraints
    if(mode == 1) {
      fu[0] = G_MOL_OPTIM->solver(xtr+1, NULL);
    }
    //compute the obj function & the nonlinear constraints & their gradients
    else if(mode == 2) {
      fu[0] = G_MOL_OPTIM->solver(xtr+1, G_WS+1);
      for(i = 1; i <= n; i++)
	fugrad[i][0] = G_WS[i];
    }
  }

  static double G_LAST_FX;
  static int G_COUNT;

  void newx (DOUBLE x[], DOUBLE u[], INTEGER itstep, 
	     DOUBLE **accinf, LOGICAL *cont) {
    double fx = accinf[itstep][2];
    if(itstep >= 20 && std::abs(fx-G_LAST_FX) < 1.0e-6) {
      G_COUNT++;
      *cont = (G_COUNT < 10);
    }
    else {
      G_LAST_FX = fx;
      G_COUNT = 0;
      *cont = 1;
    }
  }

  void setup() {
    #define  X extern
    #include "../donlp2/o8comm.h"
    #undef   X

    te0 = TRUE;
  }

  void solchk() {
    /*
    #define  X extern
    #include "../donlp2/o8comm.h"
    #undef   X
    #include "../donlp2/o8cons.h"
    */
  }

  void ef(DOUBLE x[],DOUBLE *fx) {
    /*
    #define  X extern
    #include "../donlp2/o8fuco.h"
    #undef   X
    */
  }

  void egradf(DOUBLE x[],DOUBLE gradf[]) {
    /*
    #define  X extern
    #include "../donlp2/o8fuco.h"
    #undef   X
    */
  }

  void econ(INTEGER type,INTEGER liste[],DOUBLE x[],DOUBLE con[],LOGICAL err[]) {
    /*
    #define  X extern
    #include "../donlp2/o8fuco.h"
    #undef   X
    */
  }  

  void econgrad(INTEGER liste[],INTEGER shift,DOUBLE x[],DOUBLE **grad) {
    /*
    #define  X extern
    #include "../donlp2/o8fuco.h"
    #undef   X
    */
  }

}/* extern "C" */
