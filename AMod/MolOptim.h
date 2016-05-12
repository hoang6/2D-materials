#ifndef _AMOD_MOL_OPTIM_H
#define _AMOD_MOL_OPTIM_H

#include "MolService.h"
#include "KernelBasic.h"

namespace AMod {

class MolOptim: public MolServ {
public:
  struct DonlpPara;

  MolOptim();
  MolOptim(Molecule& _mol);
  MolOptim(const MolOptim& molOptim);
  MolOptim& operator= (const MolOptim& molOptim);
  void init(PotType _potType);
  void final();
  void optimize(const DonlpPara& _para);
  /********** called by donlp2 **********/
  int nDOF() const;
  void set_x(const double* _x) const;
  void get_x(double* _x) const;
  void perturb_x(double* _x, double _amp) const;
  void get_g(double* _g) const;
  double solver(double* _x, double* _g) const; //_x and _g start from 0

  struct DonlpPara {
    std::string tag;
    bool silent;
    double amp;
    //----------
    DonlpPara() { reset(); }
    DonlpPara(const std::string& _tag, bool _silent, double _amp):
      tag(_tag), silent(_silent), amp(_amp) {}
    void reset() { tag = ""; silent = false; amp = 0.0; }
  };

private:
  KernelBasic* _kernelPtr;
};

}/* AMod */

#endif
