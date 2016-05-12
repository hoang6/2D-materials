#include "MC.h"
#include <ctime>
#include <cmath>

namespace MC {

void useNVT(long seed, 
	    AMod::Molecule& mol,
	    MCBasic& mc,
	    MCAgent& mcAgent,
	    double kT,
	    double dxmax, 
	    long rStep, 
	    long eStep, 
	    Stream& stream) {
  long kStep;
  long kkStep;
  long nStep;
  long natoms;
  double volume;
  long atomID;
  double new_x[3];
  double dE;
  double oE;
  double dEt;
  bool accept;
  long acc_count;
  long tot_count;
  double prev_tot_count;
  double acc_ratio;
  double tic;
  double toc;
  double ttic;
  double ttoc;
  bool reset_flag;
  int reset_n;

  /********** INITIALIZE **********/
  ttic = time(NULL);
  Util::RNG::seed(seed);
  mcAgent.attach(mol);
  mc.attach(mol);
  mcAgent.init();
  mc.init();
  oE = mol.totalEnergy();
  dEt = 0.0;
  natoms = mol.natoms();
  volume = mol.volume();
  acc_count = 0;
  tot_count = 0;
  prev_tot_count = 0.0;
  acc_ratio = 0.0;
  nStep = rStep+eStep;
  reset_flag = true;
  reset_n = 0;

  /********** HEADER **********/
  if(stream.showHeader()) {
    stream.info() << "MCStep" << " ";
    stream.info() << "Potential" << " ";
    stream.info() << "Acc_Ratio" << " ";
    stream.info() << "dxmax" << " ";
    stream.info() << "time" << std::endl;
    //
    stream.stat() << "N" << " " << std::noshowpos << natoms << std::endl;
    stream.stat() << "V" << " " << std::noshowpos << volume << std::endl;
    stream.stat() << "kT" << " " << std::noshowpos << kT << std::endl; 
    stream.stat() << "Potential" << std::endl;
  }
  
  /********** MC LOOP **********/
  tic = time(NULL);
  for(kStep = 0; kStep < nStep; kStep++) {
    for(kkStep = 0; kkStep < natoms; kkStep++) {
      atomID = Util::RNG::uniform_SFMT()*natoms;
      //random move
      mcAgent.pertAtom(dxmax,atomID,new_x);
      dE = mc.dEnergyMoveAtom(atomID,new_x);
      //Metropolis MC
      accept = (dE <= 0 || std::log(1.0-Util::RNG::uniform_SFMT()) <= -dE/kT);
      if(accept) {//accept
	dEt += dE;
	acc_count++;
      }
      else {//reject: undo move
	mc.unMoveAtom();
      }
      mcAgent.accPertAtom(accept);
      tot_count++;
      //record
      if(stream.dAction()) {
	if(accept) {
	  stream.action() << std::noshowpos << atomID << " ";
	  Util::vdump(stream.action(),mol.atom(atomID)->x,3) << std::endl;  
	}
	else {
	  stream.action() << "r" << std::endl;
	}
      }
      /********** ADJUST dxmax **********/
      if(tot_count == 100) {
	acc_ratio = double(acc_count)+acc_ratio*prev_tot_count;
	acc_ratio /= double(tot_count)+prev_tot_count;
	prev_tot_count += tot_count;
	acc_count = 0;
	tot_count = 0;
	if(reset_flag)
	  dxmax *= (1.0+acc_ratio-0.5);
      }
    }/* for(kkStep... */
    if(reset_flag) {
      reset_n++;
      if(reset_n >= rStep)
	reset_flag = false;
    }
    //update oE
    oE += dEt;
    dEt = 0.0;
    //dump
    if((kStep+1)%stream.dInfoStep() == 0) {
      toc = time(NULL);
      stream.info() << std::noshowpos << kStep+1 << " ";
      stream.info() << std::showpos << oE << " ";
      stream.info() << std::noshowpos << acc_ratio << " ";
      stream.info() << std::noshowpos << dxmax << " ";
      stream.info() << std::noshowpos << toc-tic << std::endl;
      tic = time(NULL);
    }//end if((kStep+1)%stream.dInfoStep() == 0)
    //dump stat header
    //dump stat
    if((kStep+1)%stream.dStatStep() == 0) {
      stream.stat() << std::showpos << oE << " ";
    }
    //dump state
    if((kStep+1)%stream.dStateStep() == 0) {
      mol.io().dumpTxt(stream.state(false));
      //add a tag (blank line) in each streams
      stream.action(false) << std::endl;
      stream.state(false) << std::endl;
      stream.stat(false) << std::endl;
      stream.info(false) << std::endl;
    }
  }/* for(kStep... */
  
  /********** FINALIZE *********/
  mc.final();
  mcAgent.final();
  ttoc = time(NULL);
  std::cout << "[Total time] " 
	    << std::noshowpos << ttoc-ttic << " seconds" << std::endl;
}

}/* MC */
