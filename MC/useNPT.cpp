#include "MC.h"
#include <ctime>
#include <cmath>

namespace MC {

void useNPT(long seed, 
	    AMod::Molecule& mol,
	    MCBasic& mc,
	    MCAgent& mcAgent,
	    double pressure,
	    double kT,
	    double dxmax, 
	    double dhmax,
	    long rStep, 
	    long eStep, 
	    long hStep, 
	    Stream& stream) {
  long kStep;
  long kkStep;
  long nStep;
  long natoms;
  long atomID;
  double new_x[3];
  double old_v;
  double new_v;
  double dW;
  double dE;
  double oE;
  double dEt;
  bool accept;
  long acc_count;
  long tot_count;
  double prev_tot_count;
  double acc_ratio;
  long h_acc_count;
  long h_tot_count;
  double prev_h_tot_count;
  double h_acc_ratio;
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
  acc_count = 0;
  tot_count = 0;
  prev_tot_count = 0.0;
  acc_ratio = 0.0;
  h_acc_count = 0;
  h_tot_count = 0;
  prev_h_tot_count = 0.0;
  h_acc_ratio = 0.0;
  nStep = rStep+eStep;
  old_v = mol.volume();
  reset_flag = true;
  reset_n = 0;

  /********** HEADER **********/
  if(stream.showHeader()) {
    stream.info() << "MCStep" << " ";
    stream.info() << "Potential" << " ";
    stream.info() << "Acc_Ratio" << " ";
    stream.info() << "Cell_Acc_Ratio" << " ";
    stream.info() << "dxmax" << " ";
    stream.info() << "dhmax" << " ";
    stream.info() << "time" << std::endl;
    //
    stream.stat() << "N" << " " << std::noshowpos << natoms << std::endl;
    stream.stat() << "P" << " " << std::showpos << pressure << std::endl;
    stream.stat() << "kT" << " " << std::noshowpos << kT << std::endl;
    stream.stat() << "Potential" << " ";
    stream.stat() << "h1x" << " ";
    stream.stat() << "h1y" << " ";
    stream.stat() << "h1z" << " ";
    stream.stat() << "h2x" << " ";
    stream.stat() << "h2y" << " ";
    stream.stat() << "h2z" << " ";
    stream.stat() << "h3x" << " ";
    stream.stat() << "h3y" << " ";
    stream.stat() << "h3z" << std::endl;
  }
  
  /********** MC LOOP **********/
  tic = time(NULL);
  for(kStep = 0; kStep < nStep; kStep++) {
    for(kkStep = 0; kkStep < natoms; kkStep++) {
      /********** TRIAL MOVE ATOM **********/
      if(hStep*natoms*Util::RNG::uniform_SFMT() >= 1.0) {
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
      }//end if(hStep*natoms*...
      /********** TRIAL MOVE BOX **********/
      else {
	mc.save();
	try {
	  mcAgent.pertAxes(dhmax);
	  dE = mc.computeEnergy()-(oE+dEt);
	  new_v = mol.volume();
	  dW = dE+pressure*(new_v-old_v)-natoms*kT*std::log(new_v/old_v);
	  accept = (dW <= 0 || 
		    std::log(1.0-Util::RNG::uniform_SFMT()) <= -dW/kT);
	}
	catch(...) {
	  stream.info(false) << "ERROR: Bad h matrix" << std::endl;
	  mc.final();
	  return;
	}
	if(accept) {//accept
	  old_v = new_v;
	  dEt += dE;
	  h_acc_count++;
	}
	else {//reject: undo move
	  mc.restore();
	}
	mcAgent.accPertAxes(accept);
	h_tot_count++;
	//record
	if(stream.dAction()) {
	  stream.action() << "h" << " ";
	  if(accept) {
	    Util::vdump(stream.action(),mol.axis(0).x,3) << " ";
	    Util::vdump(stream.action(),mol.axis(1).x,3) << " ";
	    Util::vdump(stream.action(),mol.axis(2).x,3) << std::endl;   
	    stream.action() << std::endl;
	  }
	  else {
	    stream.action() << "r" << std::endl ;
	  }
	}
      }//end else(hStep*natoms*...
      /********** ADJUST dxmax & dhmax **********/
      if(tot_count == 100) {
	acc_ratio = double(acc_count)+acc_ratio*prev_tot_count;
	acc_ratio /= double(tot_count)+prev_tot_count;
	prev_tot_count += tot_count;
	acc_count = 0;
	tot_count = 0;
	if(reset_flag)
	  dxmax *= (1.0+acc_ratio-0.5);
      }
      if(h_tot_count == 100) {
	h_acc_ratio = double(h_acc_count)+h_acc_ratio*prev_h_tot_count;
	h_acc_ratio /= double(h_tot_count)+prev_h_tot_count;
	prev_h_tot_count += h_tot_count;
	h_acc_count = 0;
	h_tot_count = 0;
	if(reset_flag)
	  dhmax *= (1.0+h_acc_ratio-0.5);
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
    //dump info
    if((kStep+1)%stream.dInfoStep() == 0) {
      toc = time(NULL);
      stream.info() << std::noshowpos << kStep+1 << " "; 
      stream.info() << std::showpos << oE << " "; 
      stream.info() << std::noshowpos << acc_ratio << " "; 
      stream.info() << std::noshowpos << h_acc_ratio << " ";
      stream.info() << std::noshowpos << dxmax << " ";
      stream.info() << std::noshowpos << dhmax << " ";
      stream.info() << std::noshowpos << toc-tic << std::endl;
      tic = time(NULL);
    }//end if((kStep+1)%stream.dInfoStep() == 0)
    //dump stat
    if((kStep+1)%stream.dStatStep() == 0) {
      stream.stat() << std::showpos << oE << " ";
      Util::vdump(stream.stat(),mol.axis(0).x,3) << " ";
      Util::vdump(stream.stat(),mol.axis(1).x,3) << " ";
      Util::vdump(stream.stat(),mol.axis(2).x,3) << std::endl;
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
  std::cout << "[Total time] " << ttoc-ttic << " seconds" << std::endl;
}

}/* MC */
