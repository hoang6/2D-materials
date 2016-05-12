#include "SWCreator.h"
#include "SWRelaxor.h"
#include "../AMod/Molecule.h"
#include "../AMod/MolAnal.h"
#include "../Util/util.h"
#include "../Util/constants.h"

namespace SW {

SWCreator::SWCreator() {}

SWCreator::SWCreator(AMod::Molecule& _mol): AMod::MolServ(_mol) {}

SWStimulus SWCreator::act(const SWStimulus& _stim) {
  SWStimulus stimF = _stim;
  SWStimulus stimB;
  AMod::Molecule& mol = molecule();
  AMod::Bond& bond = mol.bond(stimF.bondID);

  int j, k;
  int left[2], right[2];
  if(stimF.rotCCW) {
    for(k = 0; k < 2; k++)
      AMod::MolAnal::arrowNeighbs(bond.arrow(k),stimF.rotAxis,left[k],right[k]);
  }
  else {
    for(k = 0; k < 2; k++)
      AMod::MolAnal::arrowNeighbs(bond.arrow(k),stimF.rotAxis,right[k],left[k]);
  }
  /*************** reset xyz ***************/
  int id[2];
  double newx[2][3];
  for(k = 0; k < 2; k++) {
    AMod::Arrow& ar = bond.arrow(k);
    double v[3];
    for(j = 0; j < 3; j++) v[j] = -0.5*ar.dx[j];
    if(stimF.rotCCW)
      Util::rot(stimF.rotAxis, +Util::PI/2.0, v);
    else
      Util::rot(stimF.rotAxis, -Util::PI/2.0, v);
    id[k] = ar.host().id();
    for(j = 0; j < 3; j++) 
      newx[k][j] = ar.host()->x[j] + 0.5*ar.dx[j]+v[j];
  }
  mol.moveAtom(id[0], newx[0]);
  mol.moveAtom(id[1], newx[1]);
  //*************** prepare atoms & bonds ***************/
  //self 0, left 1, right 2
  AMod::Arrow* par[2][3] = {{&bond.arrow(0), NULL, NULL},
			    {&bond.arrow(1), NULL, NULL}};
  AMod::Atom* pa[2][3] = {{&bond.arrow(0).host(), NULL, NULL},
			  {&bond.arrow(1).host(), NULL, NULL}};
  for(k = 0; k < 2; k++) {
    if(left[k]) {
      par[k][1] = &par[k][0]->find(left[k]);
      pa[k][1] = &par[k][1]->moon();
    }
    if(right[k]) {
      par[k][2] = &par[k][0]->find(right[k]);
      pa[k][2] = &par[k][2]->moon();
    }
  }
  //break (left) bond
  for(k = 0; k < 2; k++) {
    if(par[k][1])
      mol.topo().deleteBond(par[k][1]->bond().id());
  }
  //new bond
  for(k = 0; k < 2; k++) {
    if(pa[1-k][1])
      mol.topo().newBond(pa[k][0]->id(),pa[1-k][1]->id());
  }
  //relax
  SWRelaxor(mol).relax(bond.id());

  stimB = SWStimulus(bond.id(), !stimF.rotCCW, stimF.rotAxis);
  return stimB;
}

SWCreator::KeyAtoms SWCreator::quad(const AMod::Molecule& _mol, 
				    const SWStimulus& _stim) {
  int k, left[2], right[2];
  KeyAtoms katoms;
  const AMod::Bond& _bond = _mol.bond(_stim.bondID);
  
  katoms.id[0][0] = _bond.arrow(0).host().id();
  katoms.id[1][0] = _bond.arrow(1).host().id();
  katoms.id[0][1] = katoms.id[1][0];
  katoms.id[1][1] = katoms.id[0][0];
  katoms.id[0][2] = katoms.id[1][0];
  katoms.id[1][2] = katoms.id[0][0];
  if(_stim.rotCCW) {
    for(k = 0; k < 2; k++)
      AMod::MolAnal::arrowNeighbs(_bond.arrow(k),_stim.rotAxis,
				  left[k],right[k]);
  }
  else {
    for(k = 0; k < 2; k++)
      AMod::MolAnal::arrowNeighbs(_bond.arrow(k),_stim.rotAxis,
				  right[k],left[k]);
  }
  for(k = 0; k < 2; k++) {
    if(left[k])
      katoms.id[k][1] = _bond.arrow(k).find(left[k]).moon().id();
    if(right[k])
      katoms.id[k][2] = _bond.arrow(k).find(right[k]).moon().id();
  }

  return katoms;
}

} /* SW */
