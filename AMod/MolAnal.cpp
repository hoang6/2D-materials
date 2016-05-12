#include "MolAnal.h"
#include "../Util/util.h"
#include "../Util/constants.h"
#include <iostream>
#include <algorithm>
#include <set>

namespace AMod {

const int MolAnal::DEFAULT_MAX_CHAIN_LENGTH = 20;

MolAnal::BondMark::BondMark() { pc[0] = pc[1]  = NULL; }

MolAnal::MolAnal(int _maxChainLen) { 
  maxChainLength(_maxChainLen); 
}

void MolAnal::reset() { maxChainLength(DEFAULT_MAX_CHAIN_LENGTH); }

bool MolAnal::arrowNorm(const Arrow& ar, double* n, const double* ref) {
  double _n[3];
  const Atom& atom = ar.host();
  int narrows = atom.arrows.size();
  const Arrow* par_i, *par_j;
  if(narrows == 0) {
    _n[0] = _n[1] = _n[2] = 0.0;
  }
  else if(narrows == 1) {
    if(ref != NULL) {
      double z[3];
      Util::cross3(ref, ar.dx, z);
      Util::cross3(ar.dx, z, _n);
    }
    else _n[0] = _n[1] = _n[2] = 0.0;
  }
  else if(narrows == 2) {
    const Arrow& _ar = ar.find(+1);
    Util::cross3(ar.dx, _ar.dx, _n);
  }
  else if(narrows == 3) {
    double oi[3], oj[3];
    par_i = &ar.find(-1);
    par_j = &ar.find(+1);
    Util::vsub(par_i->dx, ar.dx, oi, 3);
    Util::vsub(par_j->dx, ar.dx, oj, 3);
    Util::cross3(oi,oj,_n);
  }
  else {
    double z[3];
    _n[0] = _n[1] = _n[2] = 0.0;
    for(int k = 0; k < narrows; k++) {
      if(k == ar.id()) continue;
      Util::cross3(ar.dx, atom.arrows[k].dx, z);
      if(Util::norm2(z,3) < Util::EPS_DOUBLE) continue;
      Util::normalize(z,3);
      if(Util::vdotv(z,_n,3) < 0) Util::negv(z,3);
      Util::vadd(z,_n,_n,3);
    }
  }
  //post processing
  bool flag = true;
  if(Util::norm2(_n,3) < Util::EPS_DOUBLE) { 
    if(ref != NULL) Util::vcopy(ref, _n, 3);
    else _n[0] = _n[1] = _n[2] = 0.0;
  }
  if(Util::norm2(_n,3) < Util::EPS_DOUBLE) flag = false;
  else Util::normalize(_n,3);
  if(ref != NULL && Util::vdotv(ref,_n,3) < 0) Util::negv(_n,3);
  Util::vcopy(_n,n,3);
  return flag;
}

bool MolAnal::bondNorm(const Bond& bond, double* n, const double* ref) {
  int k;
  double nij[2][3], _n[3];
  for(k = 0; k < 2; k++) arrowNorm(bond.arrow(k), nij[k], ref);
  if(Util::vdotv(nij[0],nij[1],3) < 0) Util::negv(nij[1],3);
  Util::kdotvadd(0.5,nij[0],0.5,nij[1],_n,3);
  //post processing
  bool flag = true;
  if(Util::norm2(_n,3) < Util::EPS_DOUBLE) { 
    if(ref != NULL) Util::vcopy(ref, _n, 3);
    else _n[0] = _n[1] = _n[2] = 0.0;
  }
  if(Util::norm2(_n,3) < Util::EPS_DOUBLE) flag = false;
  else Util::normalize(_n,3);
  if(ref != NULL && Util::vdotv(ref,_n,3) < 0) Util::negv(_n,3);
  Util::vcopy(_n,n,3);
  return flag;
}
 
bool MolAnal::arrowNeighbs(const Arrow& ar, const double* axis, 
			   int& left, int& right) {
  double n[3], h[3], dx_o[3], dx_i[3], dx_j[3];
  double th_i, th_j;
  const Arrow* par_i, *par_j;
  int narrows = ar.host().arrows.size();
  //------ host atom has 1 arrow ------
  if(narrows < 2) {
    left = right = 0;
    return true;
  }
  else if(axis == NULL || Util::norm2(axis,3) < Util::EPS_DOUBLE) {
    left = right = 0;
    return false;
  }
  //------ host atom has 2 arrows ------
  else if(narrows == 2) {
    const Arrow& _ar = ar.find(+1);
    Util::cross3(ar.dx, _ar.dx, n);
    left = (Util::vdotv(n,axis,3)<0?0:1);
    right = 1-left;
  }
  //------ host atom has 3 arrows ------
  else if(narrows == 3) {
    arrowNorm(ar, n, axis);
    left = -1;
    right = +1;
    par_i = &ar.find(left);
    par_j = &ar.find(right);
    //flatten
    Util::project(ar.dx, n, h, 3);
    Util::vsub(ar.dx, h, dx_o, 3);
    Util::vsub(par_i->dx, h, dx_i, 3);
    Util::vsub(par_j->dx, h, dx_j, 3);
    //dx_o->x, n->z
    th_i = angle(dx_o, n, dx_i);
    th_j = angle(dx_o, n, dx_j);
    //determine left and right
    if(th_i > th_j) std::swap(left,right);
  }
  //------ host atom has more than 3 arrows ------
  else {
    std::vector< std::pair<double,int> > order;
    arrowOrder(ar,axis,order);
    left = order.front().second;
    right = order.back().second;
  }
  return true;
}

bool MolAnal::ring(const Arrow& ar, Chain& chain, const double* ref) const {
  int left, right;
  double _ref[3];
  const Arrow* par = &ar;
  if(!bondNorm(par->bond(), _ref, ref)) return 0;

  //debug
  //Count++;
  //cout << "-------------------- " << Count << " --------------------" << endl;
  chain.clear();
  while(!(chain.size() > 0 && par == &ar)) {
    if(int(chain.size()) == maxChainLength()) {
      //cerr << "Warning<MolAnal::getRing>: hit chain size limit " 
      //   << maxChainLen << endl;
      return 0;
    }
    chain.push_back(par->idpos());
    if(!bondNorm(par->bond(), _ref, _ref)) return 0;
    /*debug
    cout << "********************" << endl;
    cout << par->host().x[0] << " " 
	 << par->host().x[1] << " " 
	 << par->host().x[2] << endl;
    cout << par->dx[0] << " " 
	 << par->dx[1] << " " 
	 << par->dx[2] << endl;
    cout << n[0] << " " << n[1] << " " << n[2] << endl;
    */
    par = &par->brother();
    arrowNeighbs(*par, _ref, left, right);
    //return 0 if arrow has no right neighbs
    if(right == 0) return 0;
    par = &par->find(right);    
  }//end while(1)
  return 1;
}

void MolAnal::rings(const Molecule& mol, Chains& chains) {
  int k, nbonds = mol.nbonds();
  std::vector<int> bond_ids(nbonds);
  for(k = 0; k < nbonds; k++) bond_ids[k] = k;
  rings(mol, bond_ids.begin(), bond_ids.end(), chains);
}

void MolAnal::rings(const Molecule& mol, int bondID, 
		    Chains& chainsA, Chains& chainsB) {
  const Bond& bond = mol.bond(bondID);
  const Atom* patoms[2] = {&bond.arrow(0).host(), &bond.arrow(1).host()};
  Chains chains_tmp;
  std::set<int> bond_ids;
  int j, k, chains_sz, chain_sz;
  bool flag;

  chainsA.clear();
  chainsB.clear();

  for(k = 0; k < 2; k++)
    for(j = 0; j < int(patoms[k]->arrows.size()); j++)
      bond_ids.insert(patoms[k]->arrows[j].bond().id());
  rings(mol, bond_ids.begin(), bond_ids.end(), chains_tmp);
  chains_sz = chains_tmp.size();
  for(k = 0; k < chains_sz; k++) {
    chain_sz = chains_tmp[k].size();
    for(j = 0, flag = 0; j < chain_sz; j++) 
      if(chains_tmp[k][j].first == bond.id()) { flag = 1; break; }
    if(flag)
      chainsA.push_back(chains_tmp[k]);
    else
      chainsB.push_back(chains_tmp[k]);
  }
}

bool MolAnal::stair(const Arrow& ar, Chain& chain, const double* ref) const {
  int key, left, right;
  double _ref[3];
  const Arrow* par = &ar;
  if(!bondNorm(par->bond(), _ref, ref)) return 0;

  chain.clear();
  while(!(chain.size() > 0 && par == &ar)) {
    if(int(chain.size()) == maxChainLength()) return 1;
    chain.push_back(par->idpos());
    if(!bondNorm(par->bond(), _ref, _ref)) return 0;
    par = &par->brother();
    arrowNeighbs(*par, _ref, left, right);
    key = (chain.size()%2 ? left : right);
    //return 0 if arrow has no left/right neighbs
    if(key == 0) return 0;
    par = &par->find(key);    
  }//end while(...
  return 1;
}

void MolAnal::arrowOrder(const Arrow& ar, const double* axis,
			 std::vector< std::pair<double,int> >& order) {
  int k, n = ar.host().arrows.size();
  double th, h[3], dx_o[3], dx_k[3];
  Util::project(ar.dx, axis, h, 3);
  Util::vsub(ar.dx, h, dx_o, 3);
  order.resize(n-1);
  for(k = 1; k < n; k++) {
    const Arrow& _ar = ar.find(k);
    Util::vsub(_ar.dx, h, dx_k, 3);
    th = angle(dx_o, axis, dx_k);
    order[k-1] = std::pair<double,int>(th,k);
  }//end for(k...
  sort(order.begin(),order.end());
}

double MolAnal::angle(const double* x, const double* z, const double* t) {
  static const double pi = atan(1.0)*4;
  double y[3], ty, tx, th;
  Util::cross3(z, x, y);
  ty = Util::vdotv(t, y, 3)/sqrt(Util::vdotv(y, y, 3));
  tx = Util::vdotv(t, x, 3)/sqrt(Util::vdotv(x, x, 3));
  th = atan2(ty, tx);
  if(th < 0) th = th+2*pi;
  return th;
}

void MolAnal::regChain(Chain& chain) {
  int k, n = chain.size();
  for(k = 0; k < n; k++) {
    //debug
    //cout << chain[k].first << " " << chain[k].second << endl;
    BondMark& mark = marks[chain[k].first];
    if(mark.pc[0]) mark.pc[1] = &chain;
    else mark.pc[0] = &chain;
  }
}

bool MolAnal::checkDup(const Chain& chain) const {
  int k, n = chain.size();
  Chain* pco = NULL;
  Chain* const *pc;
  for(k = 0; k < n; k++) {
    //debug
    //cout << chain[k].first << " " << chain[k].second << endl;
    pc = marks[chain[k].first].pc;
    if(pc[0] && pc[1]) return 1;
    else if(!pc[0] && !pc[1]) return 0;
    if(k == 0) pco = (pc[0]?pc[0]:pc[1]);
    else if(pco != (pc[0]?pc[0]:pc[1])) return 0;
  }//end for(int j...
  return 1;
}

std::ostream& operator<< (std::ostream& fout, const MolAnal::Chain& chain) {
  fout << chain.size() << "\t";
  for(int k = 0; k < int(chain.size()); k++)
    fout << chain[k].first << " " << chain[k].second << "\t";
  return fout;
}

std::ostream& operator<< (std::ostream& fout, const MolAnal::Chains& chains) {
  for(int k = 0, kmax = int(chains.size()); k < kmax; k++) 
    fout << k << " " << chains[k] << std::endl;
  return fout;
}

} /* AMod */
