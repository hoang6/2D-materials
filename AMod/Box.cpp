#include "Box.h"
#include "constants.h"
#include "../Util/constants.h"
#include <cmath>
#include <algorithm>

namespace AMod {

const int Box::TOP_LEFT = 0;
const int Box::TOP_RIGHT = 1;
const int Box::BOTTOM_LEFT = 2;
const int Box::BOTTOM_RIGHT = 3;
const double Box::SCALE = CC_BOND_LENGTH*2.0*Util::SIN60;

Box::Box() {}

Box::Box(const Box& box) { *this = box; }

Box& Box::operator= (const Box& box) {
  if(this != &box) {
    info = box.info;
    ranges[0] = box.ranges[0];
    ranges[1] = box.ranges[1];
  }
  return *this;
}

void Box::prepare(const int _chiV[2], const int _heightV[2]) {
  for(int k = 0; k < 2; k++) {
    info.chiV[k] = _chiV[k];
    info.heightV[k] = _heightV[k];
  }
  
  //compute iQuad and dQuad and chi
  info.iQuad[0][0] = info.chiV[0];
  info.iQuad[0][1] = info.chiV[1];
  cr2xy(info.iQuad[0][0],info.iQuad[0][1],
	info.dQuad[0][0],info.dQuad[0][1]);
  info.chi = atan(info.dQuad[0][1]/info.dQuad[0][0]);
  
  info.iQuad[1][0] = info.heightV[0];
  info.iQuad[1][1] = info.heightV[1];
  cr2xy(info.iQuad[1][0],info.iQuad[1][1],
	info.dQuad[1][0],info.dQuad[1][1]);

  //compute colRange
  int colValues[4] = {
    0,info.iQuad[0][0],info.iQuad[1][0],
    info.iQuad[0][0]+info.iQuad[1][0]
  };
  std::sort(colValues,colValues+4);
  info.colRange[0] = colValues[0];
  info.colRange[1] = colValues[3];

  //compute rowRange
  int rowValues[4] = {
    0,info.iQuad[0][1],info.iQuad[1][1],
    info.iQuad[0][1]+info.iQuad[1][1]
  };
  std::sort(rowValues,rowValues+4);
  info.rowRange[0] = rowValues[0];
  info.rowRange[1] = rowValues[3];

  computeRanges();
}

bool Box::inQuad(int c, int r) const {
  if(r < info.rowRange[0] || r > info.rowRange[1]) return false;
  return 
    c >= ranges[0][r-info.rowRange[0]] && 
    c <= ranges[1][r-info.rowRange[0]];
}

void Box::move2quad(int c, int r, int& cp, int& rp) const {
  /*
  double beta[2];
  double denom = 
    info.iQuad[0][0]*info.iQuad[1][1]-
    info.iQuad[1][0]*info.iQuad[0][1];
  beta[0] = round(-(c*info.iQuad[1][1]+r*info.iQuad[1][0])/denom);
  beta[1] = round(+(c*info.iQuad[0][1]-r*info.iQuad[0][0])/denom);
  cp = int(c+beta[0]*info.iQuad[0][0]+beta[1]*info.iQuad[1][0]);
  rp = int(r+beta[0]*info.iQuad[0][1]+beta[1]*info.iQuad[1][1]);
  */
  int dc, dr;
  cp = c;
  rp = r;
  if(inQuad(cp,rp)) return;
  for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++) {
      if(i == 0 && j== 0) continue;
      dc = i*info.iQuad[0][0]+j*info.iQuad[1][0];
      dr = i*info.iQuad[0][1]+j*info.iQuad[1][1];
      if(inQuad(cp+dc,rp+dr)) {
	cp += dc;
	rp += dr;
	return;
      }
    }
}

void Box::neighbor(int c, int r, int pos, int& neighb_c, int& neighb_r) {
  if(pos == TOP_LEFT) {
    neighb_c = c-1;
    neighb_r = r+1;
  }
  else if(pos == TOP_RIGHT) {
    neighb_c = c;
    neighb_r = r+1;
  }
  else if(pos == BOTTOM_LEFT) {
    neighb_c = c;
    neighb_r = r-1;
  }
  else if(pos == BOTTOM_RIGHT) {
    neighb_c = c+1;
    neighb_r = r-1;
  }
}

void Box::cr2xy(int c, int r, double& x, double& y) {
  x = c+0.5*r;
  y = Util::SIN60*r;
}

void Box::xy2cr(double x, double y, int& c, int& r) {
  double dr = y/Util::SIN60;
  double dc = x-0.5*dr;
  int ir = int(dr);
  int ic = int(dc);
  double x2, y2;
  double dist2[9];
  for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++) {
      cr2xy(ic+i,ir+j,x2,y2);
      dist2[3*(i+1)+(j+1)] = distance2(x,y,x2,y2);
    }
  int indmin = std::min_element(dist2,dist2+9)-dist2;
  c = ic+(int(indmin/3)-1);
  r = ir+(indmin%3-1);
}

double Box::length(int c, int r) {
  double x,y;
  cr2xy(c,r,x,y);
  return sqrt(x*x+y*y);
}

bool Box::onBoundary(int c, int r) const {
  const int(*p)[2] = info.iQuad;
  if(c*p[0][1] == p[0][0]*(r-p[1][1])+p[0][1]*p[1][0] ||
     c*p[1][1] == p[1][0]*(r-p[0][1])+p[1][1]*p[0][0])
    return true;
  return false;
}

double Box::computeRange(int r, int id) const {
  const int(*p)[2] = info.iQuad;
  if(id == 0)
    return 
      double(p[0][0]*r)/p[0][1];
  else if(id == 1)
    return 
      double(p[0][0]*(r-p[1][1]))/p[0][1]+p[1][0];
  else if(id == 2)
    return
      double(p[1][0]*r)/p[1][1];
  else if(id == 3)
    return
      double(p[1][0]*(r-p[0][1]))/p[1][1]+p[0][0];
  return 0.0;
}

void Box::computeRanges() {
  ranges[0].clear();
  ranges[1].clear();

  int j, k;
  int ids[4], nids = 0;
  if(info.iQuad[0][1] == 0 && info.iQuad[1][1] == 0) return;
  else if(info.iQuad[0][1] == 0) { ids[0] = 2; ids[1] = 3; nids = 2; }
  else if(info.iQuad[1][1] == 0) { ids[0] = 0; ids[1] = 1; nids = 2; }
  else { for(j = 0; j < 4; j++) ids[j] = j; nids = 4; }
  
  double range[4];
  int cmin, cmax, nranges = info.rowRange[1]-info.rowRange[0]+1;
  ranges[0].resize(nranges);
  ranges[1].resize(nranges);
  info.nGrids = 0;
  for(k = info.rowRange[0]; k <= info.rowRange[1]; k++) {
    for(j = 0; j < nids; j++) range[j] = computeRange(k,ids[j]);
    cmin = cmax = 0;
    if(nids == 2) {
      std::sort(range,range+2);
      cmin = ceil(range[0]);
      cmax = floor(range[1]);
    } 
    else if(nids == 4) {
      std::sort(range,range+4);
      cmin = ceil(range[1]);
      cmax = floor(range[2]);
    }
    while(onBoundary(cmin,k) && cmin <= cmax) cmin++;
    while(onBoundary(cmax,k) && cmin <= cmax) cmax--;
    ranges[0][k-info.rowRange[0]] = cmin;
    ranges[1][k-info.rowRange[0]] = cmax;
    if(cmax >= cmin) info.nGrids += (cmax-cmin+1);
  }//end for(k...
}

double Box::distance2(double x1, double y1, double x2, double y2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

} /* AMod */
