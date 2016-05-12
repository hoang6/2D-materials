#include "theta2.h"
#include "globalConstants.h"

namespace tbtools {

void cross(const double* a, const double* b, double* c) {
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = b[0]*a[2]-b[2]*a[0];
  c[2] = a[0]*b[1]-b[0]*a[1];
}

double vdotv(const double* a, const double* b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double theta2(const double* a, const double* b, const double* r) {
  double ejik[3], eijl[3], m1, m2, theta_ij_2;
  cross(a, r, ejik);
  cross(b, r, eijl);
  m1 = vdotv(ejik, ejik);
  m2 = vdotv(eijl, eijl);
  if(m1 < EPS_DOUBLE || m2 < EPS_DOUBLE) theta_ij_2 = 0.0;
  else {
    theta_ij_2 = vdotv(ejik, eijl);
    theta_ij_2 = theta_ij_2*theta_ij_2/m1/m2;
  }
  return theta_ij_2;
}
  
void dtheta2_dr(const double* a, const double* b, const double* r, double* result) {
  double ejik[3], eijl[3], m1, m2, eedot;
  cross(a, r, ejik);
  cross(b, r, eijl);
  m1 = vdotv(ejik, ejik);
  m2 = vdotv(eijl, eijl);
  if(m1 < EPS_DOUBLE || m2 < EPS_DOUBLE) { 
    result[0] = result[1] = result[2] = 0.0;
    return;
  }
  eedot = vdotv(ejik, eijl);
  //------ d(theta2)/drx ------
  result[0] = 
    (eedot*(-2*eedot*m2*(+a[2]*ejik[1]+a[1]*a[1]*r[0]-a[0]*a[1]*r[1])-2*eedot*m1*(b[2]*eijl[1]+b[1]*b[1]*r[0]-b[0]*b[1]*r[1])+2*m1*m2*(2*a[1]*b[1]*r[0]+2*a[2]*b[2]*r[0]-a[1]*b[0]*r[1]-a[0]*b[1]*r[1]-a[2]*b[0]*r[2]-a[0]*b[2]*r[2])))/(m1*m1*m2*m2);
  //------ d(theta2)/dry ------
  result[1] = 
    (eedot*(-2*eedot*m2*(-a[2]*ejik[0]-a[0]*a[1]*r[0]+a[0]*a[0]*r[1])-2*eedot*m1*(-b[2]*eijl[0]-b[0]*b[1]*r[0]+b[0]*b[0]*r[1])+2*m1*m2*(-a[1]*b[0]*r[0]-a[0]*b[1]*r[0]+2*a[0]*b[0]*r[1]+2*a[2]*b[2]*r[1]-a[2]*b[1]*r[2]-a[1]*b[2]*r[2])))/(m1*m1*m2*m2);
  //------ d(theta2)/drz ------
  result[2] = 
    (eedot*(2*m1*m2*(-a[2]*b[0]*r[0]-a[0]*b[2]*r[0]-a[2]*b[1]*r[1]-a[1]*b[2]*r[1]+2*a[0]*b[0]*r[2]+2*a[1]*b[1]*r[2])-2*eedot*m2*(a[0]*(-a[2]*r[0]+a[0]*r[2])+a[1]*(-a[2]*r[1] + a[1]*r[2]))-2*eedot*m1*(b[0]*(-b[2]*r[0]+b[0]*r[2])+b[1]*(-b[2]*r[1]+b[1]*r[2]))))/(m1*m1*m2*m2);
}

void dtheta2_da(const double* a, const double* b, const double* r, double* result) {
  double ejik[3], eijl[3], m1, m2, eedot;
  cross(a, r, ejik);
  cross(b, r, eijl);
  m1 = vdotv(ejik, ejik);
  m2 = vdotv(eijl, eijl);
  if(m1 < EPS_DOUBLE || m2 < EPS_DOUBLE) { 
    result[0] = result[1] = result[2] = 0.0;
    return;
  }
  eedot = vdotv(ejik, eijl);
  //------ d(theta2)/dax ------
  result[0] = 
    (2*eedot*(-eedot*ejik[2]*r[1]+eijl[2]*m1*r[1]+eedot*ejik[1]*r[2]-eijl[1]*m1*r[2]))/(m1*m1*m2);
  //------ d(theta2)/day ------
  result[1] = 
    (2*eedot*(+eedot*ejik[2]*r[0]-eijl[2]*m1*r[0]-eedot*ejik[0]*r[2]+eijl[0]*m1*r[2]))/(m1*m1*m2);
  //------ d(theta2)/daz ------
  result[2] = 
    (2*eedot*(-eedot*ejik[1]*r[0]+eijl[1]*m1*r[0]+eedot*ejik[0]*r[1]-eijl[0]*m1*r[1]))/(m1*m1*m2);
}

void dtheta2_db(const double* a, const double* b, const double* r, double* result) {
  dtheta2_da(b, a, r, result);
}

}/* tbtools */
