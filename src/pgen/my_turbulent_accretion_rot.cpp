//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast.cpp
//! \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//!        cylindrical, and spherical coordinates.  Contains post-processing code
//!        to check whether blast is spherical for regression tests
//!
//! REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//!   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

int acc_mode=1; // flag for accretor setup. 0=no accretor.
Real acc_rate=1.0; // accretion efficiency for smooth accretor
Real T0=1.; // T at r>1
Real B0=1.; // ((gamma-1)/gamma * P + Ek) / Phi
bool set_B0=false; // whether set T0 or B0
bool only_damp_v=false; // damp velocity only
bool output_Mdot=false;
Real fcirc=0.0;
Real sink_radius_n_cell=4; // number of cells for sink radius
int n_cell_acc1 = 0; // cells at r<r_sink
int n_cell_acc2 = 0; // cells at r_sink - 1.5 r_sink

int RefinementCondition(MeshBlock *pmb);

void MySource(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);

Real HstOutput(MeshBlock *pmb, int iout);

Real n_cells_between_radii(Real n_in, Real n_out){
  // count how many cells are within a given radius range
  int n_max = n_out+1;
  int n_cell = 0;
  for (Real i=0; i<=n_max; i++) {
    for (Real j=0; j<=n_max; j++) {
      for (Real k=0; k<=n_max; k++) {
        Real r = std::sqrt(SQR(0.5+i)+SQR(0.5+j)+SQR(0.5+k));
        if (r>n_in && r<=n_out) n_cell += 1;
      }
    }
  }
  n_cell *= 8; // all octants
  return n_cell;
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  acc_mode = pin->GetOrAddInteger("problem","acc_mode",acc_mode);
  acc_rate = pin->GetOrAddReal("problem","acc_rate",acc_rate);
  T0 = pin->GetOrAddReal("problem","T0",T0);
  B0 = pin->GetOrAddReal("problem","B0",B0);
  set_B0 = pin->GetOrAddBoolean("problem","set_B0",set_B0);
  only_damp_v = pin->GetOrAddBoolean("problem","only_damp_v",only_damp_v);
  output_Mdot = pin->GetOrAddBoolean("problem","output_Mdot",output_Mdot);
  turb_flag = pin->GetInteger("problem","turb_flag");
  fcirc = pin->GetOrAddReal("problem","fcirc",fcirc);
  sink_radius_n_cell = pin->GetOrAddReal("problem","sink_radius_n_cell",sink_radius_n_cell);
  n_cell_acc1 = n_cells_between_radii(0, sink_radius_n_cell);
  n_cell_acc2 = n_cells_between_radii(sink_radius_n_cell, 1.5*sink_radius_n_cell);
  EnrollUserExplicitSourceFunction(MySource);
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
  }
  if (output_Mdot) {
    AllocateUserHistoryOutput(3);
    EnrollUserHistoryOutput(0, HstOutput, "Mdot"); //measured accretion rate
    EnrollUserHistoryOutput(1, HstOutput, "rho_center"); // at accretor for smooth, around accretor for vac
    EnrollUserHistoryOutput(2, HstOutput, "Jdot"); //measured angular momentum accretion rate
  }
  return;
}

namespace {
// Apply a density floor - useful for large |z| regions
Real dfloor, pfloor;
} // namespace


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  dfloor=pin->GetOrAddReal("hydro","dfloor",0.01);
  pfloor=pin->GetOrAddReal("hydro","pfloor",0.01);
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;
    //Initial Conditions
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        // density = pressure = 1
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real r = std::sqrt(SQR(x)+SQR(y)+SQR(z));
        Real vxgoal = -y*std::sqrt(fcirc);
        Real vygoal = x*std::sqrt(fcirc);
        Real rho0 = 1.0;//std::max(pow(1 + (gamma-1)/r - (gamma - 1)/2.0*fcirc/(pow(x, 2) + pow(y, 2)), 1/(gamma-1)),dfloor);
        phydro->u(IDN,k,j,i) = rho0;
        phydro->u(IM1,k,j,i) = rho0*vxgoal;
        phydro->u(IM2,k,j,i) = rho0*vygoal;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS)
          phydro->u(IEN,k,j,i) = rho0*T0/gm1 + .5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))+SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }
  }
}

void MeshBlock::UserWorkInLoop() {
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real& u_d  = phydro->u(IDN,k,j,i);
        u_d = (u_d > dfloor) ?  u_d : dfloor;
        if (NON_BAROTROPIC_EOS) {
          Real gam = 1.333333;
          Real& w_p  = phydro->w(IPR,k,j,i);
          Real& u_e  = phydro->u(IEN,k,j,i);
          Real& u_m1 = phydro->u(IM1,k,j,i);
          Real& u_m2 = phydro->u(IM2,k,j,i);
          Real& u_m3 = phydro->u(IM3,k,j,i);
          w_p = (w_p > pfloor) ?  w_p : pfloor;
          //Real di = 1.0/u_d;
          //Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
          //u_e = w_p/(gam-1.0)+ke;
        }
      }
    }
  }
  return;
}



void MySource(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar) {
  Real GM = 1.0;
  // Real T0 = 1.0; // T0 is now set by input
  Real t_cool = 1.0;
  Real gamma = pmb->peos->GetGamma();
  Real gm1 = gamma - 1.0;
  Real r_smooth; // scale below which g is smoothed; also the sink radius
  r_smooth = sink_radius_n_cell*pmb->pcoord->dx1f(0);
  
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // step 1. central point mass
        Real x = pmb->pcoord->x1v(i);
        Real y = pmb->pcoord->x2v(j);
        Real z = pmb->pcoord->x3v(k);
        Real r = std::sqrt(SQR(x)+SQR(y)+SQR(z));
        // r for on interfaces
        Real r1l = std::sqrt(SQR(pmb->pcoord->x1f(i  ))+SQR(pmb->pcoord->x2v(j))+SQR(pmb->pcoord->x3v(k)));
        Real r1r = std::sqrt(SQR(pmb->pcoord->x1f(i+1))+SQR(pmb->pcoord->x2v(j))+SQR(pmb->pcoord->x3v(k)));
        Real r2l = std::sqrt(SQR(pmb->pcoord->x1v(i))+SQR(pmb->pcoord->x2f(j  ))+SQR(pmb->pcoord->x3v(k)));
        Real r2r = std::sqrt(SQR(pmb->pcoord->x1v(i))+SQR(pmb->pcoord->x2f(j+1))+SQR(pmb->pcoord->x3v(k)));
        Real r3l = std::sqrt(SQR(pmb->pcoord->x1v(i))+SQR(pmb->pcoord->x2v(j))+SQR(pmb->pcoord->x3f(k  )));
        Real r3r = std::sqrt(SQR(pmb->pcoord->x1v(i))+SQR(pmb->pcoord->x2v(j))+SQR(pmb->pcoord->x3f(k+1)));
        // phi = -1/r at large radii, phi = const at small radii
        const Real r_smooth_3 = r_smooth*r_smooth*r_smooth;
        Real phi = (r>r_smooth) ? -1./r : -1./r_smooth + .5*(SQR(r)-SQR(r_smooth))/r_smooth_3;
        Real phi1l = (r1l>r_smooth) ? -1./r1l : -1./r_smooth + .5*(SQR(r1l)-SQR(r_smooth))/r_smooth_3;
        Real phi1r = (r1r>r_smooth) ? -1./r1r : -1./r_smooth + .5*(SQR(r1r)-SQR(r_smooth))/r_smooth_3;
        Real phi2l = (r2l>r_smooth) ? -1./r2l : -1./r_smooth + .5*(SQR(r2l)-SQR(r_smooth))/r_smooth_3;
        Real phi2r = (r2r>r_smooth) ? -1./r2r : -1./r_smooth + .5*(SQR(r2r)-SQR(r_smooth))/r_smooth_3;
        Real phi3l = (r3l>r_smooth) ? -1./r3l : -1./r_smooth + .5*(SQR(r3l)-SQR(r_smooth))/r_smooth_3;
        Real phi3r = (r3r>r_smooth) ? -1./r3r : -1./r_smooth + .5*(SQR(r3r)-SQR(r_smooth))/r_smooth_3;
        // conservative source term from self gravity module
        AthenaArray<Real> *flux = pmb->phydro->flux;
        Real dx1 = pmb->pcoord->dx1v(i);
        Real dtodx1 = dt/dx1;
        cons(IM1,k,j,i) -= dtodx1*prim(IDN,k,j,i)*(phi1r-phi1l);
        if (NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) -= dtodx1*(flux[X1DIR](IDN,k,j,i  )*(phi - phi1l) +
                                     flux[X1DIR](IDN,k,j,i+1)*(phi1r - phi));
        Real dx2 = pmb->pcoord->dx2v(j);
        Real dtodx2 = dt/dx2;
        cons(IM2,k,j,i) -= dtodx2*prim(IDN,k,j,i)*(phi2r-phi2l);
        if (NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) -= dtodx2*(flux[X2DIR](IDN,k,j  ,i)*(phi - phi2l) +
                                     flux[X2DIR](IDN,k,j+1,i)*(phi2r - phi));
        Real dx3 = pmb->pcoord->dx3v(k);
        Real dtodx3 = dt/dx3;
        cons(IM3,k,j,i) -= dtodx3*prim(IDN,k,j,i)*(phi3r-phi3l);
        if (NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) -= dtodx3*(flux[X3DIR](IDN,k  ,j,i)*(phi - phi3l) +
                                     flux[X3DIR](IDN,k+1,j,i)*(phi3r - phi));
        // now accrete... acc_mode=0 means no accretor.
        // accrete within r_smooth
        if (r<r_smooth) {
          const Real dt_dyn = dt * std::pow(r_smooth, -1.5);
          if (acc_mode==1) { // remove mass, keep velocity and specific momentum
            cons(IDN,k,j,i) -= acc_rate*dt_dyn*cons(IDN,k,j,i);
            cons(IM1,k,j,i) -= acc_rate*dt_dyn*cons(IM1,k,j,i);
            cons(IM2,k,j,i) -= acc_rate*dt_dyn*cons(IM2,k,j,i);
            cons(IM3,k,j,i) -= acc_rate*dt_dyn*cons(IM3,k,j,i);
            if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -= acc_rate*dt_dyn*cons(IEN,k,j,i);
          }
          else if (acc_mode==2) { // enforce floor density + zero velocity
            cons(IDN,k,j,i) = 1.0;
            cons(IM1,k,j,i) = 0.;
            cons(IM2,k,j,i) = 0.;
            cons(IM3,k,j,i) = 0.;
            if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) = 1.0/gm1;
          }
          else if (acc_mode==-1) { // inject mommentum (but no energy): do this on everything in sink region
            // velocity increases by r/t_dyn over one t_dyn
            const Real dt_a = dt*r/r_smooth_3;
            cons(IM1,k,j,i) += prim(IDN,k,j,i)*dt_a*x/r;
            cons(IM2,k,j,i) += prim(IDN,k,j,i)*dt_a*y/r;
            cons(IM3,k,j,i) += prim(IDN,k,j,i)*dt_a*z/r;
          }
          else if (acc_mode==-2) { // inject mommentum (but no energy): do this only to inflow
            const Real dt_a = dt*r/r_smooth_3;
            bool is_inflow = (prim(IVX,k,j,i)*x+prim(IVY,k,j,i)*y+prim(IVZ,k,j,i)*z)<0.;
            if (is_inflow) {
              cons(IM1,k,j,i) += prim(IDN,k,j,i)*dt_a*x/r;
              cons(IM2,k,j,i) += prim(IDN,k,j,i)*dt_a*y/r;
              cons(IM3,k,j,i) += prim(IDN,k,j,i)*dt_a*z/r;
            }
          }
          else if (acc_mode==-3) { // inject mommentum (but no energy): do this only to outflow
            const Real dt_a = dt*r/r_smooth_3;
            bool is_outflow = (prim(IVX,k,j,i)*x+prim(IVY,k,j,i)*y+prim(IVZ,k,j,i)*z)>0.;
            if (is_outflow) {
              cons(IM1,k,j,i) += prim(IDN,k,j,i)*dt_a*x/r;
              cons(IM2,k,j,i) += prim(IDN,k,j,i)*dt_a*y/r;
              cons(IM3,k,j,i) += prim(IDN,k,j,i)*dt_a*z/r;
            }
          }
        }
        if (r>1.) {
          // step 2. cooling: cool towards T=p/rho=T0
          if (NON_BAROTROPIC_EOS) {
            if (only_damp_v) { // only damp velocity
              Real Ek = .5*prim(IDN,k,j,i)*(
                        SQR(prim(IVX,k,j,i))+
                        SQR(prim(IVY,k,j,i))+
                        SQR(prim(IVZ,k,j,i)));
              cons(IEN,k,j,i) -= Ek * dt/t_cool * 2;
              Real P_goal = (prim(IDN,k,j,i)/r*B0 - Ek)/(gamma/gm1);
              P_goal = std::max(P_goal, 1.e-8);
              cons(IEN,k,j,i) += (P_goal-prim(IPR,k,j,i))/gm1 * (1.-std::exp(-10.*dt/t_cool));
            } else if (set_B0) { // cool towards B0
              Real T_B0 = gm1/gamma*B0/r;
              cons(IEN,k,j,i) -= dt/t_cool * (prim(IPR,k,j,i)-T_B0*prim(IDN,k,j,i))/gm1;
              Real Ek = .5*prim(IDN,k,j,i)*(
                        SQR(prim(IVX,k,j,i))+
                        SQR(prim(IVY,k,j,i))+
                        SQR(prim(IVZ,k,j,i)));
              cons(IEN,k,j,i) -= Ek * dt/t_cool * 2;
            } else { // cool towards T=p/rho=T0
              cons(IEN,k,j,i) -= dt/t_cool * (prim(IPR,k,j,i)-T0*prim(IDN,k,j,i))/gm1;
            }
          }
          // step 3. mass removal towards rho = 1.
            Real bracket = std::max(1 + (gamma-1)/r - (gamma - 1)/2.0*fcirc/(pow(x, 2) + pow(y, 2)),dfloor);
            Real rhogoal = std::max(pow(bracket, 1/(gamma-1)),dfloor); //define a goal density
          //Real new_rho_fraction = 1.*(1.-dt/t_cool/rhogoal) + 1./prim(IDN,k,j,i)*(dt/t_cool); //This doesn't work
            Real new_rho_fraction = 1.-dt/t_cool*(prim(IDN,k,j,i)-rhogoal)/prim(IDN,k,j,i); //damping towards a goal density, turns out to be equivalent
          cons(IDN,k,j,i) *= new_rho_fraction;
          cons(IM1,k,j,i) *= new_rho_fraction;
          cons(IM2,k,j,i) *= new_rho_fraction;
          cons(IM3,k,j,i) *= new_rho_fraction;
          if (NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) *= new_rho_fraction;
            
          // step 4. damp velocity
            Real damp_faction = 1.-dt/t_cool;
            //Real vxgoal = -(y / (pow(x, 2) + pow(y, 2)))*std::sqrt(fcirc); //a goal vx
            //Real vygoal = (x / (pow(x, 2) + pow(y, 2)))*std::sqrt(fcirc); //a goal vy
            Real vxgoal = -y*std::sqrt(fcirc); //a goal vx
            Real vygoal = x*std::sqrt(fcirc); //a goal vy
            // above: changed to uniform rotation
            cons(IM1,k,j,i) += prim(IDN,k,j,i)*dt/t_cool*(vxgoal-prim(IVX,k,j,i)) ;
            cons(IM2,k,j,i) += prim(IDN,k,j,i)*dt/t_cool*(vygoal-prim(IVY,k,j,i)) ;
            cons(IM3,k,j,i) *= 1.-dt/t_cool;
            
            //Real signx = std::abs(prim(IVX,k,j,i)*vxgoal)/(prim(IVX,k,j,i)*vxgoal);
            //Real signy = std::abs(prim(IVY,k,j,i)*vygoal)/(prim(IVY,k,j,i)*vygoal);
            
          // step 5. force b=0 //Bernoulli constant, optional
            
          /*if (only_damp_v) {
            Real ek = .5/SQR(cons(IDN,k,j,i))*(
                      SQR(cons(IM1,k,j,i))+
                      SQR(cons(IM2,k,j,i))+
                      SQR(cons(IM3,k,j,i)));
            Real p = std::max(1.e-8,(1./r - ek)/gamma*gm1);
            cons(IEN,k,j,i) = p/gm1 + cons(IDN,k,j,i)*ek;
          }*/
        }
        // step -1. driving turbulence: not used now
      }
    }
  }
}

// refinement condition: check the maximum pressure gradient
int RefinementCondition(MeshBlock *pmb) {
  Real P1 = pmb->pcoord->x1v(pmb->is-1) * pmb->pcoord->x1v(pmb->ie+1);
  Real P2 = pmb->pcoord->x2v(pmb->js-1) * pmb->pcoord->x2v(pmb->je+1);
  Real P3 = pmb->pcoord->x3v(pmb->ks-1) * pmb->pcoord->x3v(pmb->ke+1);
  if (P1<0. && P2<0. && P3<0.)
    return 1;
  else
    return 0;
}

// hst outputs
Real HstOutput(MeshBlock *pmb, int iout) {
  Real A = 0.;
  const Real r_sink = sink_radius_n_cell*pmb->pcoord->dx1f(0);
  const Real dV = pmb->pcoord->dx1f(0)*pmb->pcoord->dx2f(0)*pmb->pcoord->dx3f(0);
  const Real OmegaK = std::pow(r_sink, -1.5);
  // below: these are now computed during initialization
  //const Real n_cell_acc1 = 35.; //280.; // < r_sink = 4 cells
  //const Real n_cell_acc2 = 79.; //632.; // 1-1.5 r_sink
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real x = pmb->pcoord->x1v(i);
        Real y = pmb->pcoord->x2v(j);
        Real z = pmb->pcoord->x3v(k);
        Real r = std::sqrt(SQR(x)+SQR(y)+SQR(z));
        if (iout==0) { // accretion rate, then both Mdot and rho_center are displaying accretion rate
          if (acc_mode==1 && r<r_sink)
            A += dV*pmb->phydro->u(IDN,k,j,i)*OmegaK*acc_rate;
          else if (acc_mode==2 && r>r_sink && r<1.5*r_sink)
            A -= PI*r*
                 (pmb->phydro->u(IM1,k,j,i)*x + 
                  pmb->phydro->u(IM2,k,j,i)*y + 
                  pmb->phydro->u(IM3,k,j,i)*z) / n_cell_acc2;
        }
        if (iout==1) { // mean density
          if (acc_mode==1 && r<r_sink)
            A += pmb->phydro->u(IDN,k,j,i) / n_cell_acc1;
          else if (acc_mode==2 && r>r_sink && r<1.5*r_sink)
            A += pmb->phydro->u(IDN,k,j,i) / n_cell_acc2;
        }
        else if (iout==2){ //total Lz
            if (acc_mode==1 && r<r_sink)
                A+= dV*OmegaK*acc_rate*(x*pmb->phydro->u(IM2,k,j,i) - y*pmb->phydro->u(IM1,k,j,i));
        }
      }
    }
  }
  return A;
}
