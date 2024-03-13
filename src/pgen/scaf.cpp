//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scaf.cpp
//! \brief Problem generator to initialize rotational equilibrium tori SCAF (see the ref 
//! below).  Based on my_turbulenct_accretion_rot.cpp in Athena++, with edits by Wenrui Xu 
//! and Yi-Xian Chen. Modified to AthenaK by Hai-Yang Wang.
//! 2024.1.31
//!
//! References:
//!    Xu 2023 ApJ_954_180
//!    Simple Convective Accretion Flows (SCAFs)
//!    Explaining the ~-1 Density Scaling of Hot Accretion Flows around Compact Objects

#include <stdio.h>
#include <math.h>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include <algorithm>  // max(), max_element(), min(), min_element()
#include <iomanip>
#include <iostream>   // endl
#include <limits>     // numeric_limits::max()
#include <memory>
#include <sstream>    // stringstream
#include <string>     // c_str(), string
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "geodesic-grid/spherical_grid.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "radiation/radiation.hpp"

#include <Kokkos_Random.hpp>

namespace {

struct scaf_pgen{

int acc_mode; // flag for accretor setup. 0=no accretor.
Real acc_rate; // accretion efficiency for smooth accretor
Real T0; // T at r>1
Real B0; // ((gamma-1)/gamma * P + Ek) / Phi
bool set_B0; // whether set T0 or B0
bool only_damp_v; // damp velocity only
bool output_Mdot;
Real fcirc;
Real sink_radius_n_cell; // number of cells for sink radius
int n_cell_acc1; // cells at r<r_sink
int n_cell_acc2; // cells at r_sink - 1.5 r_sink
int turb_flag; // flag for turbulence

Real gamma, gm1; // ideal gas EOS data
Real dfloor, pfloor; // Apply a density floor - useful for large |z| regions

};

scaf_pgen scaf;

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

} // namespace


void RefinementCondition(MeshBlockPack* pmbp);

void AddUserSrcs(Mesh *pm, const Real bdt);

void MySource(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                   const DvceArray5D<Real> &w0, const EOS_Data &eos_data);

void HstOutput(HistoryData *pdata, Mesh *pm);


//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief Sets initial conditions for either Fishbone-Moncrief or Chakrabarti torus in GR
//! Compile with '-D PROBLEM=gr_torus' to enroll as user-specific problem generator
//!  assumes x3 is axisymmetric direction

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
// user_bcs_func = RadialBoundary;
user_srcs = true;
user_srcs_func = AddUserSrcs;
if (pmy_mesh_->adaptive) user_ref_func = RefinementCondition;
user_hist_func = HstOutput;
// user_dt_func = AccTimeStep;
// pgen_final_func = AccFinalWork;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int ie = indcs.ie, je = indcs.je, ke = indcs.ke;
  int nmb = pmbp->nmb_thispack;
  auto &coord = pmbp->pcoord->coord_data;

  // Select either Hydro or MHD
  DvceArray5D<Real> u0_, w0_;
  if (pmbp->phydro != nullptr) {
    u0_ = pmbp->phydro->u0;
    w0_ = pmbp->phydro->w0;
  } else if (pmbp->pmhd != nullptr) {
    u0_ = pmbp->pmhd->u0;
    w0_ = pmbp->pmhd->w0;
  }

 scaf.acc_mode=1; // flag for accretor setup. 0=no accretor.
 scaf.acc_rate=1.0; // accretion efficiency for smooth accretor
 scaf.T0=1.; // T at r>1
 scaf.B0=1.; // ((gamma-1)/gamma * P + Ek) / Phi
 scaf.set_B0=false; // whether set T0 or B0
 scaf.only_damp_v=false; // damp velocity only
 scaf.output_Mdot=false;
 scaf.fcirc=0.0;
 scaf.sink_radius_n_cell=4; // number of cells for sink radius
 scaf.n_cell_acc1 = 0; // cells at r<r_sink
 scaf.n_cell_acc2 = 0; // cells at r_sink - 1.5 r_sink
 scaf.turb_flag=0; // flag for turbulence

// void Mesh::InitUserMeshData(ParameterInput *pin) in Athena++
  scaf.acc_mode = pin->GetOrAddInteger("problem","acc_mode",scaf.acc_mode);
  scaf.acc_rate = pin->GetOrAddReal("problem","acc_rate",scaf.acc_rate);
  scaf.T0 = pin->GetOrAddReal("problem","T0",scaf.T0);
  scaf.B0 = pin->GetOrAddReal("problem","B0",scaf.B0);
  scaf.set_B0 = pin->GetOrAddBoolean("problem","set_B0",scaf.set_B0);
  scaf.only_damp_v = pin->GetOrAddBoolean("problem","only_damp_v",scaf.only_damp_v);
  scaf.output_Mdot = pin->GetOrAddBoolean("problem","output_Mdot",scaf.output_Mdot);
  scaf.turb_flag = pin->GetInteger("problem","turb_flag");
  scaf.fcirc = pin->GetOrAddReal("problem","fcirc",scaf.fcirc);
  scaf.sink_radius_n_cell = pin->GetOrAddReal("problem","sink_radius_n_cell",scaf.sink_radius_n_cell);
  scaf.n_cell_acc1 = n_cells_between_radii(0, scaf.sink_radius_n_cell);
  scaf.n_cell_acc2 = n_cells_between_radii(scaf.sink_radius_n_cell, 1.5*scaf.sink_radius_n_cell);

//! void MeshBlock::ProblemGenerator in Athena++

  scaf.dfloor=pin->GetOrAddReal("hydro","dfloor",0.01);
  scaf.pfloor=pin->GetOrAddReal("hydro","pfloor",0.01);
  
  // Get ideal gas EOS data
  if (pmbp->phydro != nullptr) {
    scaf.gamma = pmbp->phydro->peos->eos_data.gamma;
  } else if (pmbp->pmhd != nullptr) {
    scaf.gamma = pmbp->pmhd->peos->eos_data.gamma;
  }
  scaf.gm1 = scaf.gamma - 1.0;

    // return if restart
  if (restart) return;


  auto scaf_init = scaf;
  auto &size = pmbp->pmb->mb_size;
  Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
  const int nmkji = (pmbp->nmb_thispack)*indcs.nx3*indcs.nx2*indcs.nx1;
  const int nkji = indcs.nx3*indcs.nx2*indcs.nx1;
  const int nji = indcs.nx2*indcs.nx1;

  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;

  //Initial Conditions
  par_for("pgen_scaf_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
//   for (int k=ks; k<=ke; k++) {
//     for (int j=js; j<=je; j++) {
//       for (int i=is; i<=ie; i++) {
        // density = pressure = 1

  Real &x1min = size.d_view(m).x1min;
  Real &x1max = size.d_view(m).x1max;
  Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

  Real &x2min = size.d_view(m).x2min;
  Real &x2max = size.d_view(m).x2max;
  Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

  Real &x3min = size.d_view(m).x3min;
  Real &x3max = size.d_view(m).x3max;
  Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

  Real r = std::sqrt(SQR(x1v)+SQR(x2v)+SQR(x3v));
  Real vxgoal = -x2v*std::sqrt(scaf_init.fcirc);
  Real vygoal = x1v*std::sqrt(scaf_init.fcirc);
  Real rho0 = 1.0;//std::max(pow(1 + (gamma-1)/r - (gamma - 1)/2.0*fcirc/(pow(x, 2) + pow(y, 2)), 1/(gamma-1)),dfloor);
      
  w0_(m,IDN,k,j,i) = fmax(rho0, scaf_init.dfloor);
  w0_(m,IVX,k,j,i) = vxgoal;
  w0_(m,IVY,k,j,i) = vygoal;
  w0_(m,IVZ,k,j,i) = 0;
  w0_(m,IEN,k,j,i) = fmax(rho0*scaf_init.T0, scaf_init.pfloor)/scaf_init.gm1;
});

    if (pmbp->pmhd != nullptr) {
    std::cout << "### FATAL ERROR " << std::endl
              << "MHD not implemented for this problem" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Convert primitive to conserved 
  if (pmbp->phydro != nullptr) {
    pmbp->phydro->peos->PrimToCons(w0_, u0_, is, ie, js, je, ks, ke);
  } 
  // actually not used here
  else if (pmbp->pmhd != nullptr) {
    auto &bcc0_ = pmbp->pmhd->bcc0;
    pmbp->pmhd->peos->PrimToCons(w0_, bcc0_, u0_, is, ie, js, je, ks, ke);
  }

  return;

}


//----------------------------------------------------------------------------------------
//! \fn void AddUserSrcs()
//! \brief Add User Source Terms
// NOTE source terms must all be computed using primitive (w0) and NOT conserved (u0) vars
void AddUserSrcs(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  DvceArray5D<Real> &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  DvceArray5D<Real> &w0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  const EOS_Data &eos_data = (pmbp->pmhd != nullptr) ?
                             pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
//   if (scaf_flag) {
    //std::cout << "AddAccel" << std::endl;
    MySource(pm,bdt,u0,w0,eos_data);
//   }
  return;
}


// void MySource(MeshBlock *pmb, const Real time, const Real dt,
//               const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
//               const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
//               AthenaArray<Real> &cons_scalar) {
void MySource(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                   const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {

MeshBlockPack *pmbp = pm->pmb_pack;
auto &indcs = pmbp->pmesh->mb_indcs;
int is = indcs.is, ie = indcs.ie, nx1 = indcs.nx1;
int js = indcs.js, je = indcs.je, nx2 = indcs.nx2;
int ks = indcs.ks, ke = indcs.ke, nx3 = indcs.nx3;
const int nmb = pmbp->nmb_thispack;
auto size = pmbp->pmb->mb_size;
const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
const int nkji = nx3*nx2*nx1;
const int nji = nx2*nx1;
Real beta = bdt/pm->dt;
// Real dt = pm->dt;

// @hyw: here may cause some confusion... But bdt is actually beta[current-stage] * dt, 
// which is the time step for the current stage/(explicit part in IMEX) of the RK scheme.

// todo: check the correction and also check the IMEX scheme in im-eos

Real dt = bdt;

Real time = pm->time;
Real cfl_no = pm->cfl_no;
auto &eos = eos_data;
Real use_e = eos_data.use_e;
Real tfloor = eos_data.tfloor;
auto scaf_temp = scaf;
auto &flx_temp = pmbp->phydro->uflx; 

Real GM = 1.0;
// Real T0 = 1.0; // T0 is now set by input
Real t_cool = 1.0;
//   Real gamma = pmb->peos->GetGamma();
//   Real gm1 = gamma - 1.0;

  par_for("scaf_source", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

    auto scaf_ = scaf_temp;

    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;

    Real x1v_i = CellCenterX(i-is, nx1, x1min, x1max);
    Real x1f_ip1 = LeftEdgeX(i-is+1, nx1, x1min, x1max);
    Real x1f_i = LeftEdgeX(i-is, nx1, x1min, x1max);
    
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;

    Real x2v_j = CellCenterX(j-js, nx2, x2min, x2max);
    Real x2f_jp1 = LeftEdgeX(j-js+1, nx2, x2min, x2max);
    Real x2f_j = LeftEdgeX(j-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    Real x3v_k = CellCenterX(k-ks, nx3, x3min, x3max);
    Real x3f_kp1 = LeftEdgeX(k-ks+1, nx3, x3min, x3max);
    Real x3f_k = LeftEdgeX(k-ks, nx3, x3min, x3max);

    Real r_smooth; // scale below which g is smoothed; also the sink radius
    r_smooth = scaf_.sink_radius_n_cell*(x1f_ip1-x1f_i);
    //   r_smooth = sink_radius_n_cell*pmb->pcoord->dx1f(0);

    // step 1. central point mass
    Real x = x1v_i;
    Real y = x2v_j;
    Real z = x3v_k;

    Real r = std::sqrt(SQR(x)+SQR(y)+SQR(z)); 

    // r for on interfaces
    Real r1l = std::sqrt(SQR(x1f_i  )+SQR(x2v_j)+SQR(x3v_k));
    Real r1r = std::sqrt(SQR(x1f_ip1)+SQR(x2v_j)+SQR(x3v_k));
    Real r2l = std::sqrt(SQR(x1v_i)+SQR(x2f_j)  +SQR(x3v_k));
    Real r2r = std::sqrt(SQR(x1v_i)+SQR(x2f_jp1)+SQR(x3v_k));
    Real r3l = std::sqrt(SQR(x1v_i)+SQR(x2v_j)+SQR(x3f_k ));
    Real r3r = std::sqrt(SQR(x1v_i)+SQR(x2v_j)+SQR(x3f_kp1));
    // phi = -1/r at large radii, phi = const at small radii
    
    Real r_smooth_3 = r_smooth*r_smooth*r_smooth;

    Real phi   = (r>r_smooth)   ? -1./r   : -1./r_smooth + .5*(SQR(r)-SQR(r_smooth))/r_smooth_3;
    Real phi1l = (r1l>r_smooth) ? -1./r1l : -1./r_smooth + .5*(SQR(r1l)-SQR(r_smooth))/r_smooth_3;
    Real phi1r = (r1r>r_smooth) ? -1./r1r : -1./r_smooth + .5*(SQR(r1r)-SQR(r_smooth))/r_smooth_3;
    Real phi2l = (r2l>r_smooth) ? -1./r2l : -1./r_smooth + .5*(SQR(r2l)-SQR(r_smooth))/r_smooth_3;
    Real phi2r = (r2r>r_smooth) ? -1./r2r : -1./r_smooth + .5*(SQR(r2r)-SQR(r_smooth))/r_smooth_3;
    Real phi3l = (r3l>r_smooth) ? -1./r3l : -1./r_smooth + .5*(SQR(r3l)-SQR(r_smooth))/r_smooth_3;
    Real phi3r = (r3r>r_smooth) ? -1./r3r : -1./r_smooth + .5*(SQR(r3r)-SQR(r_smooth))/r_smooth_3;

        // conservative source term from self gravity module
        // AthenaArray<Real> *flux = pmb->phydro->flux;
        // DvceFaceFld5D<Real> uflx = pmbp->phydro->uflx;
        
        // Real dx1 = pmb->pcoord->dx1v(i);
        // Real dx2 = pmb->pcoord->dx2v(j);
        // Real dx3 = pmb->pcoord->dx3v(k);
       
        Real dx1 = CellCenterX(i-is+1, nx1, x1min, x1max) - CellCenterX(i-is, nx1, x1min, x1max);
        Real dx2 = CellCenterX(j-js+1, nx2, x2min, x2max) - CellCenterX(j-js, nx2, x2min, x2max);
        Real dx3 = CellCenterX(k-ks+1, nx3, x3min, x3max) - CellCenterX(k-ks, nx3, x3min, x3max);
        
        Real dtodx1 = dt/dx1;
        Real dtodx2 = dt/dx2;
        Real dtodx3 = dt/dx3;

        auto flx_ = flx_temp;
        auto flx1_ = flx_.x1f;
        auto flx2_ = flx_.x2f;
        auto flx3_ = flx_.x3f;

        // cons(IM1,k,j,i) -= dtodx1*prim(IDN,k,j,i)*(phi1r-phi1l);
        // cons(IM2,k,j,i) -= dtodx2*prim(IDN,k,j,i)*(phi2r-phi2l);
        // cons(IM3,k,j,i) -= dtodx3*prim(IDN,k,j,i)*(phi3r-phi3l);


    // u0(m,IM1,k,j,i) -= x*dt*w0(m,IDN,k,j,i)/(r*r*r);
    // u0(m,IM2,k,j,i) -= y*dt*w0(m,IDN,k,j,i)/(r*r*r);
    // u0(m,IM3,k,j,i) -= z*dt*w0(m,IDN,k,j,i)/(r*r*r);

        u0(m,IM1,k,j,i) -= dtodx1*w0(m,IDN,k,j,i)*(phi1r-phi1l);
        u0(m,IEN,k,j,i) -= dtodx1*(flx1_(m,IDN,k,j,i)*(phi - phi1l) +
                                   flx1_(m,IDN,k,j,i+1)*(phi1r - phi));

        u0(m,IM2,k,j,i) -= dtodx2*w0(m,IDN,k,j,i)*(phi2r-phi2l);
        u0(m,IEN,k,j,i) -= dtodx2*(flx2_(m,IDN,k,j,i)*(phi - phi2l) +
                                   flx2_(m,IDN,k,j+1,i)*(phi2r - phi));

        u0(m,IM3,k,j,i) -= dtodx3*w0(m,IDN,k,j,i)*(phi3r-phi3l);
        u0(m,IEN,k,j,i) -= dtodx3*(flx3_(m,IDN,k,j,i)*(phi - phi3l) +
                                   flx3_(m,IDN,k+1,j,i)*(phi3r - phi));  
        
        // if (NON_BAROTROPIC_EOS) {
        //   cons(IEN,k,j,i) -= dtodx1*(flux[X1DIR](IDN,k,j,i  )*(phi - phi1l) +
        //                              flux[X1DIR](IDN,k,j,i+1)*(phi1r - phi));
        //   cons(IEN,k,j,i) -= dtodx2*(flux[X2DIR](IDN,k,j  ,i)*(phi - phi2l) +
        //                                flux[X2DIR](IDN,k,j+1,i)*(phi2r - phi));
        //   cons(IEN,k,j,i) -= dtodx3*(flux[X3DIR](IDN,k  ,j,i)*(phi - phi3l) +
        //                                flux[X3DIR](IDN,k+1,j,i)*(phi3r - phi));  
        // }

        // now accrete... acc_mode=0 means no accretor.
        // accrete within r_smooth
        if (r<r_smooth) {
          Real dt_dyn = dt * std::pow(r_smooth, -1.5);
          if (scaf_.acc_mode==1) { // remove mass, keep velocity and specific momentum
            u0(m,IDN,k,j,i) -= scaf_.acc_rate*dt_dyn*u0(m,IDN,k,j,i);
            u0(m,IM1,k,j,i) -= scaf_.acc_rate*dt_dyn*u0(m,IM1,k,j,i);
            u0(m,IM2,k,j,i) -= scaf_.acc_rate*dt_dyn*u0(m,IM2,k,j,i);
            u0(m,IM3,k,j,i) -= scaf_.acc_rate*dt_dyn*u0(m,IM3,k,j,i);
            // if (NON_BAROTROPIC_EOS) 
            u0(m,IEN,k,j,i) -= scaf_.acc_rate*dt_dyn*u0(m,IEN,k,j,i);
          }
          else if (scaf_.acc_mode==2) { // enforce floor density + zero velocity
            u0(m,IDN,k,j,i) = 1.0;
            u0(m,IM1,k,j,i) = 0.;
            u0(m,IM2,k,j,i) = 0.;
            u0(m,IM3,k,j,i) = 0.;
            // if (NON_BAROTROPIC_EOS) 
            u0(m,IEN,k,j,i) = 1.0/scaf_.gm1;
          }
          else if (scaf_.acc_mode==-1) { // inject mommentum (but no energy): do this on everything in sink region
            // velocity increases by r/t_dyn over one t_dyn
            Real dt_a = dt*r/r_smooth_3;
            
            u0(m,IM1,k,j,i) += w0(m,IDN,k,j,i)*dt_a*x/r;
            u0(m,IM2,k,j,i) += w0(m,IDN,k,j,i)*dt_a*y/r;
            u0(m,IM3,k,j,i) += w0(m,IDN,k,j,i)*dt_a*z/r;
          }
          else if (scaf_.acc_mode==-2) { // inject mommentum (but no energy): do this only to inflow
            Real dt_a = dt*r/r_smooth_3;
            bool is_inflow = (w0(m,IVX,k,j,i)*x+w0(m,IVY,k,j,i)*y+w0(m,IVZ,k,j,i)*z)<0.;
            if (is_inflow) {
              u0(m,IM1,k,j,i) += w0(m,IDN,k,j,i)*dt_a*x/r;
              u0(m,IM2,k,j,i) += w0(m,IDN,k,j,i)*dt_a*y/r;
              u0(m,IM3,k,j,i) += w0(m,IDN,k,j,i)*dt_a*z/r;
            }
          }
          else if (scaf_.acc_mode==-3) { // inject mommentum (but no energy): do this only to outflow
            Real dt_a = dt*r/r_smooth_3;
            bool is_outflow = (w0(m,IVX,k,j,i)*x+w0(m,IVY,k,j,i)*y+w0(m,IVZ,k,j,i)*z)>0.;
            if (is_outflow) {
              u0(m,IM1,k,j,i) += w0(m,IDN,k,j,i)*dt_a*x/r;
              u0(m,IM2,k,j,i) += w0(m,IDN,k,j,i)*dt_a*y/r;
              u0(m,IM3,k,j,i) += w0(m,IDN,k,j,i)*dt_a*z/r;
            }
          }
        }


        if (r>1.) {
          // step 2. cooling: cool towards T=p/rho=T0
          // if (NON_BAROTROPIC_EOS) {
            if (scaf_.only_damp_v) { // only damp velocity
              Real Ek = .5*w0(m,IDN,k,j,i)*(SQR(w0(m,IVX,k,j,i))+SQR(w0(m,IVY,k,j,i))+SQR(w0(m,IVZ,k,j,i)));
              
              u0(m,IEN,k,j,i) -= Ek * dt/t_cool * 2;
              
              Real P_goal = (w0(m,IDN,k,j,i)/r*scaf_.B0 - Ek)/(scaf_.gamma/scaf_.gm1);
              P_goal = std::max(P_goal, 1.e-8);
              
              u0(m,IEN,k,j,i) += (P_goal-w0(m,IEN,k,j,i)*scaf_.gm1)/scaf_.gm1 * (1.-std::exp(-10.*dt/t_cool));
            } 
            else if (scaf_.set_B0) { // cool towards B0
              Real T_B0 = scaf_.gm1/scaf_.gamma*scaf_.B0/r;
              
              u0(m,IEN,k,j,i) -= dt/t_cool * (w0(m,IEN,k,j,i)*scaf_.gm1-T_B0*w0(m,IDN,k,j,i))/scaf_.gm1;
              
              Real Ek = .5*w0(m,IDN,k,j,i)*(SQR(w0(m,IVX,k,j,i))+SQR(w0(m,IVY,k,j,i))+SQR(w0(m,IVZ,k,j,i)));
              
              u0(m,IEN,k,j,i) -= Ek * dt/t_cool * 2;
            } else { // cool towards T=p/rho=T0
              u0(m,IEN,k,j,i) -= dt/t_cool * (w0(m,IEN,k,j,i)*scaf_.gm1-scaf_.T0*w0(m,IDN,k,j,i))/scaf_.gm1;
            }
          // } // (NON_BAROTROPIC_EOS)


          // step 3. mass removal towards rho = 1.
            Real bracket = std::max(1 + (scaf_.gamma-1)/r - (scaf_.gamma - 1)/2.0*scaf_.fcirc/(pow(x, 2) + pow(y, 2)),scaf_.dfloor);
            Real rhogoal = std::max(pow(bracket, 1/(scaf_.gamma-1)),scaf_.dfloor); //define a goal density
          //Real new_rho_fraction = 1.*(1.-dt/t_cool/rhogoal) + 1./w0(m,IDN,k,j,i)*(dt/t_cool); //This doesn't work
            Real new_rho_fraction = 1.-dt/t_cool*(w0(m,IDN,k,j,i)-rhogoal)/w0(m,IDN,k,j,i); //damping towards a goal density, turns out to be equivalent
          u0(m,IDN,k,j,i) *= new_rho_fraction;
          u0(m,IM1,k,j,i) *= new_rho_fraction;
          u0(m,IM2,k,j,i) *= new_rho_fraction;
          u0(m,IM3,k,j,i) *= new_rho_fraction;
        //   if (NON_BAROTROPIC_EOS)
          u0(m,IEN,k,j,i) *= new_rho_fraction;
            
          // step 4. damp velocity
            Real damp_faction = 1.-dt/t_cool;
            //Real vxgoal = -(y / (pow(x, 2) + pow(y, 2)))*std::sqrt(fcirc); //a goal vx
            //Real vygoal = (x / (pow(x, 2) + pow(y, 2)))*std::sqrt(fcirc); //a goal vy
            Real vxgoal = -y*std::sqrt(scaf_.fcirc); //a goal vx
            Real vygoal = x*std::sqrt(scaf_.fcirc); //a goal vy
            // above: changed to uniform rotation
            u0(m,IM1,k,j,i) += w0(m,IDN,k,j,i)*dt/t_cool*(vxgoal-w0(m,IVX,k,j,i)) ;
            u0(m,IM2,k,j,i) += w0(m,IDN,k,j,i)*dt/t_cool*(vygoal-w0(m,IVY,k,j,i)) ;
            u0(m,IM3,k,j,i) *= 1.-dt/t_cool;
            
            //Real signx = std::abs(w0(m,IVX,k,j,i)*vxgoal)/(w0(m,IVX,k,j,i)*vxgoal);
            //Real signy = std::abs(w0(m,IVY,k,j,i)*vygoal)/(w0(m,IVY,k,j,i)*vygoal);
            
          // step 5. force b=0 //Bernoulli constant, optional
            
          /*if (only_damp_v) {
            Real ek = .5/SQR(u0(m,IDN,k,j,i))*(
                      SQR(u0(m,IM1,k,j,i))+
                      SQR(u0(m,IM2,k,j,i))+
                      SQR(u0(m,IM3,k,j,i)));
            Real p = std::max(1.e-8,(1./r - ek)/scaf_.gamma*scaf_.gm1);
            u0(m,IEN,k,j,i) = p/scaf_.gm1 + u0(m,IDN,k,j,i)*ek;
          }*/
        }
        // step -1. driving turbulence: not used now
        
         
  });

    return;
}



// refinement condition: check the maximum pressure gradient
// int RefinementCondition(MeshBlock *pmb) {
void RefinementCondition(MeshBlockPack* pmbp) {
  // capture variables for kernels
  Mesh *pm = pmbp->pmesh;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie, nx1 = indcs.nx1;
  int js = indcs.js, je = indcs.je, nx2 = indcs.nx2;
  int ks = indcs.ks, ke = indcs.ke, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;

  // check (on device) Hydro/MHD refinement conditions over all MeshBlocks
  auto refine_flag_ = pm->pmr->refine_flag;

  int nmb = pmbp->nmb_thispack;
  int mbs = pm->gids_eachrank[global_variable::my_rank];

  par_for_outer("ConsRefineCond",DevExeSpace(), 0, 0, 0, (nmb-1),
  KOKKOS_LAMBDA(TeamMember_t tmember, const int m) {
  // original Athena++ refinement condition:
  // Real P1 = pmb->pcoord->x1v(pmb->is-1) * pmb->pcoord->x1v(pmb->ie+1);
  // Real P2 = pmb->pcoord->x2v(pmb->js-1) * pmb->pcoord->x2v(pmb->je+1);
  // Real P3 = pmb->pcoord->x3v(pmb->ks-1) * pmb->pcoord->x3v(pmb->ke+1);
  // if (P1<0. && P2<0. && P3<0.) return 1;

  Real &x1min = size.d_view(m).x1min;
  Real &x1max = size.d_view(m).x1max;
  // Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
  
  Real &x2min = size.d_view(m).x2min;
  Real &x2max = size.d_view(m).x2max;
  // Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

  Real &x3min = size.d_view(m).x3min;
  Real &x3max = size.d_view(m).x3max;
  // Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

  Real ax1min = x1min*x1max>0.0? fmin(fabs(x1min), fabs(x1max)) : 0.0;
  Real ax2min = x2min*x2max>0.0? fmin(fabs(x2min), fabs(x2max)) : 0.0;
  Real ax3min = x3min*x3max>0.0? fmin(fabs(x3min), fabs(x3max)) : 0.0;
  Real rad_min = sqrt(SQR(ax1min)+SQR(ax2min)+SQR(ax3min));

  Real P1 = CellCenterX(is-2, nx1, x1min, x1max) * CellCenterX(ie+2, nx1, x1min, x1max);
  Real P2 = CellCenterX(js-2, nx2, x2min, x2max) * CellCenterX(je+2, nx2, x2min, x2max);
  Real P3 = CellCenterX(ks-2, nx3, x3min, x3max) * CellCenterX(ke+2, nx3, x3min, x3max);

// std::cout << "P1: " << P1 << " P2: " << P2 << " P3: " << P3 << std::endl;
  // if (P1<0. && P2<0. && P3<0.) refine_flag_.d_view(m+mbs)=1;   
  // else refine_flag_.d_view(m+mbs)=0;   

  if (rad_min<0.005) refine_flag_.d_view(m+mbs)=1;   
  else refine_flag_.d_view(m+mbs)=0;   

  // if (fabs(x1v-r_refine) && x2v<0. && x3v<0.) refine_flag_.d_view(m)=1;   
  // else refine_flag_.d_view(m)=0;   

  });
  


}

// hst outputs
void HstOutput(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &eos_data = pm->pmb_pack->phydro->peos->eos_data;
  int &nhydro_ = pm->pmb_pack->phydro->nhydro;
  
  int nhistvar = 13;
  pdata->nhist = 13;

  if (pdata->nhist > NHISTORY_VARIABLES) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "User history function specified pdata->nhist larger than"
              << " NHISTORY_VARIABLES" << std::endl;
    exit(EXIT_FAILURE);
  }
  pdata->label[0] = "mdot_";
  pdata->label[1] = "rho__";
  pdata->label[2] = "Lz___";
  pdata->label[3] = "Ly___";
  pdata->label[4] = "Lx___";
  pdata->label[5] = "Ltot_";
  pdata->label[6] = "Lz2__";
  pdata->label[7] = "Ly2__";
  pdata->label[8] = "Lx2__";
  pdata->label[9] = "lz___";
  pdata->label[10] = "ly___";
  pdata->label[11] = "lx___";
  

  for (int n=0; n<nhistvar; ++n) {
    pdata->hdata[n] = 0.0;
  }

  // capture class variabels for kernel
  auto &u0_ = pm->pmb_pack->phydro->u0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;

  // extract grids, number of radii, number of fluxes, and history appending index

  // set number of and names of history variables for hydro or mhd
  //  (1) accretion rate, then both Mdot and rho_center are displaying accretion rate
  //  (2) mean density
  //  (3) total Lz

  auto scaf_temp = scaf;
  Real A = 0.;

  // for (int n=0; n<NREDUCTION_VARIABLES; ++n) {
  //   sum_this_mb.the_array[n] = 0.0;
  // }

  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    auto scaf_ = scaf_temp;

    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    Real x = x1v;
    Real y = x2v;
    Real z = x3v;
    Real r = std::sqrt(SQR(x)+SQR(y)+SQR(z));
    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

    Real dV = vol;
    Real r_sink = scaf_.sink_radius_n_cell*size.d_view(m).dx1;
    Real OmegaK = std::pow(r_sink, -1.5);

    // Hydro conserved variables:
    array_sum::GlobalSum hvars;

    Real den = u0_(m,IDN,k,j,i);
    Real px = u0_(m,IM1,k,j,i);
    Real py = u0_(m,IM2,k,j,i);
    Real pz = u0_(m,IM3,k,j,i);
    Real Ek = .5*den*(SQR(px)+SQR(py)+SQR(pz));

    // (1) accretion rate, then both Mdot and rho_center are displaying accretion rate
    if (scaf_.acc_mode==1 && r<r_sink)
      hvars.the_array[0] = dV*u0_(m,IDN,k,j,i)*OmegaK*scaf_.acc_rate;
    else if (scaf_.acc_mode==2 && r>r_sink && r<1.5*r_sink)
      hvars.the_array[0] = -M_PI*r*
                           (u0_(m,IM1,k,j,i)*x + 
                            u0_(m,IM2,k,j,i)*y + 
                            u0_(m,IM3,k,j,i)*z) / scaf_.n_cell_acc2;

    // (2) mean density

    if (scaf_.acc_mode==1 && r<r_sink)
      hvars.the_array[1] = u0_(m,IDN,k,j,i) / scaf_.n_cell_acc1;
    else if (scaf_.acc_mode==2 && r>r_sink && r<1.5*r_sink)
      hvars.the_array[1] = u0_(m,IDN,k,j,i) / scaf_.n_cell_acc2;

    // (3) total Lz
    if (scaf_.acc_mode==1 && r<r_sink)
      hvars.the_array[2] = dV*OmegaK*scaf_.acc_rate*(x*u0_(m,IM2,k,j,i) - y*u0_(m,IM1,k,j,i));

    // (4) total Ly
    if (scaf_.acc_mode==1 && r<r_sink)
      hvars.the_array[3] = dV*OmegaK*scaf_.acc_rate*(z*u0_(m,IM1,k,j,i) - x*u0_(m,IM3,k,j,i));

    // (5) total Lx
    if (scaf_.acc_mode==1 && r<r_sink)
      hvars.the_array[4] = dV*OmegaK*scaf_.acc_rate*(y*u0_(m,IM3,k,j,i) - z*u0_(m,IM2,k,j,i));

    // (6) total angular momentum (at the sink cell)
    if (scaf_.acc_mode==1 && r<r_sink)
      hvars.the_array[5] = dV*OmegaK*scaf_.acc_rate*(SQR(x)*u0_(m,IM1,k,j,i) + SQR(y)*u0_(m,IM2,k,j,i) + SQR(z)*u0_(m,IM3,k,j,i));

    // (7,8,9) total angular momentum (in the domain)
    // @hyw: todo: different radial bins
    // current version: 10 times the radius of the sink cell 
    if (scaf_.acc_mode==1 && r>=r_sink && r<10*r_sink) {
      hvars.the_array[6] = dV*(x*py - y*px); // Lz
      hvars.the_array[7] = dV*(z*px - x*pz); // Ly
      hvars.the_array[8] = dV*(y*pz - z*py); // Lx
    }

    // (10,11,12) specific angular momentum (total/mass)
    if (scaf_.acc_mode==1 && r>=r_sink && r<10*r_sink) {
      hvars.the_array[9] = (x*py - y*px)/den; // Lz
      hvars.the_array[10] = (z*px - x*pz)/den; // Ly
      hvars.the_array[11] = (y*pz - z*py)/den; // Lx
    }

    // fill rest of the_array with zeros, if nhist < NREDUCTION_VARIABLES
    for (int n=nhist_; n<NREDUCTION_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }
  // for (int n=pdata->nhist; n<NREDUCTION_VARIABLES; ++n) {
  //   pdata->hdata[n] = 0.0;
  // }
  return;

}

