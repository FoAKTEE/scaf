//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ion-eos.cpp
//  \brief
// created: 2024.3.
// @hyw IMEX sink (done) & cooling (possibly in the todo list)

#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "im_eos.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"

namespace im_eos {
//----------------------------------------------------------------------------------------
// constructor, parses input file and initializes data structures and parameters

ImEos::ImEos(MeshBlockPack *pp, ParameterInput *pin) : 
  pmy_pack(pp) {

  // @hyw: basically copy from the pgen-scaf
  
  // Read various coefficients
  params.acc_mode=1; // flag for accretor setup. 0=no accretor.
  params.acc_rate=1.0; // accretion efficiency for smooth accretor
  params.T0=1.; // T at r>1
  params.B0=1.; // ((gamma-1)/gamma * P + Ek) / Phi
  params.set_B0=false; // whether set T0 or B0
  params.only_damp_v=false; // damp velocity only
  params.output_Mdot=false;
  params.fcirc=0.0;
  params.sink_radius_n_cell=4; // number of cells for sink radius
  params.n_cell_acc1 = 0; // cells at r<r_sink
  params.n_cell_acc2 = 0; // cells at r_sink - 1.5 r_sink
  params.turb_flag=0; // flag for turbulence

  // void Mesh::InitUserMeshData(ParameterInput *pin) in Athena++
  params.acc_mode = pin->GetOrAddInteger("problem","acc_mode",params.acc_mode);
  params.acc_rate = pin->GetOrAddReal("problem","acc_rate",params.acc_rate);
  params.T0 = pin->GetOrAddReal("problem","T0",params.T0);
  params.B0 = pin->GetOrAddReal("problem","B0",params.B0);
  params.set_B0 = pin->GetOrAddBoolean("problem","set_B0",params.set_B0);
  params.only_damp_v = pin->GetOrAddBoolean("problem","only_damp_v",params.only_damp_v);
  params.output_Mdot = pin->GetOrAddBoolean("problem","output_Mdot",params.output_Mdot);
  params.turb_flag = pin->GetInteger("problem","turb_flag");
  params.fcirc = pin->GetOrAddReal("problem","fcirc",params.fcirc);
  params.sink_radius_n_cell = pin->GetOrAddReal("problem","sink_radius_n_cell",params.sink_radius_n_cell);
  params.n_cell_acc1 = n_cells_between_radii(0, params.sink_radius_n_cell);
  params.n_cell_acc2 = n_cells_between_radii(params.sink_radius_n_cell, 1.5*params.sink_radius_n_cell);

//! void MeshBlock::ProblemGenerator in Athena++

  params.dfloor=pin->GetOrAddReal("hydro","dfloor",0.01);
  params.pfloor=pin->GetOrAddReal("hydro","pfloor",0.01);

  params.gamma = pp->phydro->peos->eos_data.gamma;
  params.gm1 = params.gamma - 1.0;
  
}

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

//----------------------------------------------------------------------------------------
//! \fn  void IonNeutral::AssembleIonNeutralTasks
//  \brief Adds tasks for ion-neutral (two-fluid) mhd to stage start/run/end task lists
//  Called by MeshBlockPack::AddPhysics() function directly after MHD constrctr

void ImEos::AssembleImEosTasks(TaskList &start, TaskList &run, TaskList &end) {
  TaskID none(0);

  using namespace hydro;  // NOLINT(build/namespaces)
  Hydro *phyd = pmy_pack->phydro;


  id.n_irecv = start.AddTask(&Hydro::InitRecv, phyd, none);
  // assemble run task list
  id.impl_2x = run.AddTask(&ImEos::FirstTwoImpRKImEos, this, none);  // does CopyCons

  // id.i_flux  = run.AddTask(&MHD::Fluxes, pmhd, id.impl_2x);
  // id.i_sendf = run.AddTask(&MHD::SendFlux, pmhd, id.i_flux);
  // id.i_recvf = run.AddTask(&MHD::RecvFlux, pmhd, id.i_sendf);
  // id.i_expl  = run.AddTask(&MHD::ExpRKUpdate, pmhd, id.i_recvf);

  id.n_flux = run.AddTask(&Hydro::Fluxes, phyd, id.impl_2x);
  id.n_sendf = run.AddTask(&Hydro::SendFlux, phyd, id.n_flux);
  id.n_recvf = run.AddTask(&Hydro::RecvFlux, phyd, id.n_sendf);
  id.n_expl = run.AddTask(&Hydro::ExpRKUpdate, phyd, id.n_recvf);

  id.impl = run.AddTask(&ImEos::ImpRKUpdateImEos, this, id.n_expl);
  // id.i_restu = run.AddTask(&MHD::RestrictU, pmhd, id.impl);
  id.n_restu = run.AddTask(&Hydro::RestrictU, phyd, id.impl);

  // id.i_sendu = run.AddTask(&MHD::SendU, pmhd, id.n_restu);
  id.n_sendu = run.AddTask(&Hydro::SendU, phyd, id.n_restu);
  // id.i_recvu = run.AddTask(&MHD::RecvU, pmhd, id.i_sendu);
  id.n_recvu = run.AddTask(&Hydro::RecvU, phyd, id.n_sendu);

  // id.efld  = run.AddTask(&MHD::CornerE, pmhd, id.i_recvu);
  // id.ct    = run.AddTask(&MHD::CT, pmhd, id.efld);
  // id.restb = run.AddTask(&MHD::RestrictB, pmhd, id.ct);
  // id.sendb = run.AddTask(&MHD::SendB, pmhd, id.restb);
  // id.recvb  = run.AddTask(&MHD::RecvB, pmhd, id.sendb);

  // id.i_bcs   = run.AddTask(&MHD::ApplyPhysicalBCs, pmhd, id.recvb);
  id.n_bcs   = run.AddTask(&Hydro::ApplyPhysicalBCs, phyd, id.n_recvu);
  // id.i_prol  = run.AddTask(&MHD::Prolongate, pmhd, id.i_bcs);
  id.n_prol  = run.AddTask(&Hydro::Prolongate, phyd, id.n_bcs);
  // id.i_c2p   = run.AddTask(&MHD::ConToPrim, pmhd, id.i_prol);
  id.n_c2p   = run.AddTask(&Hydro::ConToPrim, phyd, id.n_prol);
  // id.i_newdt = run.AddTask(&MHD::NewTimeStep, pmhd, id.i_c2p);
  id.n_newdt = run.AddTask(&Hydro::NewTimeStep, phyd, id.n_c2p);

  // assemble end task list
  // id.i_clear = end.AddTask(&MHD::ClearSend, pmhd, none);
  id.n_clear = end.AddTask(&Hydro::ClearSend, phyd, none);
  

  return;
}
//----------------------------------------------------------------------------------------
//! \fn ImEos::FirstTwoImpRK
//  \brief Executes first two implicit stages of the ImEx integrator for ion-neutral
//  drag term.  Should be the first task called in TaskList.

TaskStatus ImEos::FirstTwoImpRKImEos(Driver *pdrive, int stage) {
  if (stage != 1) {return TaskStatus::complete;}  // only execute on first stage

  // copy conserved hydro and MHD variables
  hydro::Hydro *phyd = pmy_pack->phydro;
  Kokkos::deep_copy(DevExeSpace(), phyd->u1, phyd->u0);

  // Solve implicit equations first time (nexp_stage = -1)
  auto status = ImpRKUpdateImEos(pdrive, -1);

  // Solve implicit equations second time (nexp_stage = 0)
  status = ImpRKUpdateImEos(pdrive, 0);

  // update primitive variables for both hydro and MHD
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  
  phyd->peos->ConsToPrim(phyd->u0, phyd->w0, false, 0, n1m1, 0, n2m1, 0, n3m1);
  
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn  void IonNeutral::ImpRKUpdate
//  \brief Implicit RK update of ion-neutral drag term. Used as part of ImEx RK integrator
//  This function should be added AFTER the explicit updates in the task list, so that
//  source terms are evaluated using partially updated values (including explicit terms
//  such as flux divergence).  This means soure terms must only be evaluated using
//  conserved variables (u0), as primitives (w0) are not updated until end of TaskList.
//
//  Note indices of source term array correspond to:
    // ru(0) -> u(IM1)
    // ru(1) -> u(IM2)
    // ru(2) -> u(IM3)
    // ru(3) -> u(IDN)
    // ru(4) -> u(IEN)
    // @hyw: all HYDRO ***



TaskStatus ImEos::ImpRKUpdateImEos(Driver *pdriver, int estage) {
  // # of implicit stage (1,2,3,4,[5]).  Note estage=(# of explicit stage)=(1,2,[3])
  // estage <= 0 corresponds to first two fully implicit stages
  int istage = estage + 2;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  // @hyw
  auto &size = pmy_pack->pmb->mb_size;
  int is = indcs.is, ie = indcs.ie, nx1 = indcs.nx1;
  int js = indcs.js, je = indcs.je, nx2 = indcs.nx2;
  int ks = indcs.ks, ke = indcs.ke, nx3 = indcs.nx3;

  // Add stiff source term (ion-neutral drag) evaluated with values from previous stages,
  // i.e. the R(U^1), R(U^2), etc. terms, to partially updated conserved variables.
  // Only required for istage = (2,3,4,[5])

  hydro::Hydro *phyd = pmy_pack->phydro;
  if (istage > 1) {
    int scr_level = 0;
    size_t scr_size = 0;

    auto u = phyd->u0;

    auto &a_twid = pdriver->a_twid;
    Real dt = pmy_pack->pmesh->dt;
    auto ru_ = pdriver->impl_src_imeos;

    auto params_temp = params;

    par_for_outer("imex_exp",DevExeSpace(),scr_size,scr_level,0,nmb1,0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
      for (int s=0; s<=(istage-2); ++s) {
        Real adt = a_twid[istage-2][s]*dt;
        par_for_inner(member, 0, (n1-1), [&](const int i) {

          auto params_ = params_temp;     
                 
          if (params.acc_mode==2) {
            u(m,IDN,k,j,i) = 1.0;
            u(m,IM1,k,j,i) = 0.;
            u(m,IM2,k,j,i) = 0.;
            u(m,IM3,k,j,i) = 0.;
            u(m,IEN,k,j,i) = 1.0/params_.gm1;
          } 
          else {
            u(m,IDN,k,j,i) += adt*ru_(s,m,0,k,j,i);
            u(m,IM1,k,j,i) += adt*ru_(s,m,1,k,j,i);
            u(m,IM2,k,j,i) += adt*ru_(s,m,2,k,j,i);
            u(m,IM3,k,j,i) += adt*ru_(s,m,3,k,j,i);
            u(m,IEN,k,j,i) += adt*ru_(s,m,4,k,j,i);
          }
          


        });
      }
    });
  }



  // Update ion/neutral momentum equations with analytic solution of implicit difference
  // equations for ion-neutral drag.
  // Only required for istage = (1,2,3,[4])

  if (estage < pdriver->nexp_stages) {

auto u = phyd->u0;
auto params_temp = params;

Real time =pmy_pack->pmesh->time;
Real dt = pmy_pack->pmesh->dt;

par_for("imex_imp",DevExeSpace(),0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

auto params_ = params_temp;


Real &x1min = size.d_view(m).x1min;
Real &x1max = size.d_view(m).x1max;
Real x1v = CellCenterX(i, nx1, x1min, x1max);

Real &x2min = size.d_view(m).x2min;
Real &x2max = size.d_view(m).x2max;
Real x2v = CellCenterX(j, nx2, x2min, x2max);

Real &x3min = size.d_view(m).x3min;
Real &x3max = size.d_view(m).x3max;
Real x3v = CellCenterX(k, nx3, x3min, x3max);

Real x1v_i = CellCenterX(i-is, nx1, x1min, x1max);
Real x1f_ip1 = LeftEdgeX(i-is+1, nx1, x1min, x1max);
Real x1f_i = LeftEdgeX(i-is, nx1, x1min, x1max);

Real r = std::sqrt(x1v*x1v + x2v*x2v + x3v*x3v);
Real x = x1v;
Real y = x2v;
Real z = x3v;

Real denom;
Real r_smooth = params_.sink_radius_n_cell*(x1f_ip1-x1f_i);
Real r_smooth_3 = r_smooth*r_smooth*r_smooth;

Real rhom = u(m,IDN,k,j,i);
Real uux = u(m,IM1,k,j,i);
Real uuy = u(m,IM2,k,j,i);
Real uuz = u(m,IM3,k,j,i);
Real etot = u(m,IEN,k,j,i);

Real vx = u(m,IM1,k,j,i)/u(m,IDN,k,j,i);
Real vy = u(m,IM2,k,j,i)/u(m,IDN,k,j,i);
Real vz = u(m,IM3,k,j,i)/u(m,IDN,k,j,i);


// @hyw
// start doing the sink (inflow; outflow etc)

if (r<r_smooth) {

Real dt_dyn = dt * std::pow(r_smooth, -1.5);
Real dt_a = dt*r/r_smooth_3;
  
if (params_.acc_mode==1) { // remove mass, keep velocity and specific momentum

denom = 1.0 + params_.acc_rate*dt_dyn;

u(m,IDN,k,j,i) = u(m,IDN,k,j,i)/denom;
u(m,IM1,k,j,i) = u(m,IM1,k,j,i)/denom;
u(m,IM2,k,j,i) = u(m,IM2,k,j,i)/denom;
u(m,IM3,k,j,i) = u(m,IM3,k,j,i)/denom;
u(m,IEN,k,j,i) = u(m,IEN,k,j,i)/denom;
  
}
  
else if (params_.acc_mode==2) { // enforce floor density + zero velocity
    
u(m,IDN,k,j,i) = 1.0;
u(m,IM1,k,j,i) = 0.;
u(m,IM2,k,j,i) = 0.;
u(m,IM3,k,j,i) = 0.;
u(m,IEN,k,j,i) = 1.0/params_.gm1;
  
}
  
else if (params_.acc_mode==-1) { // inject mommentum (but no energy): do this on everything in sink region

// velocity increases by r/t_dyn over one t_dyn

u(m,IM1,k,j,i) += u(m,IDN,k,j,i)*dt_a*x/r;
u(m,IM2,k,j,i) += u(m,IDN,k,j,i)*dt_a*y/r;
u(m,IM3,k,j,i) += u(m,IDN,k,j,i)*dt_a*z/r;
  
}
  
else if (params_.acc_mode==-2) { // inject mommentum (but no energy): do this only to inflow

bool is_inflow = (vx*x+vy*y+vz*z)<0.;

if (is_inflow) {
  u(m,IM1,k,j,i) += u(m,IDN,k,j,i)*dt_a*x/r;
  u(m,IM2,k,j,i) += u(m,IDN,k,j,i)*dt_a*y/r;
  u(m,IM3,k,j,i) += u(m,IDN,k,j,i)*dt_a*z/r;
}

}

else if (params_.acc_mode==-3) { // inject mommentum (but no energy): do this only to outflow

bool is_outflow = (vx*x+vy*y+vz*z)>0.;

if (is_outflow) {
  u(m,IM1,k,j,i) += u(m,IDN,k,j,i)*dt_a*x/r;
  u(m,IM2,k,j,i) += u(m,IDN,k,j,i)*dt_a*y/r;
  u(m,IM3,k,j,i) += u(m,IDN,k,j,i)*dt_a*z/r;
}

}
}

    });
  }

  // Compute stiff source term (ion-neutral drag) using variables updated in this stage,
  // i.e R(U^n), for use in later stages.  Only required for istage = (1,2,3,[4])
  if (estage < pdriver->nexp_stages) {

auto params_temp = params;
// @hyw
Real time = pmy_pack->pmesh->time;
Real dt = pmy_pack->pmesh->dt;

int s = istage-1;
auto u = phyd->u0;
auto ru_ = pdriver->impl_src_imeos;

par_for("imex_rup",DevExeSpace(),0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

auto params_ = params_temp;
// @hyw

Real &x1min = size.d_view(m).x1min;
Real &x1max = size.d_view(m).x1max;
Real x1v = CellCenterX(i, nx1, x1min, x1max);

Real &x2min = size.d_view(m).x2min;
Real &x2max = size.d_view(m).x2max;
Real x2v = CellCenterX(j, nx2, x2min, x2max);

Real &x3min = size.d_view(m).x3min;
Real &x3max = size.d_view(m).x3max;
Real x3v = CellCenterX(k, nx3, x3min, x3max);

Real x1v_i = CellCenterX(i-is, nx1, x1min, x1max);
Real x1f_ip1 = LeftEdgeX(i-is+1, nx1, x1min, x1max);
Real x1f_i = LeftEdgeX(i-is, nx1, x1min, x1max);

Real r = std::sqrt(x1v*x1v + x2v*x2v + x3v*x3v);
Real x = x1v;
Real y = x2v;
Real z = x3v;

Real sinkcoeff;
Real r_smooth = params_.sink_radius_n_cell*(x1f_ip1-x1f_i);
Real r_smooth_3 = r_smooth*r_smooth*r_smooth;

Real rhom = u(m,IDN,k,j,i);
Real uux = u(m,IM1,k,j,i);
Real uuy = u(m,IM2,k,j,i);
Real uuz = u(m,IM3,k,j,i);
Real etot = u(m,IEN,k,j,i);

Real vx = u(m,IM1,k,j,i)/u(m,IDN,k,j,i);
Real vy = u(m,IM2,k,j,i)/u(m,IDN,k,j,i);
Real vz = u(m,IM3,k,j,i)/u(m,IDN,k,j,i);

// @hyw
// adding the **implicit** source term

if (r<r_smooth) {

Real dt_dyn = dt * std::pow(r_smooth, -1.5);
Real dt_a = dt*r/r_smooth_3;
  
if (params_.acc_mode==1) { // remove mass, keep velocity and specific momentum

sinkcoeff = params_.acc_rate*dt_dyn;

ru_(s,m,0,k,j,i) = -u(m,IDN,k,j,i)*sinkcoeff;
ru_(s,m,1,k,j,i) = -u(m,IM1,k,j,i)*sinkcoeff;
ru_(s,m,2,k,j,i) = -u(m,IM2,k,j,i)*sinkcoeff;
ru_(s,m,3,k,j,i) = -u(m,IM3,k,j,i)*sinkcoeff;
ru_(s,m,4,k,j,i) = -u(m,IEN,k,j,i)*sinkcoeff;
  
}
  
else if (params_.acc_mode==2) { // enforce floor density + zero velocity
    
ru_(s,m,0,k,j,i) = -u(m,IDN,k,j,i)*sinkcoeff;
ru_(s,m,1,k,j,i) = -u(m,IM1,k,j,i)*sinkcoeff;
ru_(s,m,2,k,j,i) = -u(m,IM2,k,j,i)*sinkcoeff;
ru_(s,m,3,k,j,i) = -u(m,IM3,k,j,i)*sinkcoeff;
ru_(s,m,4,k,j,i) = -u(m,IEN,k,j,i)*sinkcoeff;
  
// u(m,IDN,k,j,i) = 1.0;
// u(m,IM1,k,j,i) = 0.;
// u(m,IM2,k,j,i) = 0.;
// u(m,IM3,k,j,i) = 0.;
// u(m,IEN,k,j,i) = 1.0/params_.gm1;
  
}
  
else if (params_.acc_mode==-1) { // inject mommentum (but no energy): do this on everything in sink region

// velocity increases by r/t_dyn over one t_dyn

sinkcoeff = params_.acc_rate*dt_dyn;

ru_(s,m,0,k,j,i) = 0;
ru_(s,m,1,k,j,i) = u(m,IDN,k,j,i)*dt_a*x/r;
ru_(s,m,2,k,j,i) = u(m,IDN,k,j,i)*dt_a*y/r;
ru_(s,m,3,k,j,i) = u(m,IDN,k,j,i)*dt_a*z/r;
ru_(s,m,4,k,j,i) = 0;
  
}
  
else if (params_.acc_mode==-2) { // inject mommentum (but no energy): do this only to inflow

bool is_inflow = (vx*x+vy*y+vz*z)<0.;

if (is_inflow) {

ru_(s,m,0,k,j,i) = 0;
ru_(s,m,1,k,j,i) = u(m,IDN,k,j,i)*dt_a*x/r;
ru_(s,m,2,k,j,i) = u(m,IDN,k,j,i)*dt_a*y/r;
ru_(s,m,3,k,j,i) = u(m,IDN,k,j,i)*dt_a*z/r;
ru_(s,m,4,k,j,i) = 0;

}

}

else if (params_.acc_mode==-3) { // inject mommentum (but no energy): do this only to outflow

bool is_outflow = (vx*x+vy*y+vz*z)>0.;

if (is_outflow) {

ru_(s,m,0,k,j,i) = 0;
ru_(s,m,1,k,j,i) = u(m,IDN,k,j,i)*dt_a*x/r;
ru_(s,m,2,k,j,i) = u(m,IDN,k,j,i)*dt_a*y/r;
ru_(s,m,3,k,j,i) = u(m,IDN,k,j,i)*dt_a*z/r;
ru_(s,m,4,k,j,i) = 0;

}

}
}


 


  });
  }

  return TaskStatus::complete;
}

} // namespace ion_neutral
