#ifndef IM_EOS_IM_EOS_HPP_
#define IM_EOS_IM_EOS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file im-eos.hpp
//  \brief definitions for ImEos class

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "driver/driver.hpp"

//----------------------------------------------------------------------------------------
//! \struct ImEosTaskIDs
//  \brief container to hold TaskIDs of all im-eos tasks

struct ImEosTaskIDs {
  TaskID i_irecv;
  TaskID n_irecv;
  TaskID impl_2x;
  TaskID i_flux;
  TaskID i_sendf;
  TaskID i_recvf;
  TaskID i_expl;
  TaskID i_restu;
  TaskID n_flux;
  TaskID n_sendf;
  TaskID n_recvf;
  TaskID n_expl;
  TaskID n_restu;
  TaskID impl;
  TaskID i_sendu;
  TaskID i_recvu;
  TaskID n_sendu;
  TaskID n_recvu;
  TaskID efld;
  TaskID ct;
  TaskID restb;
  TaskID sendb;
  TaskID recvb;
  TaskID i_bcs;
  TaskID n_bcs;
  TaskID i_prol;
  TaskID n_prol;
  TaskID i_c2p;
  TaskID n_c2p;
  TaskID i_newdt;
  TaskID n_newdt;
  TaskID i_clear;
  TaskID n_clear;
};

namespace im_eos {

//----------------------------------------------------------------------------------------
//! \class ImEos

class ImEos {
 public:
  ImEos(MeshBlockPack *ppack, ParameterInput *pin);
  ~ImEos();

  // @hyw: add sink and cooling related params (for the implicit solvers)
  struct imexparam {

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

  } params;

  Real n_cells_between_radii(Real n_in, Real n_out);

  // container to hold names of TaskIDs
  ImEosTaskIDs id;

  // functions
  void AssembleImEosTasks(TaskList &start, TaskList &run, TaskList &end);
  TaskStatus FirstTwoImpRKImEos(Driver* pdrive, int stage);
  TaskStatus ImpRKUpdateImEos(Driver* pdrive, int stage);

 private:
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Hydro
};

} // namespace im_eos
#endif // IM_EOS_IM_EOS_HPP_
