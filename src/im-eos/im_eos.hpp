#ifndef IM_EOS_IM_EOS_HPP_
#define IM_EOS_IM_EOS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ion-neutral.hpp
//  \brief definitions for IonNeutral class

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "driver/driver.hpp"

//----------------------------------------------------------------------------------------
//! \struct IonNeutralTaskIDs
//  \brief container to hold TaskIDs of all ion-neutral tasks

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

  // F = - gamma rho_i rho_n (u_i - u_n) + xi rho_n u_n - alpha rho_i^2 u_i
  // G = xi rho_n - alpha rho_i^2
  // @hyw: add sink and cooling related params (for the implicit solvers)
  struct imexparam {

  } params;

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
