#ifndef OPTPARAM_H_
#define OPTPARAM_H_

#include "model.hpp"
#include "MCseq.hpp"
#include "PartFunc.hpp"
#include "AnaTraj.hpp"

struct OptParam {
  int Length;
  double alpha_H,alpha_J,alpha_K;
  double beta_H,beta_J, beta_K;
  double eps_H, eps_J, eps_K;
  int nsweeps, nsweeps_eq;
  bool use_WL;
  void set_default_alpha(const Model& model);
  double update_H(Model& model);
  double update_J(Model& model);
  double update_K(Model& model);
  void update_all(Model& model, const int iter, const int nsamples);
  double update_iala(Model& model,const int iala);
  
  void optimize(Model& model, MCseq& mcs, const int niter);
  void optimize_mul(Model& model, MCseq& mcs, const int niter);
  void WLoptimize(Model& model, MCseq& mcs, const int niter);
  void optimize_ala(Model& model, MCseq& mcs, const int iala, const int niter);

  void optimize_plm(Model& model, const int niter, const string param_out);

  void optimize_traj(Model& model, AnaTraj& atraj);
  
  // "momentum" terms;
  darray2 H_O, H_I;
  darray3 J_OO, J_OI, J_IO, J_II; // short-range coupling
  darray4 K_OO; // long-range coupling

  void write_momentum(Model& model, const string filename);
  void read_momentum(Model& model, const string filename);
  OptParam(const Model& model);
};

#endif  //OPTPARAM_H_
