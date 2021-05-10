#ifndef _PARTFUNC_H_
#define _PARTFUNC_H_

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "model.hpp"
#include "InvertMatrix.hpp"

struct PartFunc {
  PartFunc(Model& model);

  int model_length;
  double beta;
  darray1 lambda,looplength;
  darray2 Kmf_O; // "frozen" fields
  darray2 n0_O; // density of a natural sequence with pseudocounts
  darray3 T_OO, T_OI, T_IO, T_II; // transfer matrices.
  darray3 IT_II; // (I - T_II)^{-1}
  darray2 specT_II; // spectral radii of T_II matrices.
  darray3 U; // integrated transfer matrices.
  darray2 lnZf, lnZb; 
  darray2 lnZf_I, lnZb_I;
  darray2 en_O,en_I;
  darray3 ens_OO, ens_OI, ens_IO, ens_II;
  darray4 enl_OO;
  map< int, vector<int> > IntList; // interactions list for fast look-up.

  void init_expected_number_density(Model& model);
  void set_expected_number_density(void);
  void set_expected_number_density_nocorr(void);
  
  void set_expected_number_density_nonbonded(Model& model);

  void add_expected_number_density_plm(Model& model, const int iseq);
  void copy_expected_number_density(Model& model);
  
  void set_transfer_matrices(Model& model);

  void set_IT_II(const int i, darray2& IT);
  void set_IT_IIgen(const int i, darray2& IT);// generalized inverse version.
  void set_specT_II(void); 
  void set_U(void);

  void set_partial_partition_function(void);
  void set_lnZ_I(void);
  void set_mean_field_zero(void);
  void set_mean_field_obs(Model& model);
  void set_mean_field(Model& model, const int iseq);
  void update_mean_field(Model& model, const double alpha);

  
  void init_looplength(Model& model);
  void init_lambda(void);
  double check_lambda(void); // returns the magnitude of lambda
  double update_lambda(const double alpha);


  void make_scf_insert(Model& model);
  void make_field(Model& model, const bool lambda_p=false);
  double get_free_energy(void);
  
  double update_J(Model& m, const double sca);
  double learn_J(Model& model, const int max_steps, const double step_size);

  double run(Model& model);
  double runSCF(Model& model, const int niter=500, const double eps=1e-10, const double alpha=0.3);
  
  void PLM_sample(Model& model);
  
};

#endif // _PARTFUNC_H_
