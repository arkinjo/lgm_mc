#include "PartFunc.hpp"


static const double ascale_min = 0.001;

void my_gsl_handler (const char * reason, 
		     const char * file, 
		     int line, 
		     int gsl_errno)
{
  cout << "gsl: " << file << ":" << line << ":" << reason << endl;
  return;
}

/* generalized inverse of square matrix */

static gsl_matrix *svdX = gsl_matrix_alloc(nstate_insert,nstate_insert);
static gsl_matrix *svdV = gsl_matrix_alloc(nstate_insert,nstate_insert);
static gsl_vector *svdS = gsl_vector_alloc(nstate_insert);
static gsl_vector *svdwork = gsl_vector_alloc(nstate_insert);

PartFunc::PartFunc(Model& model)
{
  model_length = model.length();
  beta = model.beta;
  
  darray1::extent_gen ext1;
  darray2::extent_gen ext2;
  darray3::extent_gen ext3;
  darray4::extent_gen ext4;
  
  lambda.resize(ext1[model_length+1]);
  looplength.resize(ext1[model_length+1]);
  // partition function etc.
  lnZf.resize(ext2[model_length+2][nstate_core]);
  lnZb.resize(ext2[model_length+2][nstate_core]);
  lnZf_I.resize(ext2[model_length+1][nstate_insert]);
  lnZb_I.resize(ext2[model_length+1][nstate_insert]);
  // transfer matrices etc.
  specT_II.resize(ext2[model_length+1][nstate_insert]);
  T_OO.resize(ext3[model_length+1][nstate_core][nstate_core]);
  T_OI.resize(ext3[model_length+1][nstate_core][nstate_insert]);
  T_IO.resize(ext3[model_length+1][nstate_insert][nstate_core]);
  T_II.resize(ext3[model_length+1][nstate_insert][nstate_insert]);
  IT_II.resize(ext3[model_length+1][nstate_insert][nstate_insert]);
  U.resize(ext3[model_length+1][nstate_core][nstate_core]);
  // "Frozen" fields
  Kmf_O.resize(ext2[model_length+2][nstate_core]);

  n0_O.resize(ext2[model_length+2][nstate_core]);
  en_O.resize(ext2[model_length+2][nstate_core]);
  en_I.resize(ext2[model_length+2][nstate_insert]);

  ens_OO.resize(ext3[model_length+2][nstate_core][nstate_core]);
  ens_OI.resize(ext3[model_length+2][nstate_core][nstate_insert]);
  ens_IO.resize(ext3[model_length+2][nstate_insert][nstate_core]);
  ens_II.resize(ext3[model_length+2][nstate_insert][nstate_insert]);

  enl_OO.resize(ext4[model_length+2][model_length+2][nstate_core][nstate_core]);

  IntPairs &intpairs = model.get_intpairs();
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    const int i = ((it->i) < (it->j)) ? it->i : it->j;
    const int j = ((it->i) < (it->j)) ? it->j : it->i;
    if(IntList.count(i) == 0) {
      vector<int> v;
      IntList[i] = v;
    }
    IntList[i].push_back(j);
  }
}

void PartFunc::set_specT_II(void)
{
  for(int i = 1; i < model_length; ++i) {
    double tmaxi = neg_infinity;
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	if(tmaxi < T_II[i][a][b]) tmaxi = T_II[i][a][b];
      }
    }
    double norm = 0.0;
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	double it = exp(T_II[i][a][b] - tmaxi);
	if(it < 1.0e-10) it = 0;
	norm += it*it;
	gsl_matrix_set(svdX,a,b,it);
      }
    }
    gsl_set_error_handler(&my_gsl_handler);

    int svret = gsl_linalg_SV_decomp(svdX,svdV,svdS,svdwork);
    if(svret) {
      cout << "SV_decomp returned " << svret << endl;
      for(int a = 0; a < nstate_insert; ++a) {
	for(int b = 0; b < nstate_insert; ++b) {
	  cout << "T_II " << i << "\t" << a << "\t" << b << "\t" 
	       << exp(T_II[i][a][b] - tmaxi) << endl;
	}
      }
      abort();
    }

    for(int a = 0; a < nstate_insert; ++a) {
      specT_II[i][a] = log(gsl_vector_get(svdS,a)) + tmaxi;
    }
  }

  return;
}

void PartFunc::set_IT_IIgen(const int i, darray2& IT)
{
  double it;

  for(int a = 0; a < nstate_insert; ++a) {
    for(int b = 0; b < nstate_insert; ++b) {
      it = ((a == b) ? 1.0 : 0) - exp(T_II[i][a][b] - lambda[i]);
      gsl_matrix_set(svdX,a,b,it);
      IT[a][b] = 0.0;
    }
  }

  gsl_linalg_SV_decomp(svdX,svdV,svdS,svdwork);

  double smax = gsl_vector_get(svdS,0);
  for(int c = 0; c < nstate_insert; ++c) {
    double s = gsl_vector_get(svdS,c);
    if (s/smax > 1.0e-8) {
      for(int a = 0; a < nstate_insert; ++a) {
	for(int b = 0; b < nstate_insert; ++b) {
	  IT[a][b] += gsl_matrix_get(svdV,a,c) * gsl_matrix_get(svdX,b,c) / s;
	}
      }
    }
  }

  int negv = 0;
  for(int a = 0; a < nstate_insert; ++a) {
    for(int b = 0; b < nstate_insert; ++b) {
      if(IT[a][b] < 0.0 && fabs(IT[a][b]) > 1.0e-8) negv++;
      double it = (IT[a][b] > 0.0) ? log(IT[a][b]) : neg_infinity;
      IT[a][b] = IT_II[i][a][b] = it;
    }
  }

  return;
}


static ublas::matrix<double> lu_input(nstate_insert,nstate_insert);
static ublas::matrix<double> lu_inv(nstate_insert,nstate_insert);

void PartFunc::set_IT_II(const int i, darray2& IT)
{
  double it;

  for(int a = 0; a < nstate_insert; ++a) {
    for(int b = 0; b < nstate_insert; ++b) {
      lu_input(a,b) = ((a == b) ? 1.0 : 0) - exp(T_II[i][a][b] - lambda[i]);
    }
  }

  InvertMatrix(lu_input,lu_inv);

  int negv = 0;
  for(int a = 0; a < nstate_insert; ++a) {
    for(int b = 0; b < nstate_insert; ++b) {
      if(lu_inv(a,b) < 0.0 && fabs(lu_inv(a,b)) > 1.0e-8) negv++;
      double it = (lu_inv(a,b) > 0.0) ? log(lu_inv(a,b)) : neg_infinity;
      IT[a][b] = IT_II[i][a][b] = it;
    }
  }
  return;
}


void PartFunc::set_U(void)
{
  darray2 IT(boost::extents[nstate_insert][nstate_insert]);
  darray2 tIO(boost::extents[nstate_insert][nstate_core]); // = T_IO
  darray2 tOI(boost::extents[nstate_core][nstate_insert]); // = T_OI
  darray2 tOO(boost::extents[nstate_core][nstate_core]); // = T_OI*IT*T_IO
  darray2 tIIO(boost::extents[nstate_insert][nstate_core]); // = IT*T_IO

  int a = terminal_state;
  for(int b = 0; b < nstate_core; ++b) {
    U[0][terminal_state][b] = T_OO[0][terminal_state][b];
  }

  for(int a = 0; a < nstate_core; ++a) {
    U[model_length][a][terminal_state] = T_OO[model_length][a][terminal_state];
  }

  for(int i = 1; i < model_length; ++i) {
    if(lambda[i] > specT_II[i][0]) {
      set_IT_II(i,IT);
    } else {
      set_IT_IIgen(i,IT);
    }

    std::fill(tIO.data(),tIO.data() + tIO.num_elements(), neg_infinity);
    std::fill(tOI.data(),tOI.data() + tOI.num_elements(), neg_infinity);
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	tIO[a][b] = T_IO[i][a][b];
      }
    }
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	tOI[a][b] = T_OI[i][a][b];
      }
    }

    log_mat_mult(IT, tIO, tIIO, nstate_insert, nstate_insert, nstate_core);
    log_mat_mult(tOI, tIIO, tOO, nstate_core, nstate_insert, nstate_core);

    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	U[i][a][b] = log_add(T_OO[i][a][b], tOO[a][b]);
      }
    }
  }
  
  return;
}

void PartFunc::set_transfer_matrices(Model& m)
{
  for(int i = 0; i <= model_length; ++i) {
    if(i == 0) {
      int a = terminal_state;
      for(int b = 0; b < nstate_core; ++b) {
	T_OO[i][a][b] = beta*(m.H_O[i+1][b] + Kmf_O[i+1][b]);
      }
    }
    else if(i == model_length) {
      int b = terminal_state;
      for(int a = 0; a < nstate_core; ++a) {
	T_OO[i][a][b] = 0.0;
      }
    }
    else {
      for(int a = 0; a < nstate_insert; ++a) {
	for(int b = 0; b < nstate_insert; ++b) {
	  T_II[i][a][b] =  beta*(m.J_II[i][a][b] + m.H_I[i][b]);
	}
      }
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  T_OO[i][a][b] = beta*(m.J_OO[i][a][b] + m.H_O[i+1][b] + Kmf_O[i+1][b]);
	}
	for(int b = 0; b < nstate_insert; ++b) {
	  T_OI[i][a][b] = beta*(m.J_OI[i][a][b] + m.H_I[i][b]);
	}
      }
      for(int a = 0; a < nstate_insert; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  T_IO[i][a][b] = beta*(m.J_IO[i][a][b] + m.H_O[i+1][b] + Kmf_O[i+1][b]);
	}
      }
    }
  }

  return;
}

void PartFunc::set_partial_partition_function(void)
{
  darray2 ut(boost::extents[nstate_core][nstate_core]);

  darray2 zp(boost::extents[1][nstate_core]);
  darray2 zn(boost::extents[1][nstate_core]);

  std::fill(lnZf.data(),lnZf.data() + lnZf.num_elements(), neg_infinity);

  // Zf{i} = Z_{i-1}U_{i-1,i}
  lnZf[0][terminal_state] = 0.0;
  for(int i = 1; i <= model_length+1; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      zp[0][a] = lnZf[i-1][a];
      for(int b = 0; b < nstate_core; ++b) {
	ut[a][b] = U[i-1][a][b];
      }
    }
    log_mat_mult(zp, ut, zn, 1, nstate_core, nstate_core);
    for(int a = 0; a < nstate_core; ++a) {
      lnZf[i][a] = zn[0][a];
    }
  }

  darray2 wp(boost::extents[nstate_core][1]);
  darray2 wn(boost::extents[nstate_core][1]);

  std::fill(lnZb.data(),lnZb.data() + lnZb.num_elements(), neg_infinity);
  lnZb[model_length+1][terminal_state] = 0.0;

  for(int i = model_length; i >= 0; --i) {
    for(int a = 0; a < nstate_core; ++a) {
      wp[a][0] = lnZb[i+1][a];
      for(int b = 0; b < nstate_core; ++b) {
	ut[a][b] = U[i][a][b];
      }
    }
    log_mat_mult(ut, wp, wn, nstate_core, nstate_core, 1);
    for(int a = 0; a < nstate_core; ++a) {
      lnZb[i][a] = wn[a][0];
    }
  }

  set_lnZ_I();

  return;
}

void PartFunc::set_lnZ_I(void)
{
  darray2 tIO(boost::extents[nstate_insert][nstate_core]);
  darray2 tOI(boost::extents[nstate_core][nstate_insert]);
  darray2 tIT(boost::extents[nstate_insert][nstate_insert]);
  darray2 tOII(boost::extents[nstate_core][nstate_insert]);
  darray2 tIIO(boost::extents[nstate_insert][nstate_core]);
  darray1 vO(boost::extents[nstate_core]);
  darray1 vI(boost::extents[nstate_core]);

  for(int i = 0; i <= model_length; ++i) 
    {
      std::fill(tOI.data(), tOI.data() + tOI.num_elements(),neg_infinity);
      std::fill(tIO.data(), tIO.data() + tIO.num_elements(),neg_infinity);
      std::fill(tIIO.data(), tIIO.data() + tIIO.num_elements(),neg_infinity);
      std::fill(tOII.data(), tOII.data() + tOII.num_elements(),neg_infinity);

      // T_{O_iI_i}(I - T_{I_iI_i})^{-1} -> tOII
      for(int a = 0; a < nstate_core; ++a) 
	{
	  for(int b = 0; b < nstate_insert; ++b) 
	    {
	      tOI[a][b] = T_OI[i][a][b];
	    }
	}
      for(int a = 0; a < nstate_insert; ++a) 
	{
	  for(int b = 0; b < nstate_insert; ++b) 
	    {
	      tIT[a][b] = IT_II[i][a][b];
	    }
	}
      log_mat_mult(tOI,tIT, tOII, nstate_core, nstate_insert,nstate_insert);

      // Z_{O_i} [T_{O_iI_i}(I - T_{I_iI_i})^{-1}] in logarithm.
      for(int a = 0; a < nstate_core; ++a) 
	{
	  if(i > 0) {
	    vO[a] = lnZf[i][a];
	  }
	  else {
	    vO[a] = (a == terminal_state) ? 0 : neg_infinity;
	  }
	}
      for(int b = 0; b < nstate_insert; ++b) 
	{
	  for(int a = 0; a < nstate_core; ++a) 
	    {
	      vI[a] = tOII[a][b];
	    }
	  lnZf_I[i][b] = log_dot(vO,vI,nstate_core);
	}

      // (I - T_{I_iI_i})^{-1}T_{I_iO_{i+1}}
      for(int a = 0; a < nstate_insert; ++a)
	{
	  for(int b = 0; b < nstate_insert; ++b) 
	    {
	      tIT[a][b] = IT_II[i][a][b];
	    }
	}
      for(int a = 0; a < nstate_insert; ++a) 
	{
	  for(int b = 0; b < nstate_core; ++b) 
	    {
	      tIO[a][b] = T_IO[i][a][b];
	    }
	}
      log_mat_mult(tIT,tIO, tIIO, nstate_insert, nstate_insert,nstate_core);
      // (I - T_{I_iI_i})^{-1}T_{I_iO_{i+1}} W_{O_{i+1}}
      for(int b = 0; b < nstate_core; ++b) 
	{
	  if(i < model_length) {
	    vO[b] = lnZb[i+1][b];
	  }
	  else {
	    vO[b] = (b == terminal_state) ? 0.0 : neg_infinity;
	  }
	}
      for(int a = 0; a < nstate_insert; ++a) 
	{
	  for(int b = 0; b < nstate_core; ++b) 
	    {
	      vI[b] = tIIO[a][b];
	    }
	  lnZb_I[i][a] = log_dot(vO,vI,nstate_core);
	}
    }

  return;
}

double PartFunc::get_free_energy(void)
{
  return -lnZf[model_length+1][terminal_state]/beta;
}

void PartFunc::init_expected_number_density(Model& m)
{
  std::fill(en_O.data(), en_O.data() + en_O.num_elements(), 0.0);
  en_O[0][terminal_state] = en_O[model_length+1][terminal_state] = 1.0;
  std::fill(en_I.data(), en_I.data() + en_I.num_elements(), 0.0);

  std::fill(ens_OO.data(), ens_OO.data() + ens_OO.num_elements(), 0.0);
  std::fill(ens_OI.data(), ens_OI.data() + ens_OI.num_elements(), 0.0);
  std::fill(ens_IO.data(), ens_IO.data() + ens_IO.num_elements(), 0.0);
  std::fill(ens_II.data(), ens_II.data() + ens_II.num_elements(), 0.0);

  std::fill(enl_OO.data(), enl_OO.data() + enl_OO.num_elements(), 0.0);
  return;
}

void PartFunc::set_mean_field_zero(void)
{
  std::fill(Kmf_O.data(), Kmf_O.data() + Kmf_O.num_elements(), 0.0);
  return;
}

void PartFunc::set_mean_field_obs(Model& model)
{
  set_mean_field_zero();
  for(int i = 1; i <= model_length; ++i) {
    if(IntList.count(i) == 0) continue;
    int n = IntList[i].size();
    for(int k = 0; k < n; ++k) {
      int j = IntList[i][k];
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  Kmf_O[i][a] += model.K_OO[i][j][a][b]*model.n_O[j][b];
	  Kmf_O[j][b] += model.K_OO[i][j][a][b]*model.n_O[i][a];
	}
      }
    }
  }
}

void PartFunc::set_mean_field(Model& model, const int iseq)
{
  modelSeq &aseq = model.get_PreModel(iseq);
  const double pc = model.lambdaH_O;
  double d;

  set_mean_field_zero();
  for(int i = 1; i <= model_length; ++i) {
    if(IntList.count(i) == 0) continue;
    int n = IntList[i].size();
    for(int k = 0; k < n; ++k) {
      int j = IntList[i][k];
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  Kmf_O[i][a] += model.K_OO[i][j][a][b]*en_O[j][b];
	  Kmf_O[j][b] += model.K_OO[i][j][a][b]*en_O[i][a];
	}
      }
    }
  }

  return;
}

void PartFunc::update_mean_field(Model& model, const double alpha)
{
  for(int i = 1; i <= model_length; ++i) {
    if(IntList.count(i) == 0) continue;
    int n = IntList[i].size();
    for(int a = 0; a < nstate_core; ++a) {
      double mf=0.0;
      for(int k = 0; k < n; ++k) {
	int j = IntList[i][k];
  	for(int b = 0; b < nstate_core; ++b) {
	  mf += model.K_OO[i][j][a][b]*en_O[j][b];
	}
      }
      double delta = mf - Kmf_O[i][a];
      Kmf_O[i][a] += alpha*delta;
    }
  }

  return;
}

void PartFunc::init_looplength(Model& model)
{
  for(int i = 1; i < model_length; ++i) {
    looplength[i] = 0;
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	looplength[i] += model.ns_II[i][a][b];
      }
    }
  }
  return;
}

double PartFunc::check_lambda(void)
{
  double dlam = 0.0;

  for(int i = 1; i < model_length; ++i) {
    if(dlam < lambda[i]) dlam = lambda[i];
  }

  return dlam;
}

void PartFunc::init_lambda(void)
{
  // initialize lambda.
  const double eps = 1;

  for(int i = 1; i < model_length; ++i) {
    double norm = 0.0;
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	double it = exp(T_II[i][a][b]);
	norm += it*it;
      }
    }
    norm = sqrt(norm);
    lambda[i] = specT_II[i][0] + log(2);
  }

  return;
}

double PartFunc::update_lambda(const double alpha)
{
  double del = 0.0;
  double en;
  double lnZ = -get_free_energy()*beta;
  darray2 It(boost::extents[1][nstate_insert]), 
    fi(boost::extents[1][nstate_insert]),
    bi(boost::extents[nstate_insert][1]), 
    ItI(boost::extents[1][1]),
    T(boost::extents[nstate_insert][nstate_insert]);

  for(int i = 1; i < model_length; ++i) {
    double len = 0.0;
    for(int a = 0; a < nstate_insert; ++a) {
      fi[0][a] = lnZf_I[i][a];
      bi[a][0] = lnZb_I[i][a];
      for(int b = 0; b < nstate_insert; ++b) {
	T[a][b] = T_II[i][a][b];
       }
    }

    log_mat_mult(fi,T,It, 1, nstate_insert, nstate_insert);
    log_mat_mult(It,bi, ItI, 1, nstate_insert, 1);

    double nlambda = ItI[0][0] - lnZ - log(looplength[i]);

    if(std::isnan(nlambda) || std::isinf(nlambda)) {
      double looplen = exp(ItI[0][0] - lnZ - lambda[i]);
      cerr << "lambda is nan/inf! " << i << "\t"
	   << lambda[i] << "\t"
	   << looplen << "\t"
	   <<  ItI[0][0] << "\t" 
	   << lnZ << "\t" << log(looplength[i]) << endl;
      double s = looplength[i];
      double fs = 2*s/(2*s + 1 - sqrt(4*s + 1));
      nlambda = specT_II[i][0] + log(fs);
    }

    double olambda = lambda[i];
    double d = olambda - nlambda;
    double tlambda = olambda - alpha*d;
    const double eps = 1;

    lambda[i] = tlambda;
    del += d*d;
  }
  
  return sqrt(del/(model_length+1));
}

void PartFunc::make_scf_insert(Model& model) 
{
  const int niter = 300;
  double dlam0 = 1.0e+10;
  const double eps = 1e-14;
  double alpha = 0.1;

  set_specT_II();
  init_lambda();

  int iter;
  for(iter=1; iter <= niter; ++iter) {
    set_U();
    set_partial_partition_function();

    double dlam = update_lambda(alpha);
    if (dlam > dlam0) {
      alpha *= 0.9;
    }
    dlam0 = dlam;
    if(dlam < eps) break;
  }

  if(dlam0 > eps) {
    double mag = check_lambda();
    cerr << "scf_insert: " << iter << " dlam= " << dlam0
	 << "\t||lambda||= " << mag << endl;
  }

  return;
}

void PartFunc::make_field(Model& model, const bool lambda_p)
{

  set_transfer_matrices(model);
  if(lambda_p) {
      std::fill(lambda.data(),lambda.data() + lambda.num_elements(), 0);
      make_scf_insert(model);
  }
  set_U();
  set_partial_partition_function();
}

void PartFunc::copy_expected_number_density(Model& m)
{
  std::copy(en_O.data(), en_O.data() + en_O.num_elements(), m.en_O.data());
  std::copy(en_I.data(), en_I.data() + en_I.num_elements(), m.en_I.data());

  std::copy(ens_OO.data(), ens_OO.data() + ens_OO.num_elements(), m.ens_OO.data());
  std::copy(ens_OI.data(), ens_OI.data() + ens_OI.num_elements(), m.ens_OI.data());
  std::copy(ens_IO.data(), ens_IO.data() + ens_IO.num_elements(), m.ens_IO.data());
  std::copy(ens_II.data(), ens_II.data() + ens_II.num_elements(), m.ens_II.data());

  set_expected_number_density_nonbonded(m);
}

void PartFunc::set_expected_number_density(void)
{

  const double lnZ_total = -get_free_energy()*beta;
  
  for(int i = 0; i <= model_length; ++i) {
    // single site
    if(i > 0) {
      for(int a = 0; a < nstate_core; ++a) {
	en_O[i][a] = exp(lnZf[i][a] + lnZb[i][a] - lnZ_total);
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      en_I[i][a] = exp(lnZf_I[i][a] + lnZb_I[i][a] - lnZ_total);
    }

    // short-range pairs
    for(int a = 0; a < nstate_core; ++a) {
      if(i == 0 && a != terminal_state) continue;
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	ens_OO[i][a][b] = exp(lnZf[i][a] + T_OO[i][a][b] + lnZb[i+1][b] - lnZ_total);
      }
      for(int b = 0; b < nstate_insert; ++b) {
	ens_OI[i][a][b] = exp(lnZf[i][a] + T_OI[i][a][b] + lnZb_I[i][b] - lnZ_total);
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	ens_IO[i][a][b] = exp(lnZf_I[i][a] + T_IO[i][a][b] + lnZb[i+1][b] - lnZ_total);
      }
      for(int b = 0; b < nstate_insert; ++b) {
	ens_II[i][a][b] = exp(lnZf_I[i][a]  + T_II[i][a][b] + lnZb_I[i][b] - lnZ_total - lambda[i]);
      }
    }
  }

  return;
}

void PartFunc::set_expected_number_density_nonbonded(Model& m)
{
  for(int i = 0; i <= model_length; ++i) {
    if(i > 0) {
      for(int a = 0; a < nstate_core; ++a) {
	m.enl_OO[i][i][a][a] = en_O[i][a];
	for(int j = i; j <= model_length; ++j) {
	  if(j > i) {
	    for(int b = 0; b < nstate_core; ++b) {
	      m.enl_OO[j][i][b][a] = m.enl_OO[i][j][a][b] = enl_OO[i][j][a][b];
	    }
	  }
	  for(int b = 0; b < nstate_insert; ++b) {
	    m.enl_IO[j][i][b][a] = m.enl_OI[i][j][a][b] = en_O[i][a]*en_I[j][b];
	  }
	}
      }
    }
    for(int j = i; j <= model_length; ++j) {
      for(int a = 0; a < nstate_insert; ++a) {
	for(int b = 0; b < nstate_insert; ++b) {
	  m.enl_II[j][i][b][a] = m.enl_II[i][j][a][b] = en_I[i][a]*en_I[j][b];
	}
      }
    }
  }
  return;
}

void PartFunc::set_expected_number_density_nocorr()
{
  for(int i = 1; i <= model_length-2; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      for(int j = i+2; j <= model_length; ++j) {
	for(int b = 0; b < nstate_core; ++b) {
	  enl_OO[i][j][a][b] = en_O[i][a]*en_O[j][b];
	  enl_OO[j][i][b][a] = enl_OO[i][j][a][b];
	}
      }
    }
  }
  return;
}

void PartFunc::add_expected_number_density_plm(Model& m, const int iseq)
{
  const double w = m.get_weight(iseq);
  const double w2 = 0.5*w;
  
  for(int i = 0; i <= model_length; ++i) {
    if(i > 0) {
      for(int a = 0; a < nstate_core; ++a) {
	m.en_O[i][a] += w*en_O[i][a];
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      m.en_I[i][a] += w*en_I[i][a];
    }

    for(int a = 0; a < nstate_core; ++a) {
      if(i == 0 && a != terminal_state) continue;
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	m.ens_OO[i][a][b] += w*ens_OO[i][a][b];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	m.ens_OI[i][a][b] += w*ens_OI[i][a][b];
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	m.ens_IO[i][a][b] += w*ens_IO[i][a][b];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	m.ens_II[i][a][b] += w*ens_II[i][a][b];
      }
    }
  }

  for(int i = 1; i <= model_length; ++i) {
    if(IntList.count(i) == 0) continue;
    int n = IntList[i].size();
    for(int k = 0; k < n; ++k) {
      int j = IntList[i][k];
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  double c = w2*(en_O[i][a]*n0_O[j][b] + n0_O[i][a]*en_O[j][b]);
	  m.enl_OO[i][j][a][b] += c;
	  m.enl_OO[j][i][b][a] += c;
	}
      }
    }
  }
  return;
}

double PartFunc::update_J(Model& m, const double sca)
{
  const double lnZ = -get_free_energy()*beta;
  const double alpha = sca;
  double nval;
  double delta = 0.0,adel,del;
  darray2 joo(boost::extents[nstate_core][nstate_core]);
  darray2 joi(boost::extents[nstate_core][nstate_insert]);
  darray2 jio(boost::extents[nstate_insert][nstate_core]);
  darray2 jii(boost::extents[nstate_insert][nstate_insert]);
  const int c = nstate_match; // "deletion"
  double gauge = 0.0;
      
  // Assuming all H_O & H_I are ZERO during learning!!
  for(int i = 1; i < model_length; i++) {
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	nval = log(m.ns_II[i][a][b]) + lnZ 
	  - lnZf_I[i][a] - lnZb_I[i][b] - m.H_I[i][b];
	del = alpha*(nval - m.J_II[i][a][b]);
	adel=fabs(del);
	if(delta<adel) delta=adel;
	m.J_II[i][a][b] += del;
      }
    }

    gauge = 0.0;
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	joo[a][b] = log(m.ns_OO[i][a][b]) + lnZ - lnZf[i][a] - lnZb[i+1][b] 
	  - Kmf_O[i+1][b] - m.H_O[i+1][b];
      }
    }
    
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	joi[a][b] = log(m.ns_OI[i][a][b]) + lnZ - lnZf[i][a] - lnZb_I[i][b] - m.H_I[i][b];
      }
    }
    
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	jio[a][b] = log(m.ns_IO[i][a][b]) + lnZ - lnZf_I[i][a] - lnZb[i+1][b] 
	  -Kmf_O[i+1][b] - m.H_O[i+1][b];
      }
    }

    gauge = -1.0 * joo[c][c];
    
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	del= alpha*(joo[a][b] - m.J_OO[i][a][b] + gauge);
	adel=fabs(del);
	if(delta<adel) delta=adel;
	m.J_OO[i][a][b] += del;
      }
    }

    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	del = alpha*(joi[a][b] - m.J_OI[i][a][b] + 0.5*gauge);
	adel=fabs(del);
	if(delta<adel) delta=adel;
	m.J_OI[i][a][b] += del;
      }
    }

    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	del = alpha*(jio[a][b] - m.J_IO[i][a][b] + 0.5*gauge);
	adel=fabs(del);
	if(delta<adel) delta=adel;
	m.J_IO[i][a][b] += del;
      }
    }
  }

  return delta;
}

double PartFunc::run(Model& model)
{
  model.initialize_expected_densities();
  
  init_expected_number_density(model);
  init_looplength(model);
  set_mean_field_obs(model);
  std::fill(lambda.data(),lambda.data() + lambda.num_elements(), 0);
  make_field(model);
  set_expected_number_density();
  set_expected_number_density_nocorr();
  copy_expected_number_density(model);

  return get_free_energy();
}

double PartFunc::runSCF(Model& model, const int niter, const double eps, const double alpha)
{
  model.initialize_expected_densities();
  std::fill(lambda.data(),lambda.data() + lambda.num_elements(), 0);
  init_expected_number_density(model);
  init_looplength(model);
  set_mean_field_obs(model);
  double fe = 1e+10;
  for(int iter = 1; iter <= niter; ++iter) {
    make_field(model);
    set_expected_number_density();
    update_mean_field(model, alpha);
    const double fe_new = get_free_energy();
    const double del = fabs(fe_new - fe);
    cerr << "PartFunc::runSCF: FreeEnergy= " << iter << "\t"
	 << fe_new << "\t" << del << endl;
    if(del < eps) break;
    fe = fe_new;
  }
  
  copy_expected_number_density(model);
  return get_free_energy();
}

double PartFunc::learn_J(Model& model, const int max_steps, const double step_size)
{
  double diff, dlam, delta;
  const double epsJ = 1e-12;
  const double epsI = 1e-8;
  double jscale = step_size;

  model.initialize_expected_densities();
  
  init_expected_number_density(model);
  init_looplength(model);
  set_mean_field_obs(model);  

  for(int istep = 1; istep <= max_steps; ++istep) {
    make_field(model, true);
    double fe = get_free_energy();
    dlam = check_lambda();
    delta = update_J(model, jscale);
    cout << istep << " " << "delJ1= "
	 << delta << "\t"
	 << dlam << "\t"
      	 << jscale << "\t"
	 << fe << endl;
    if(dlam < epsI) break;
  }

  double delta0 = 10000;
  if(dlam < epsI) {
    cout << "# Leaning with lambda = 0." << endl;
    for(int istep = 1; istep <= max_steps; ++istep) {
      make_field(model); // make_scf_insert is not called!
      double fe = get_free_energy();
      double delta = update_J(model, jscale);
      cout << istep << " " << "delJ2= "
	   << delta << "\t"
	   << jscale << "\t"
	   << fe << endl;
      if(delta < epsJ && istep > 5) break;

      if(delta>delta0) jscale *= 0.9;
      delta0 = delta;
    }
  }

  set_expected_number_density();
  set_expected_number_density_nocorr();
  copy_expected_number_density(model);
  return diff;
}

void PartFunc::PLM_sample(Model& model)
{
  const int Nseqs = model.get_num_seq();
  double TotalWeight = 0.0;
  model.initialize_expected_densities();
  for(int im = 0; im < Nseqs; ++im) {
    TotalWeight += model.get_weight(im);
    if((im+1) % 1000 == 0) {
      cerr << "PLM_sample: " << im + 1 << " sequences processed..." << endl;
    }
    set_mean_field(model, im);
    set_transfer_matrices(model);
    set_U();
    set_partial_partition_function();
    set_expected_number_density();
    add_expected_number_density_plm(model, im);
  }
  model.normalize_expected_densities(TotalWeight);
  return;
}
