#include "OptParam.hpp"

OptParam::OptParam(const Model& model)
{
  Length = model.length();
  set_default_alpha(model);
  use_WL = false;
  beta_H = 0.0;
  beta_J = 0.95;
  beta_K = 0.95;

  eps_H = 1.0e-3;
  eps_J = 1.0e-3;
  eps_K = 1.0e-3;
  
  nsweeps = 100;
  nsweeps_eq = 1000;

  darray3::extent_gen ext2;
  H_O.resize(ext2[Length+2][nstate_core]);
  std::fill(H_O.data(), H_O.data() + H_O.num_elements(), 0.0);
  
  H_I.resize(ext2[Length+2][nstate_insert]);
  std::fill(H_I.data(), H_I.data() + H_I.num_elements(), 0.0);
  
  darray3::extent_gen ext3;
  J_OO.resize(ext3[Length+2][nstate_core][nstate_core]);
  std::fill(J_OO.data(), J_OO.data() + J_OO.num_elements(), 0.0);

  J_OI.resize(ext3[Length+2][nstate_core][nstate_insert]);
  std::fill(J_OI.data(), J_OI.data() + J_OI.num_elements(), 0.0);

  J_IO.resize(ext3[Length+2][nstate_insert][nstate_core]);
  std::fill(J_IO.data(), J_IO.data() + J_IO.num_elements(), 0.0);

  J_II.resize(ext3[Length+2][nstate_insert][nstate_insert]);
  std::fill(J_II.data(), J_II.data() + J_II.num_elements(), 0.0);
  
  darray4::extent_gen ext4;
  K_OO.resize(ext4[Length+2][Length+2][nstate_core][nstate_core]);
  std::fill(K_OO.data(), K_OO.data() + K_OO.num_elements(), 0.0);

}

void OptParam::set_default_alpha(const Model& model)
{
  alpha_H = 0.0;
  alpha_J = 0.01;
  alpha_K = 0.01;

  return;
}

double OptParam::update_H(Model& model)
{
  double delta = 0.0;
  double del,diff;
  double cnt=0.0;

  for(int i = 1; i <= Length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      diff = model.n_O[i][a] - model.en_O[i][a]
	; //- model.lambdaH_O * model.H_O[i][a];
      del = alpha_H*diff + beta_H*H_O[i][a];
      model.H_O[i][a] += del;
      H_O[i][a] = del;
      delta += diff*diff;
      cnt++;
    }

    if(i == Length) continue;
    for(int a = 0; a < nstate_insert; ++a) {
      diff = model.n_I[i][a] - model.en_I[i][a]
	; //- model.lambdaH_I * model.H_I[i][a];
      del = alpha_H*diff + beta_H*H_I[i][a];
      model.H_I[i][a] += del;
      H_I[i][a] = del;
      delta += diff*diff;
      cnt++;
    }
  }

  return sqrt(delta/cnt);
}

double OptParam::update_J(Model& model)
{
  double del,diff;
  double cnt = 0.0;
  double delta = 0.0;
  const int c = terminal_state;
  for(int i = 1; i < Length; ++i) {
    const double lenI = model.len_I[i];
    const double gauge = -(model.ns_OO[i][c][c] - model.ens_OO[i][c][c]);
    const double gauge2 = 0.5*gauge;
    for(int a = 0; a < nstate_core; ++a) {
      const double na = model.n_O[i][a];
      for(int b = 0; b < nstate_core; ++b) {
	if(a == terminal_state && b == terminal_state) continue;
	const double corr = fabs(model.ns_OO[i][a][b] - na*model.n_O[i+1][b]
				 - model.cs0_OO);
	diff = gauge;
	if(corr > model.NS_CUTOFF) 
	  diff += model.ns_OO[i][a][b] - model.ens_OO[i][a][b];
	del = alpha_J*diff + beta_J*J_OO[i][a][b];
	model.J_OO[i][a][b] += del;
	J_OO[i][a][b] = del;
	delta += diff*diff;
	cnt++;
      }
      for(int b = 0; b < nstate_insert; ++b) {
	const double corr = fabs(model.ns_OI[i][a][b] - na*model.n_Int[i][b]
				 - model.cs0_OI);
	diff = gauge2;
	if(corr > model.NS_CUTOFF) 
	  diff += model.ns_OI[i][a][b] - model.ens_OI[i][a][b];
	del = alpha_J*diff + beta_J*J_OI[i][a][b];
	model.J_OI[i][a][b] += del;
	J_OI[i][a][b] = del;
	delta += diff*diff;
	cnt++;
      }
    }

    for(int a = 0; a < nstate_insert; ++a) {
      double na = model.n_Ict[i][a];
      for(int b = 0; b < nstate_core; ++b) {
	const double corr = fabs(model.ns_IO[i][a][b] - na*model.n_O[i+1][b]
				 - model.cs0_IO);
	diff = gauge2;
	if(corr > model.NS_CUTOFF) 
	  diff += model.ns_IO[i][a][b] - model.ens_IO[i][a][b];
	del = alpha_J*diff + beta_J*J_IO[i][a][b];
	model.J_IO[i][a][b] += del;
	J_IO[i][a][b] = del;
	delta += diff*diff;
	cnt++;
      }
      na = model.n_I[i][a];
      for(int b = 0; b < nstate_insert; ++b) {
      	const double corr = fabs(model.ns_II[i][a][b] - na*model.n_I[i][b]
      				 - model.cs0_II);
	if(corr <= model.NS_CUTOFF) continue;
      	diff = model.ns_II[i][a][b] - model.ens_II[i][a][b];
      	del = alpha_J*diff + beta_J*J_II[i][a][b];
      	if(model.J_II[i][a][b] + del >= 0) {
	  model.J_II[i][a][b] = -1e-5;
	} else {
	  model.J_II[i][a][b] += del;
	}
      	J_II[i][a][b] = del;
      	delta += diff*diff;
      	cnt++;
      }
    }
  }

  return sqrt(delta/cnt);
}

double OptParam::update_K(Model& model)
{
  double delta = 0.0;
  double del, diff;
  IntPairs& intpairs = model.get_intpairs();
  int i,j;
  const int c = terminal_state;
  
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    i = it->i; j = it->j;
    const double gauge = model.nl_OO[i][j][c][c] - model.enl_OO[i][j][c][c];
    for (int a = 0; a < nstate_core; ++a) {
      const double na = model.n_O[i][a];
      for (int b = 0; b < nstate_core; ++b) {
	if(a == terminal_state && b == terminal_state) continue;
	const double nOO = model.nl_OO[i][j][a][b];
	const double corr = fabs(nOO - (na * model.n_O[j][b]) - model.cl0_OO);
	diff = -gauge;
	if(corr > model.NL_CUTOFF) 
	  diff += (nOO - model.enl_OO[i][j][a][b]);
	del = alpha_K*diff + beta_K*K_OO[i][j][a][b];
	model.K_OO[i][j][a][b] += del;
	model.K_OO[j][i][b][a] = model.K_OO[i][j][a][b];
	K_OO[i][j][a][b] = del;
	delta += diff*diff;
      }
    }
  }

  return sqrt(delta/intpairs.size());
}

void OptParam::update_all(Model& model, const int iter, const int nsamples)
{
  cout.setf(std::ios::scientific);

  double delH = update_H(model);
  double delJ = update_J(model);
  double delK = update_K(model);
  //  model.set_Ising_gauge();
  const double eobs = model.get_energy();
  const double eexp = model.get_energy_exp();
  const double slen = model.get_seqlength_exp();
  cout << "OptParam::update: " << iter
       << "\t" << nsamples
       << "\t" << delH
       << "\t" << delJ
       << "\t" << delK
       << "\t" << eobs
       << "\t" << eexp - eobs
       << "\t" << slen
       << endl;
}

double OptParam::update_iala(Model& model, const int i)
{
  double delta = 0.0;
  double del,diff;
  double cnt=0.0;
  const int ALA = find_amino1_core('A');
  
  for(int a = 0; a < nstate_core; ++a) {
    const double tn = (a == ALA) ? 0.95 : (0.05/20);
    diff = tn - model.en_O[i][a];
    del = alpha_H*diff + beta_H*H_O[i][a];
    model.H_O[i][a] += del;
    H_O[i][a] = del;
    delta += diff*diff;
    cnt++;
  }
  return sqrt(delta/cnt);
}

void OptParam::optimize(Model& model, MCseq& mcs, const int niter)
{
  int iter;

  for(iter = 1; iter <= niter; ++iter) {
    mcs.get_energy(model);
    mcs.run(model, nsweeps_eq, false);
    mcs.init_stats(model);
    int nsamples = mcs.run(model, nsweeps, true);
    mcs.finish_stats(model);
    assert(nsamples == nsweeps);
    update_all(model, iter, nsamples);
  }
}

void OptParam::optimize_mul(Model& model, MCseq& mcs, const int niter)
{
  int iter,nsamples;
  bool gibbs = mcs.gibbs_sampler_p;

  for(iter = 1; iter <= niter; ++iter) {
    mcs.init_stats(model);
    for(int irun=1; irun <= 10; ++irun) {
      if(irun%2) 
	mcs.set_natural_sequence(model);
      else
	mcs.set_random_sequence(model,10);
      mcs.gibbs_sampler_p = true;
      mcs.run(model, nsweeps_eq, false);
      mcs.gibbs_sampler_p = gibbs;
      nsamples = mcs.run(model, nsweeps, true);
    }
    mcs.finish_stats(model);
    update_all(model, iter, nsamples);
  }
}

void OptParam::optimize_ala(Model& model, MCseq& mcs, const int iala,
			    const int niter)
{
  int iter,nsamples;
  bool gibbs = mcs.gibbs_sampler_p;

  for(iter = 1; iter <= niter; ++iter) {
    mcs.init_stats(model);
    for(int irun=1; irun <= 10; ++irun) {
      if(irun%2) 
	mcs.set_natural_sequence(model);
      else
	mcs.set_random_sequence(model,10);
      mcs.gibbs_sampler_p = true;
      mcs.run(model, nsweeps_eq, false);
      mcs.gibbs_sampler_p = gibbs;
      nsamples = mcs.run(model, nsweeps, true);
    }
    mcs.finish_stats(model);
    const double del = update_iala(model, iala);
    cout << "OptParam::update: " << iter
	 << "\t" << nsamples
	 << "\t" << del << endl;

  }
}

void OptParam::WLoptimize(Model& model, MCseq& mcs, const int niter)
{
  int iter;

  mcs.WLprec = 1e-5;
  for(iter = 1; iter <= niter; ++iter) {
    mcs.scan_natural_sequences(model);
    double DENE = sqrt(mcs.evar_natural);
    double min_energy = ceil(mcs.eave_natural - 3*DENE);
    double max_energy = floor(mcs.eave_natural + 3*DENE);
    // double tmax_energy = ceil(mcs.emax_natural + 10.0);
    // if(max_energy < tmax_energy) max_energy = tmax_energy;

    while(1) {
      mcs.set_natural_sequence(model);
      //mcs.set_random_sequence(model);
      if(mcs.Energy > min_energy && mcs.Energy < max_energy) break;
    }
    cerr << "nat.seq. " << "\t"
	 << mcs.Energy << "\t"
	 << mcs.imin_natural << "\t"
	 << mcs.emin_natural << "\t"
      	 << mcs.imax_natural << "\t"
      	 << mcs.emax_natural << "\t"
	 << mcs.eave_natural << "\t"
      	 << sqrt(mcs.evar_natural) << endl;
    
    mcs.select_type = MULTI_CANONICAL;
    mcs.setup_WL(min_energy, max_energy, 0.1);

    int nwl=0;
    while(1) {
      cerr << "#Determining WL DOS... " << ++nwl << endl;
      const int isweeps = mcs.run(model, nsweeps, false);
      mcs.save_DOS("hoge.dos");
      if(isweeps < nsweeps) break;
    }

    while(1) {
      mcs.set_natural_sequence(model);
      if(mcs.Energy > mcs.MIN_ENERGY
	 && mcs.Energy < mcs.MAX_ENERGY) break;
    }
    mcs.STATS_MULTICANONICAL = true;
    mcs.init_stats(model);
    int nsamples = mcs.run(model, nsweeps, true);
    mcs.finish_stats(model);
    assert(nsamples == nsweeps);
    update_all(model, iter, nsamples);
    model.write_parameters("keke.param0");
  }
}

void OptParam::write_momentum(Model& model, const string filename)
{
  ofstream ofs(filename.c_str());

  ofs.setf(std::ios::scientific);
  ofs.precision(15);
  const int model_length = model.length();
  ofs << "#SingleSite" << endl;
  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      ofs << "O\t" << i << "\t" << amino1_core[a] << "\t"
	  << H_O[i][a] << endl;
    }
    if(i == model_length) continue;
    for(int a = 0; a < nstate_insert; ++a) {
      ofs << "I\t" << i << "\t" << amino1_insert[a] << "\t"
	  << H_I[i][a] << endl;
    }
  }

  ofs << "#Short-RangePair" << endl;
  for(int i = 1; i < model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ofs << "OI_s\t" << i << "\t" 
	    << amino1_core[a] << "\t" << amino1_insert[b] << "\t"
	    << J_OI[i][a][b] << endl;
      }

      for(int b = 0; b < nstate_core; ++b) {
	ofs << "OO_s\t" << i << "\t"
	    << amino1_core[a] << "\t" << amino1_core[b] << "\t"
	    << J_OO[i][a][b] << endl;
      }
    }

    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ofs << "II_s\t" << i << "\t" 
	    << amino1_insert[a] << "\t" << amino1_insert[b] << "\t"
	    << J_II[i][a][b] << endl;
      }
      for(int b = 0; b < nstate_core; ++b) {
	ofs << "IO_s\t" << i << "\t" 
	    << amino1_insert[a] << "\t" << amino1_core[b] << "\t"
	    << J_IO[i][a][b] << endl;
      }
    }
  }

  IntPairs& intpairs = model.get_intpairs();
  ofs << "#Long-RangePair\t" << intpairs.size() << endl;
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    const int i = it->i, j = it->j;
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	ofs << "OO_l\t" << i << "\t" << j << "\t" 
	    << amino1_core[a] << "\t" << amino1_core[b] << "\t"
	    << K_OO[i][j][a][b] << endl;
      }
    }
  }

  ofs.close();
  return;
}

void OptParam::read_momentum(Model& model, const string filename)
{
  ifstream ifs(filename.c_str());
  string header;
  int ii,jj;
  char aa, bb;
  const int model_length = model.length();
  
  ifs >> header;
  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      ifs >> header >> ii >> aa >> H_O[i][a];
    }
    if(i== model_length) continue;
    for(int a = 0; a < nstate_insert; ++a) {
      ifs >> header >> ii >> aa >> H_I[i][a];
    }
  }

  ifs >> header;
  for(int i = 1; i < model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ifs >> header >> ii >> aa >> bb >> J_OI[i][a][b];
      }
      for(int b = 0; b < nstate_core; ++b) {
	ifs >> header >> ii >> aa >> bb >> J_OO[i][a][b];
      }
    }

    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ifs >> header >> ii >> aa >> bb >> J_II[i][a][b];
      }
      for(int b = 0; b < nstate_core; ++b) {
	ifs >> header >> ii >> aa >> bb >> J_IO[i][a][b];
      }
    }
  }

  IntPairs& intpairs = model.get_intpairs();
  int nint;
  ifs >> header >> nint;
  for(int k = 0; k < nint; ++k) {
    int i,j;
    double x,y,z;
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	ifs >> header >> i >> j >> aa >> bb >> z;
	K_OO[i][j][a][b] = K_OO[j][i][b][a] = z;
      }
    }
  }

  ifs.close();
  return;
}

void OptParam::optimize_traj(Model& model, AnaTraj& atraj)
{
  MCseq mcs(model);
  mcs.STATS_MULTICANONICAL = false;
  mcs.stats_full_nonbonded = false;
  mcs.init_stats(model);
  int Nseqs=0;
  while(atraj.pick_seq(model, mcs)) {
    const double e = mcs.Energy;
    mcs.count_stats(model);
    Nseqs++;
    if(Nseqs % 10000 == 0) cerr << Nseqs << " sequences processed." << endl;
  }
  mcs.finish_stats(model);
  update_all(model, 1, Nseqs);

  return;
}

void OptParam::optimize_plm(Model& model, const int niter, const string param_out)
{
  const int Nseqs = model.get_num_seq();
  PartFunc pfunc(model);

  cerr << "OptParam::optimize_plm: " << Nseqs << " sequences." << endl;
  for(int iter = 1; iter <= niter; ++iter) {
    pfunc.PLM_sample(model);
    update_all(model,iter,Nseqs);
    if(iter % 10 == 0) model.write_parameters(param_out);
  }
  model.write_parameters(param_out);
}
