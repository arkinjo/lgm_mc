#include "model.hpp"

void Model::check_const_site(void)
{
  cout.setf(std::ios::scientific);
  for(int i = 1; i <= model_length; ++i) {
    double t = 0.0;
    for(int a = 0; a < nstate_core; ++a) {
      t += n_O[i][a];
    }
    cout << "check_const_site: " << i << "\t" << t << endl;
  }
  return;
}

void Model::check_const_site_exp(void)
{
  cout.setf(std::ios::scientific);
  for(int i = 1; i <= model_length; ++i) {
    double t = 0.0;
    for(int a = 0; a < nstate_core; ++a) {
      t += en_O[i][a];
    }
    cout << "check_const_site_exp: " << i << "\t" << t << endl;
  }
  return;
}

void Model::check_const_bonded(void) {
  cout.setf(std::ios::scientific);
  for(int i = 0; i <= model_length; ++i) {
    double t;
    for(int a = 0; a < nstate_core; ++a) {
      if(i == 0 && a != terminal_state) continue;
      t = 0.0;
      for(int b = 0; b < nstate_core; ++b) {
	t += ns_OO[i][a][b];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	t += ns_OI[i][a][b];
      }
      cout << "check_const_bonded(O+): " << i << "\t" << amino1_core[a]
	   << "\t" << t << "\t" << n_O[i][a]
	   << "\t" << t - n_O[i][a] << endl;
    }

    for(int b = 0; b < nstate_core; ++b) {
      if(i == model_length && b != terminal_state) continue;
      t = 0.0;
      for(int a = 0; a < nstate_core; ++a) {
	if(i == 0 && a != terminal_state) continue;
	t += ns_OO[i][a][b];
      }
      for(int a = 0; a < nstate_insert; ++a) {
	t += ns_IO[i][a][b];
      }
      cout << "check_const_bonded(O-): " << i << "\t" << amino1_core[b]
	   << "\t" << t << "\t" << n_O[i+1][b]
	   << "\t" << t - n_O[i+1][b] << endl;
    }


    for(int a = 0; a < nstate_insert; ++a) {
      t = 0.0;
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	t += ns_IO[i][a][b];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	t += ns_II[i][a][b];
      }
      cout << "check_const_bonded(I+): " << i << "\t" << amino1_insert[a]
	   << "\t" << t << "\t" << n_I[i][a]
	   << "\t" << t - n_I[i][a] << endl;
    }

    for(int b = 0; b < nstate_insert; ++b) {
      t = 0.0;
      for(int a = 0; a < nstate_core; ++a) {
	if(i == 0 && a != terminal_state) continue;
	t += ns_OI[i][a][b];
      }
      for(int a = 0; a < nstate_insert; ++a) {
	t += ns_II[i][a][b];
      }
      cout << "check_const_bonded(I-): " << i << "\t" << amino1_insert[b]
	   << "\t" << t << "\t" << n_I[i][b]
	   << "\t" << t - n_I[i][b] << endl;
    }
    cout << "//" << endl;
  }

  return;
}

void Model::check_const_bonded_exp(void) {
  cout.setf(std::ios::scientific);
  double t,t1,t2;
  
  for(int i = 0; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      if(i == 0 && a != terminal_state) continue;
      t1 = t2 = 0.0;
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	t1 += ens_OO[i][a][b];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	t2 += ens_OI[i][a][b];
      }
      t = t1 + t2;
      cout << "check_const_bonded_exp(O+): " << i << "\t" << amino1_core[a] << "\t" 
	   << t << "\t" << en_O[i][a] 
	   << "\t" << t - en_O[i][a]
	   << endl;
    }

    for(int b = 0; b < nstate_core; ++b) {
      if(i == model_length && b != terminal_state) continue;
      t1 = t2 = 0.0;
      for(int a = 0; a < nstate_core; ++a) {
	if(i == 0 && a != terminal_state) continue;
	t1 += ens_OO[i][a][b];
      }
      for(int a = 0; a < nstate_insert; ++a) {
	t2 += ens_IO[i][a][b];
      }
      t = t1 + t2;
      cout << "check_const_bonded_exp(O-): " << i << "\t" << amino1_core[b] << "\t" 
	   << t << "\t" << en_O[i+1][b] 
	   << "\t" << t - en_O[i+1][b]
	   << endl;
    }

    for(int a = 0; a < nstate_insert; ++a) {
      t1 = t2 = 0.0;
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	t1 += ens_IO[i][a][b];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	t2 += ens_II[i][a][b];
      }
      t = t1 + t2;
      cout << "check_const_bonded_exp(I+): " << i << "\t" << amino1_insert[a] << "\t" 
	   << t << "\t" << en_I[i][a]
	   << "\t" << t - en_I[i][a]
	   << endl;
    }

    for(int b = 0; b < nstate_insert; ++b) {
      t1 = t2 = 0.0;
      for(int a = 0; a < nstate_core; ++a) {
	if(i == 0 && a != terminal_state) continue;
	t1 += ens_OI[i][a][b];
      }
      for(int a = 0; a < nstate_insert; ++a) {
	t2 += ens_II[i][a][b];
      }
      t = t1 + t2;
      cout << "check_const_bonded_exp(I-): " << i << "\t" << amino1_insert[b] << "\t" 
	   << t << "\t" << en_I[i][b]
	   << "\t" << t - en_I[i][b]
	   << endl;
    }
    cout << "//" << endl;
  }

  return;
}

void Model::check_const_nonbonded(void)
{
  cout.setf(std::ios::scientific);
  for(int i = 1; i <= model_length; ++i) {
    for(int j = 0; j <= model_length; ++j) {
      for(int b = 0; b < nstate_insert; ++b) {
	double t = 0.0;
	for(int a = 0; a < nstate_core; ++a) {
	  t += nl_OI[i][j][a][b];
	}
	cout << "check_const_nonbonded(I+) "
	     << i << "\t" << j << "\t" << amino1_insert[b] << "\t"
	     << t << "\t" << n_I[j][b]
	     << "\t" << t - n_I[j][b] << endl;

      }
    }
    for(int j = i; j <= model_length; ++j) {
      for(int a = 0; a < nstate_core; ++a) {
	double t = 0.0;
	for(int b = 0; b < nstate_core; ++b) {
	  t += nl_OO[i][j][a][b];
	}
	cout << "check_const_nonbonded(O+) "
	     << i << "\t" << j << "\t" << amino1_core[a] << "\t"
	     << t << "\t" << n_O[i][a]
	     << "\t" << t - n_O[i][a] << endl;
      }
      for(int b = 0; b < nstate_core; ++b) {
	double t = 0.0;
	for(int a = 0; a < nstate_core; ++a) {
	  t += nl_OO[i][j][a][b];
	}
	cout << "check_const_nonbonded(O-) "
	     << i << "\t" << j << "\t" << amino1_core[b] << "\t"
	     << t << "\t" << n_O[j][b]
	     << "\t" << t - n_O[j][b] << endl;
      }
    }
  }
  cout << "//" << endl;
  return;
}

void Model::check_const_nonbonded_exp(void)
{
  cout.setf(std::ios::scientific);
  for(int i = 1; i <= model_length; ++i) {
    for(int j = 0; j <= model_length; ++j) {
      for(int b = 0; b < nstate_insert; ++b) {
	double t = 0.0;
	for(int a = 0; a < nstate_core; ++a) {
	  t += enl_OI[i][j][a][b];
	}
	cout << "check_const_nonbonded(I+) "
	     << i << "\t" << j << "\t" << amino1_insert[b] << "\t"
	     << t << "\t" << en_I[j][b]
	     << "\t" << t - en_I[j][b] << endl;

      }
    }
    for(int j = i; j <= model_length; ++j) {
      for(int a = 0; a < nstate_core; ++a) {
	double t = 0.0;
	for(int b = 0; b < nstate_core; ++b) {
	  t += enl_OO[i][j][a][b];
	}
	cout << "check_const_nonbonded(O+) "
	     << i << "\t" << j << "\t" << amino1_core[a] << "\t"
	     << t << "\t" << en_O[i][a]
	     << "\t" << t - en_O[i][a] << endl;
      }
      for(int b = 0; b < nstate_core; ++b) {
	double t = 0.0;
	for(int a = 0; a < nstate_core; ++a) {
	  t += enl_OO[i][j][a][b];
	}
	cout << "check_const_nonbonded(O-) "
	     << i << "\t" << j << "\t" << amino1_core[b] << "\t"
	     << t << "\t" << en_O[j][b]
	     << "\t" << t - en_O[j][b] << endl;
      }
    }
  }
  cout << "//" << endl;
  return;
}

double Model::check_site(const bool show)
{
  cout.setf(std::ios::scientific);
  double tot=0;
  double n,en,del;
  double ntot = 0;
  double max_del = -1, adel;
  for(int i = 0; i <= model_length; ++i) {
    if(i>0) {
      for(int a = 0; a < nstate_core; ++a) {
	n = n_O[i][a];
	en = en_O[i][a];
	del = en - n;
	adel = fabs(del);
	max_del = (adel > max_del) ? adel : max_del;
	tot += del*del;
	ntot += n*n;
	if(show) {
	  cerr << "O\t" << i << "\t" << a << "\t" 
	       << n_O[i][a] << "\t" 
	       << en << "\t" 
	       << del << "\t" 
	       << H_O[i][a]
	       << endl;
	}
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      n = n_I[i][a];
      en = en_I[i][a];
      del = en - n;
      adel = fabs(del);
      max_del = (adel > max_del) ? adel : max_del;
      tot += del*del;
      ntot += n*n;
      if(show) {
	cerr << "I\t" << i << "\t" << a << "\t"
	     << n_I[i][a] << "\t" 
	     << en << "\t" 
	     << del << "\t" 
	     << H_I[i][a]
	     << endl;
      }
    }
  }
  return max_del;
}

double Model::check_bonded(const bool show)
{
  cout.setf(std::ios::scientific);
  double tot=0.0;
  double n,en, del;
  double ntot = 0.0;
  double max_del=-1, adel;
  for(int i = 0; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      if(i == 0 && a != terminal_state) continue;

      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	n = ns_OO[i][a][b];
	en = ens_OO[i][a][b];
	del = en - n;
	adel = fabs(del);
	if(max_del < adel) max_del = adel; 
	tot += del*del;
	ntot += n*n;
	if(show) {
	  cerr << "OO_s\t" << i << "\t" << a << "," << b << "\t"
	       << n << "\t"
	       << en << "\t"
	       << del << "\t"
	       << J_OO[i][a][b]
	       << endl;
	}
      }
      for(int b = 0; b < nstate_insert; ++b) {
	n = ns_OI[i][a][b];
	en = ens_OI[i][a][b];
	del = en - n;
	adel = fabs(del);
	if(max_del < adel) max_del = adel; 
	tot += del*del;
	ntot += n*n;
	if(show) {
	  cerr << "OI_s " << i << "\t" << a << "," << b << "\t"
	       << n << "\t"
	       << en << "\t"
	       << del << "\t"
	       << J_OI[i][a][b]
	       << endl;
	}
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	n = ns_IO[i][a][b];
	en = ens_IO[i][a][b];
	del = en - n;
	adel = fabs(del);
	if(max_del < adel) max_del = adel; 
	tot += del*del;
	ntot += n*n;
	if(show) {
	  cerr << "IO_s\t" << i << "\t" << a << "," << b << "\t"
	       << n << "\t"
	       << en << "\t"
	       << del << "\t"
	       << J_IO[i][a][b]
	       << endl;
	}
      }
      for(int b = 0; b < nstate_insert; ++b) {
	n = ns_II[i][a][b];
	en = ens_II[i][a][b];
	del = en - n;
	adel = fabs(del);
	if(max_del < adel) max_del = adel; 
	tot += del*del;
	ntot += n*n;
	if(show) {
	  cerr << "II_s\t" << i << "\t" << a << "," << b << "\t"
	       << n << "\t"
	       << en << "\t"
	       << del << "\t"
	       << J_II[i][a][b]
	       << endl;
	}
      }
    }
  }

  // cout << "site_nn total " << tot << " / " << ntot << endl;
  return max_del;
}

double Model::check_nonbonded(const bool show)
{
  cout.setf(std::ios::scientific);
  double tot=0.0;
  double n,en, del;
  double ntot = 0.0;
  double max_del=-1, adel;
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    const int i = it->i, j = it->j;
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	n = nl_OO[i][j][a][b];
	en = enl_OO[i][j][a][b];
	del = en - n;
	adel = fabs(del);
	if(max_del < adel) max_del = adel; 
	tot += del*del;
	ntot += n*n;
	if(show) {
	  cerr << "OO_l\t" << i << "\t" << j
	       << "\t" << a << "\t" << b << "\t"
	       << nl_OO[i][j][a][b] << "\t"
	       << en << "\t"
	       << del << "\t"
	       << K_OO[i][j][a][b]
	       << endl;
	}
      }
    }
  }

  // cout << "site_nn total " << tot << " / " << ntot << endl;
  return max_del;
}


