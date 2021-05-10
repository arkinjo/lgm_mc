#include "model.hpp"

void Model::set_null_sequence(void)
{
  std::fill(n_O.data(), n_O.data() + n_O.num_elements(), 0.0);
  std::fill(n_I.data(), n_I.data() + n_I.num_elements(), 0.0);

  std::fill(ns_OO.data(), ns_OO.data() + ns_OO.num_elements(), 0.0);
  std::fill(ns_OI.data(), ns_OI.data() + ns_OI.num_elements(), 0.0);
  std::fill(ns_IO.data(), ns_IO.data() + ns_IO.num_elements(), 0.0);
  std::fill(ns_II.data(), ns_II.data() + ns_II.num_elements(), 0.0);

  std::fill(nl_OO.data(), nl_OO.data() + nl_OO.num_elements(), 0.0);
  std::fill(nl_OI.data(), nl_OI.data() + nl_OI.num_elements(), 0.0);
  std::fill(nl_IO.data(), nl_IO.data() + nl_IO.num_elements(), 0.0);
  std::fill(nl_II.data(), nl_II.data() + nl_II.num_elements(), 0.0);

  for(int i = 1; i <= model_length; ++i) {
    n_O[i][terminal_state] = 1;
    for(int j = 1; j <= model_length; ++j) {
      nl_OO[i][j][terminal_state][terminal_state] = 1;
    }
  }

  for(int i = 1; i < model_length; ++i) {
    ns_OO[i][terminal_state][terminal_state] = 1;
  }

  add_pseudo_counts(n_O,n_I,ns_OO,ns_OI,ns_IO, ns_II,
  		    nl_OO,nl_OI,nl_IO, nl_II, true);

  return;
}

void Model::set_uniform_sequence(void)
{
  std::fill(n_O.data(), n_O.data() + n_O.num_elements(), 0.0);
  std::fill(n_I.data(), n_I.data() + n_I.num_elements(), 0.0);

  std::fill(ns_OO.data(), ns_OO.data() + ns_OO.num_elements(), 0.0);
  std::fill(ns_OI.data(), ns_OI.data() + ns_OI.num_elements(), 0.0);
  std::fill(ns_IO.data(), ns_IO.data() + ns_IO.num_elements(), 0.0);
  std::fill(ns_II.data(), ns_II.data() + ns_II.num_elements(), 0.0);

  std::fill(nl_OO.data(), nl_OO.data() + nl_OO.num_elements(), 0.0);
  std::fill(nl_OI.data(), nl_OI.data() + nl_OI.num_elements(), 0.0);
  std::fill(nl_IO.data(), nl_IO.data() + nl_IO.num_elements(), 0.0);
  std::fill(nl_II.data(), nl_II.data() + nl_II.num_elements(), 0.0);

  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      n_O[i][a] = 1.0/nstate_core;
      for(int j = 1; j <= model_length; ++j) {
	for(int b = 0; b < nstate_core; ++b) {
	  nl_OO[i][j][a][b] = 1.0/(nstate_core*nstate_core);
	}
      }
    }
  }

  for(int i = 1; i < model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	ns_OO[i][a][b] = 1.0/(nstate_core*nstate_core);
      }
    }
  }

  add_pseudo_counts(n_O,n_I,ns_OO,ns_OI,ns_IO, ns_II,
  		    nl_OO,nl_OI,nl_IO, nl_II, true);

  return;
}

double Model::get_seqlength_generic(darray2& no, darray2& ni,
				    double& lO, double& lI)
{
  lO = lI = 0.0;
  for(int i = 1; i <= model_length; i++) {
    for(int a = 0; a < nstate_match; ++a) {
      lO += no[i][a];
    }

    if(i == model_length) break;
    for(int a = 0; a < nstate_insert; ++a) {
      lI += ni[i][a];
    }
  }
  return (lO + lI);
}

double Model::get_seqlength(double& lO, double& lI)
{
  return get_seqlength_generic(n_O,n_I, lO, lI);
}

double Model::get_seqlength(void)
{
  double lO,lI;
  return get_seqlength(lO, lI);
}

double Model::get_seqlength_ref(double& lO, double& lI)
{
  return get_seqlength_generic(n0_O,n0_I, lO, lI);
}

double Model::get_seqlength_ref(void)
{
  double lO,lI;
  return get_seqlength_ref(lO, lI);
}

double Model::get_seqlength_exp(double& lO, double& lI)
{
  return get_seqlength_generic(en_O, en_I, lO, lI);
}

double Model::get_seqlength_exp(void)
{
  double lO, lI;
  return get_seqlength_exp(lO, lI);
}

double Model::get_divergence(void)
{
  std::fill(div_O.data(), div_O.data() + div_O.num_elements(),0.0);
  std::fill(ddiv_O.data(), ddiv_O.data() + ddiv_O.num_elements(),0.0);
  std::fill(mdiv_O.data(), mdiv_O.data() + mdiv_O.num_elements(),0.0);
  
  double divO = 0.0;
  const double b2 = beta*beta;
  const double ene = get_energy_exp();
  
  for(int i = 1; i <= model_length; i++) {
    for(int a = 0; a < nstate_core; ++a) {
      const double ena = en_O[i][a];
      if(ena > 0.0) {
	const double wa = log(ena/n_O[i][a]);
	div_O[i] += ena*wa;
	ddiv_O[i] += wa*(nE_O[i][a] - ena*ene)*b2;
	mdiv_O[i] += ena*(1 - ena);
      }
    }
    divO += div_O[i];
  }
  divO /= model_length;


  std::fill(div_OO.data(), div_OO.data() + div_OO.num_elements(), 0.0);
  std::fill(div_OOs.data(), div_OOs.data() + div_OOs.num_elements(), 0.0);
  std::fill(ddiv_OOs.data(), ddiv_OOs.data() + ddiv_OOs.num_elements(), 0.0);
  
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    const int i = it->i, j = it->j;
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	const double enab = enl_OO[i][j][a][b];
	if(enab > 0.0) {
	  const double wa = log(enab/nl_OO[i][j][a][b]);
	  const double div = enab*wa;
	  div_OO[i][j] += div;
	  div_OOs[i] += div;
	  div_OOs[j] += div;
	  const double ddiv = wa*(nE_OO[i][j][a][b] - enab*ene)*b2;
	  ddiv_OOs[i] += ddiv;
	  ddiv_OOs[j] += ddiv;
	}
      }
    }
    cerr << "Div_OO\t"
	 << i << "\t"
	 << j << "\t"
	 << div_OO[i][j] << endl;
  }

  return divO;
}

double Model::get_divergence_from_ref(void)
{
  return get_divergence(n0_O, div_R);
}

double Model::get_divergence(darray2& tn_O, darray1& tdiv)
{
  std::fill(tdiv.data(), tdiv.data() + tdiv.num_elements(),0.0);

  double divO = 0.0;
  const double b2 = beta*beta;
  const double ene = get_energy_exp();
  
  for(int i = 1; i <= model_length; i++) {
    for(int a = 0; a < nstate_core; ++a) {
      if(en_O[i][a] > 0.0) {
	tdiv[i] += en_O[i][a]*log(en_O[i][a]/tn_O[i][a]);
      }
    }
    divO += tdiv[i];
  }

  divO /= model_length;
  return divO;
}

double Model::get_deviation(double& devO, double& devI)
{
  std::fill(dev_O.data(), dev_O.data() + dev_O.num_elements(),0.0);
  std::fill(dev_I.data(), dev_I.data() + dev_I.num_elements(),0.0);

  devO = devI = 0.0;
  double q=0.0;  
  for(int i = 1; i <= model_length; i++) {
    q = 0.0;
    for(int a = 0; a < nstate_core; ++a) {
      double d = (en_O[i][a] - n_O[i][a]);
      q += d*d;
    }
    dev_O[i] += sqrt(q);
    devO += dev_O[i];

    if(i == model_length) break;
    q=0.0;
    for(int a = 0; a < nstate_insert; ++a) {
      double d = en_I[i][a] - n_I[i][a];
      q += d*d;
    }
    dev_I[i] += sqrt(q);
    devI += dev_I[i];
  }

  double tdev = (devO + devI) / (2.0*model_length - 1.0);
  devO /= model_length;
  devI /= (model_length-1);

  return tdev;
}

double Model::get_site_entropy_generic(darray2& tn_O)
{
  std::fill(ent_O.data(), ent_O.data() + ent_O.num_elements(), 0.0);

  double entO = 0.0;
  
  for(int i = 1; i <= model_length; i++) {
    for(int a = 0; a < nstate_core; ++a) {
      if(tn_O[i][a] > 0.0) {
	ent_O[i] += -tn_O[i][a]*log(tn_O[i][a]);
      }
    }
    entO += ent_O[i];
  }
  
  return (entO/model_length);
}

double Model::get_site_entropy(void)
{
  return get_site_entropy_generic(n_O);
}

double Model::get_site_entropy_ref(void)
{
  return get_site_entropy_generic(n0_O);
}

double Model::get_site_entropy_exp(void)
{
  return get_site_entropy_generic(en_O);
}

double Model::get_energy_generic(darray2& tn_O, darray2& tn_I,
				 darray3& tns_OO, darray3& tns_OI,
				 darray3& tns_IO, darray3& tns_II,
				 darray4& tnl_OO,
				 double &eshort, double &elong, double &egibbs,
				 darray1& tsite_energy)
{
  eshort = elong = egibbs = 0.0;
  std::fill(tsite_energy.data(), tsite_energy.data() + tsite_energy.num_elements(),0.0);

  for(int a = 0; a < nstate_core; ++a) {
    egibbs -= H_O[model_length][a]*tn_O[model_length][a];
    tsite_energy[model_length] -= H_O[model_length][a]*tn_O[model_length][a];
  }

  double e;
  for(int i = 1; i < model_length; i++) {
    for(int a = 0; a < nstate_core; ++a) {
      e = -H_O[i][a]*tn_O[i][a];
      egibbs += e;
      tsite_energy[i] += e;
      for(int b = 0; b < nstate_core; ++b) {
	e = -J_OO[i][a][b]*tns_OO[i][a][b];
	eshort += e;
	tsite_energy[i] += e;
	tsite_energy[i+1] += e;
      }
      for(int b = 0; b < nstate_insert; ++b) {
	e = -J_OI[i][a][b]*tns_OI[i][a][b];
	eshort += e;
	tsite_energy[i] += e;
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      egibbs -= H_I[i][a]*tn_I[i][a];
      for(int b = 0; b < nstate_insert; ++b) {
	eshort -= J_II[i][a][b]*tns_II[i][a][b];
      }
      for(int b = 0; b < nstate_core; ++b) {
	if(i == model_length && b != terminal_state) continue;
	e = -J_IO[i][a][b]*tns_IO[i][a][b];
	eshort += e;
	tsite_energy[i+1] += e;
      }
    }
  }

  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    int i = it->i, j = it->j;
    for(int a = 0; a < nstate_core; ++a)
      for(int b = 0; b < nstate_core; ++b) {
	e = -K_OO[i][j][a][b]*tnl_OO[i][j][a][b];
	elong += e;
	tsite_energy[i] += e;
	tsite_energy[j] += e;
      }
  }

  return (eshort + elong + egibbs);
}
  

double Model::get_energy(void)
{
  double b,n,s;
  
  return get_energy(b,n,s);
}

double Model::get_energy(double& eb, double &enb, double &eg)
{
  return get_energy_generic(n_O,n_I,ns_OO,ns_OI,ns_IO,ns_II, nl_OO,
			    eb, enb, eg, site_energy);
}

double Model::get_energy_ref(void)
{
  double b,n,s;
  return get_energy_ref(b,n,s);
}

double Model::get_energy_ref(double& eb, double &enb, double &eg)
{
  eb = enb = eg = 0.0;
  return get_energy_generic(n0_O,n0_I,ns0_OO,ns0_OI,ns0_IO,ns0_II, nl0_OO,
			    eb, enb, eg, site_energy_ref);
}


double Model::get_energy_exp(void)
{
  double b,n,g;
  return get_energy_exp(b,n,g);
}

double Model::get_energy_exp(double &eb, double &enb, double &egibbs)
{
  eb = enb = egibbs = 0.0;
  return get_energy_generic(en_O,en_I,ens_OO,ens_OI,ens_IO, ens_II, enl_OO,
			    eb,enb,egibbs, site_energy_exp);
}

double Model::chem_pot_work(void)
{
  double work = 0.0;
  
  for(int i = 0; i <= model_length; i++) {
    for(int a = 0; a < nstate_insert; ++a) {
      work += (en_I[i][a] - n0_I[i][a])*H_I[i][a];
    }
    if(i>0) {
      for(int a = 0; a < nstate_core; ++a) {
	work += (en_O[i][a] - n0_O[i][a])*H_O[i][a];
      }
    }
  }

  return work;
}
