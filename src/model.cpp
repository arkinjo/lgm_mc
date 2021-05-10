#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "msa.hpp"
#include "model.hpp"


static const int min_seq_sep = 6;
// used in set_intpairs().
// Those pairs with |i-j| >= min_seq_sep are included in "intpairs".

static const int ADD_PSEUDOCOUNT = 1;

SSite::SSite() 
{
  O = '-';
  I="";
}

IPair::IPair(int i0, int j0)
{
  i = i0; j = j0;
}

Model::Model() : temperature(1.0), beta(1.0), Mu_len(0.0) {}
Model::~Model() {}

char normal_amino1(const char a)
{
  switch (a) {
  case 'B':
    return 'N';
  case 'Z':
    return 'Q';
  case 'J':
    return 'L';
  case 'X':
    return 'A';
  case 'U':
    return 'C';
  case 'O':
    return 'K';
  case 'b':
    return 'n';
  case 'z':
    return 'q';
  case 'j':
    return 'l';
  case 'x':
    return 'a';
  case 'u':
    return 'c';
  case 'o':
    return 'k';
  default:
    return a;
  }
}

int find_amino1_core(const char a)
{
  char aa = normal_amino1(a);
    
  string::size_type i = amino1_core.find(aa);
  if(i == string::npos) {
    cerr << "find_amino1_core(" << a << ")=" << i << endl;
    i = 0; // arbitrary
  }
  return i;
}

int find_amino1_insert(const char a)
{
  char aa = normal_amino1(a);
  
  string::size_type i = amino1_insert.find(aa);
  if(i == string::npos) {
    cerr << "find_amino1_insert(" << a << ")=" << i << endl;
    i = 0; // arbitrary
  }
  return i;
}

bool amino1_core_p(const char a)
{
  char aa = normal_amino1(a);
  string::size_type i = amino1_core.find(aa);
  return ((i == string::npos) ? false : true);
}

bool amino1_insert_p(const char a)
{
  char aa = normal_amino1(a);
   string::size_type i = amino1_insert.find(aa);
   return ((i == string::npos) ? false : true);
}

void Model::assign_sites_1st(const Msa& msa) 
{
  vector<double> count(msa_length,0);

  string seq = msa.sequence(0);
  for (int i = 0; i < msa_length; ++i) {
    if(seq[i] != '-') {
      count[i] += 1;
    }
  }
  msa_sites.clear();
  for(int i = 0; i < msa_length; ++i) {
    msa_sites.push_back(((double)count[i] >0));
  }

}

void Model::assign_sites(const Msa& msa) 
{
  vector<double> count(msa_length,0);

  double tot_w = 0.0;
  for (int j = 0; j < num_seq; ++j) {
    string seq = msa.sequence(j);
    double w = msa.get_weight(j);
    tot_w += w;
    for (int i = 0; i < msa_length; ++i) {
      if(seq[i] != '-') {
	count[i] += w;
      }
    }
  }

  double maxcnt = 0.0;
  for (int i = 0; i < msa_length; ++i) {
    count[i] /= tot_w;
    maxcnt = (maxcnt < count[i]) ? count[i] : maxcnt;
  }
  cerr << "assign_sites max count = " << maxcnt << " / " << num_seq << endl;

  msa_sites.clear();
  for(int i = 0; i < msa_length; ++i) {
    msa_sites.push_back(((double)count[i] >= match_fraction));
  }

}

void Model::assign_sites_hmm(const Msa& msa) 
{
  msa_sites.clear();
  for (int i = 0; i < msa_length; ++i) {
    msa_sites.push_back(false);
  }
  for (int j = 0; j < num_seq; ++j) {
    string seq = msa.sequence(j);
    for (int i = 0; i < msa_length; ++i) {
      string::size_type ind = amino1_core.find(seq[i]);
      if(ind != string::npos) {
	msa_sites[i] = true;
      }
    }
  }

}

modelSeq Model::seq_to_modelSeq(const string &seq)
{
  modelSeq ms(model_length+2);

  for(int i = 0; i <= model_length + 1; ++i) {
    ms[i].O = '-';
    ms[i].I = "";
  }

  int isite = 0;
  for(int i = 0; i < msa_length; ++i) {
    if(msa_sites[i]) {
      ++isite;
      ms[isite].O = seq[i];
    }
    else if(seq[i] != '.' && seq[i] !='-') {
      if(isite > 0 && isite < model_length) {
	// N & C terminal insertions are trimmed!!!!
	//	if(ms[isite].I.size() < 10) // long insertions are truncated!!!
	  ms[isite].I.push_back(tolower(seq[i]));
      }
    }
  }

  ms[0].O = ms[model_length+1].O = terminal_residue;

  return ms;
}

modelSeq Model::seq_to_modelSeq_blast(const Msa& msa,const int iseq)
{
  modelSeq ms(model_length+2);

  for(int i = 0; i <= model_length + 1; ++i) {
    ms[i].O = '-';
    ms[i].I = "";
  }

  const int qs = msa.qstart(iseq);
  const int qe = msa.qend(iseq);
  const string qseq = msa.qseq(iseq);
  const string sseq = msa.sseq(iseq);
  const int len = qseq.length();
  int isite=qs-1;
  for(int i = 0; i < len; ++i) {
    if(qseq[i] != '-') {
      isite++;
      ms[isite].O = sseq[i];
    } else {
      ms[isite].I.push_back(tolower(sseq[i]));
    }
  }
  ms[0].O = ms[model_length+1].O = terminal_residue;

  return ms;
}

void Model::init_params(const double ep)
{
  std::fill(H_O.data(), H_O.data() + H_O.num_elements(), 0.0);
  std::fill(H_I.data(), H_I.data() + H_I.num_elements(), 0.0);
  std::fill(J_OO.data(), J_OO.data() + J_OO.num_elements(), 0.0);
  std::fill(J_OI.data(), J_OI.data() + J_OI.num_elements(), 0.0);
  std::fill(J_IO.data(), J_IO.data() + J_IO.num_elements(), 0.0);
  std::fill(J_II.data(), J_II.data() + J_II.num_elements(), 0.0);
  std::fill(K_OO.data(), K_OO.data() + K_OO.num_elements(), 0.0);
  return;
}

void Model::allocate(const int mlen)
{
  model_length = mlen;
  darray1::extent_gen ext1;
  site_energy.resize(ext1[model_length+1]);
  site_energy_exp.resize(ext1[model_length+1]);
  site_energy_ref.resize(ext1[model_length+1]);
  
  div_O.resize(ext1[model_length+1]);
  ddiv_O.resize(ext1[model_length+1]);

  div_OOs.resize(ext1[model_length+1]);
  ddiv_OOs.resize(ext1[model_length+1]);

  mdiv_O.resize(ext1[model_length+1]);
  div_R.resize(ext1[model_length+1]);

  dev_O.resize(ext1[model_length+1]);
  dev_I.resize(ext1[model_length+1]);

  ent_O.resize(ext1[model_length+1]);
  len_I.resize(ext1[model_length+1]);
  darray2::extent_gen ext2;

  // observed number densities
  n_O.resize(ext2[model_length+2][nstate]);
  n_I.resize(ext2[model_length+2][nstate]);
  n_Int.resize(ext2[model_length+2][nstate]);
  n_Ict.resize(ext2[model_length+2][nstate]);
  n0_O.resize(ext2[model_length+2][nstate]);
  n0_I.resize(ext2[model_length+2][nstate]);
  // log of expected number densities
  en_O.resize(ext2[model_length+2][nstate]);
  en_I.resize(ext2[model_length+2][nstate]);
  // 1-site parameters (chemical potentials)
  H_O.resize(ext2[model_length+2][nstate]);
  H_I.resize(ext2[model_length+2][nstate]);

  nE_O.resize(ext2[model_length+2][nstate]);
  nE_I.resize(ext2[model_length+2][nstate]);

  div_OO.resize(ext2[model_length+1][model_length+1]);

  darray4::extent_gen ext4;
  // number densities
  ns_OO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ns_OI.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ns_IO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ns_II.resize(ext4[model_length+2][model_length+2][nstate][nstate]);

  // reference densities
  ns0_OO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ns0_OI.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ns0_IO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ns0_II.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  // expected densities
  ens_OO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ens_OI.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ens_IO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  ens_II.resize(ext4[model_length+2][model_length+2][nstate][nstate]);

  // coupling constants.
  J_OO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  J_OI.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  J_IO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  J_II.resize(ext4[model_length+2][model_length+2][nstate][nstate]);


  nl_OO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  nl_OI.resize(ext4[model_length+2][model_length+1][nstate][nstate]);
  nl_IO.resize(ext4[model_length+1][model_length+2][nstate][nstate]);
  nl_II.resize(ext4[model_length+1][model_length+1][nstate][nstate]);
  nl0_OO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  nl0_OI.resize(ext4[model_length+2][model_length+1][nstate][nstate]);
  nl0_IO.resize(ext4[model_length+1][model_length+2][nstate][nstate]);
  nl0_II.resize(ext4[model_length+1][model_length+1][nstate][nstate]);

  enl_OO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);
  enl_OI.resize(ext4[model_length+2][model_length+1][nstate][nstate]);
  enl_IO.resize(ext4[model_length+1][model_length+2][nstate][nstate]);
  enl_II.resize(ext4[model_length+1][model_length+1][nstate][nstate]);
  K_OO.resize(ext4[model_length+2][model_length+2][nstate][nstate]);

  nE_OO.resize(ext4[model_length+1][model_length+1][nstate][nstate]);
}

void Model::sum_I(void)
{
  std::fill(len_I.data(), len_I.data() + len_I.num_elements(), 0.0);
  std::fill(n_Int.data(), n_Int.data() + n_Int.num_elements(), 0.0);
  std::fill(n_Ict.data(), n_Ict.data() + n_Ict.num_elements(), 0.0);
  for(int i = 1; i < model_length; ++i) {
    for(int a = 0; a < nstate; ++a) {
      len_I[i] += n_I[i][a];
      for(int b = 0; b < nstate; ++b) {
	n_Int[i][a] += ns_OI[i][b][a];
	n_Ict[i][a] += ns_IO[i][a][b];
      }
    }
    cerr << "Insert_length: " << i << "\t" << len_I[i] << endl;
  }
  return;
}

void Model::set_temperature(const double t)
{
  temperature = t;
  beta = 1.0/temperature;
  return;
}

void Model::set_lambda(void)
{
  cerr << "# set_lambda called..." << endl;
  lambdaH_O = GAMMA/nstate;
  lambdaH_I = GAMMA/nstate;
  
  lambdaJ_OO = GAMMA/(2*nstate*nstate);
  lambdaJ_II = GAMMA/(2*nstate*nstate);
  lambdaJ_OI = lambdaJ_IO = GAMMA/(2*nstate*nstate);
  

  lambdaK_OO = GAMMA / (nstate*nstate);
  lambdaK_OI = GAMMA / (nstate*nstate);
  lambdaK_II = GAMMA / (nstate*nstate);

  const double C0 = 1.0/(GAMMA + 1.0);
  cs0_OO = C0*(lambdaJ_OO - C0*lambdaH_O*lambdaH_O);
  cs0_OI = C0*(lambdaJ_OI - C0*lambdaH_O*lambdaJ_OI*nstate);
  cs0_IO = cs0_OI;
  cs0_II = C0*(lambdaJ_II - C0*lambdaH_I*lambdaH_I);

  cl0_OO = C0*(lambdaK_OO - C0*lambdaH_O*lambdaH_O);
  cl0_OI = C0*(lambdaK_OI - C0*lambdaH_O*lambdaH_I);
  cl0_IO = cl0_OI;
  cl0_II = C0*(lambdaK_II - C0*lambdaH_I*lambdaH_I);
  
  return;
}

Model::Model(const Msa& msa, const bool hmm, const bool blast)
{
  msa_length = msa.length();
  num_seq = msa.num_seq();
  temperature = beta = 1.0;
  cerr << "Number of sequences in MSA: " << num_seq << endl;

  //default regularization parameters
  GAMMA = DEFAULT_GAMMA;

  set_lambda();
  // NS_CUTOFF= 0.001;  NL_CUTOFF= 0.001;
  NS_CUTOFF= 5.0;  NL_CUTOFF= 5.0;

  if(hmm) {
    assign_sites_hmm(msa);
  } else if(blast) {
    ; // do nothing! See seq_to_modelSeq_blast
  } else {
    assign_sites(msa);
  }

  if(blast) {
    model_length= msa.get_base_length();
  } else {
    model_length = 0;
    for(int i = 0; i < msa_length; ++i)
      if(msa_sites[i]) ++model_length;
  }
  
  cerr << "Model length: " << model_length << endl;

  allocate(model_length);
  set_mu_len(0.0);

  PreModel.resize(num_seq);
  weights.resize(num_seq);
  if(blast) {
    for(int m = 0; m < num_seq; ++m) {
      PreModel[m] = seq_to_modelSeq_blast(msa,m);
    }
  } else {
    for(int m = 0; m < num_seq; ++m) {
      PreModel[m] = seq_to_modelSeq(msa.sequence(m));
    }
  }
  set_weights();

  stat_all();
  sum_I();
  init_params();

  SeqLength = get_seqlength();
}

void Model::set_nl_cutoff(const int ninc)
{
  // NL_CUTOFF = 0.001;
  // return;
  ////////////////////////////////////////////////////
  vector<double> v;
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    const int i = it->i, j = it->j;
    for(int a = 0; a < nstate; ++a) {
      for(int b = 0; b < nstate; ++b) {
	double corr = fabs(nl_OO[i][j][a][b] - n_O[i][a]*n_O[j][b]);
	v.push_back(corr);
      }
    }
  }
  std::sort(v.begin(), v.end());
  const int n = v.size();
  cerr << "set_nl_cutoff: n,ninc= " << n << "\t"
       << ninc << "\t" << v[n-1] << endl;
  if(ninc < n) {
    NL_CUTOFF = v[n-ninc-1];
  } else {
    NL_CUTOFF = 2*lambdaK_OO/(1.0 + GAMMA);
  }
  cerr << "# NL_CUTOFF = " << NL_CUTOFF << endl;
}

void Model::set_mu_len(const double mu)
{
  Mu_len = mu;
  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate_match; ++a) {
      H_O[i][a] += Mu_len;
    }
    if(i== model_length) break;
    for(int a = 0; a < nstate; ++a) {
      H_I[i][a] += Mu_len;
    }
  }
  return;
}

void Model::set_mu_lenI(const double mu)
{
  for(int i = 1; i < model_length; ++i) {
    for(int a = 0; a < nstate; ++a) {
      H_I[i][a] += mu;
    }
  }
  return;
}

void Model::set_weights(void)
{
  /*
    Setting weights according to:
    *** Henikoff & Henikoff, J. Mol. Biol. (1994) 243:574-578. ***
    This paper is a gem!
   */
  darray1 nrestyp(boost::extents[model_length+1]);
  darray2 nres(boost::extents[model_length+1][nstate]);
  for(int k = 0; k < num_seq; ++k) {
    for(int i = 1; i <= model_length; ++i) {
      if(PreModel[k][i].O == '-') continue;
      int a = find_amino1_core(PreModel[k][i].O);

      if (nres[i][a] == 0.0) nrestyp[i] += 1.0;
      nres[i][a] += 1.0;
    }
  }
  
  double total_weight = 0.0;
  for(int k = 0; k < num_seq; ++k) {
    double w=0.0;
    for(int i = 1; i <= model_length; ++i) {
      if(PreModel[k][i].O == '-') continue;
      int a = find_amino1_core(PreModel[k][i].O);
      w += 1.0/(nrestyp[i]*nres[i][a]);
    }
    weights[k] = w/model_length;
    total_weight += weights[k];
  }

  //total_weight should be 1.0! Nevertheless, let's normalize.
  for(int k = 0; k < num_seq; ++k) {
    weights[k] /= total_weight;
  }
  return; 
}

void Model::stat_all(void)
{
  stat_site(weights);
  stat_site_bonded(weights);
  stat_site_nonbonded(weights);
  if(ADD_PSEUDOCOUNT) {
    add_pseudo_counts(n_O,n_I,ns_OO,ns_OI,ns_IO, ns_II,
		      nl_OO,nl_OI,nl_IO, nl_II, true);
  }
  return;
}

void Model::stat_site(const vector<double>& weights)
{
  std::fill(n_O.data(), n_O.data() + n_O.num_elements(), 0.0);
  std::fill(n_I.data(), n_I.data() + n_I.num_elements(), 0.0);

  darray1 ins_count(boost::extents[nstate]);
  for(int m = 0; m < num_seq; ++m) {
    double w = weights[m];
    modelSeq mseq = PreModel[m];

    for(int i = 0; i <= model_length; ++i) {
      // Matching / Deletion
      if(i > 0) {
	string::size_type ind = find_amino1_core(mseq[i].O);
	n_O[i][ind] += w;
      }
      // Insertion
      std::fill(ins_count.data(), ins_count.data() + ins_count.num_elements(),0.0);
      for(int j = 0; j < mseq[i].I.length(); ++j) {
	string::size_type ind = find_amino1_insert(mseq[i].I[j]);
	ins_count[ind] += 1;
      }
      for(int a = 0; a < nstate; ++a) {
	n_I[i][a] += w*ins_count[a];
      }
    }
  }

  n0_O = n_O;
  n0_I = n_I;
  return;
}

void Model::stat_site_bonded(const vector<double>& weights)
{
  std::fill(ns_OO.data(), ns_OO.data() + ns_OO.num_elements(), 0.0);
  std::fill(ns_OI.data(), ns_OI.data() + ns_OI.num_elements(), 0.0);
  std::fill(ns_IO.data(), ns_IO.data() + ns_IO.num_elements(), 0.0);
  std::fill(ns_II.data(), ns_II.data() + ns_II.num_elements(), 0.0);
  
  darray2 ins_count(boost::extents[nstate][nstate]);

  for(int m = 0; m < num_seq; ++m) {
    double w = weights[m];
    modelSeq mseq = PreModel[m];    
    string::size_type a;
    string::size_type b;

    for(int i = 0; i <= model_length; ++i) {
      if(mseq[i].I.length() > 1) {
	std::fill(ins_count.data(), 
		  ins_count.data() + ins_count.num_elements(),0.0);
	for(int j = 0; j < mseq[i].I.length() - 1; ++j) {
	  a = find_amino1_insert(mseq[i].I[j]);
	  b = find_amino1_insert(mseq[i].I[j+1]);
	  ins_count[a][b] += 1;
	}
	for(int a = 0; a < nstate; ++a) {
	  for(int b = 0; b < nstate; ++b) {
	    ns_II[i][a][b] += w*ins_count[a][b];
	  }
	}
      }

      if(mseq[i].I != "") {
	int li = mseq[i].I.length() - 1;

	a = find_amino1_insert(mseq[i].I[li]);
	b = find_amino1_core(mseq[i+1].O);
	ns_IO[i][a][b] += w;

	a = find_amino1_core(mseq[i].O);
	b = find_amino1_insert(mseq[i].I[0]);
	ns_OI[i][a][b] += w;
      }
      else { //if(mseq[i].I == "")
	a = find_amino1_core(mseq[i].O);
	b = find_amino1_core(mseq[i+1].O);
	ns_OO[i][a][b] += w;
      }
    }
  }

  return;
}

void Model::stat_site_nonbonded(const vector<double>& weights)
{
  std::fill(nl_OO.data(), nl_OO.data() + nl_OO.num_elements(), 0.0);
  std::fill(nl_OI.data(), nl_OI.data() + nl_OI.num_elements(), 0.0);
  std::fill(nl_IO.data(), nl_IO.data() + nl_IO.num_elements(), 0.0);
  std::fill(nl_II.data(), nl_II.data() + nl_II.num_elements(), 0.0);

  for(int m = 0; m < num_seq; ++m) {
    double w = weights[m];
    modelSeq mseq = PreModel[m];
    for(int i = 0; i <= model_length; ++i) {
      string insa = mseq[i].I;
      for(int j = i; j <= model_length; ++j) {
	string insb = mseq[j].I;
	for(string::iterator isa = insa.begin(); isa != insa.end(); ++isa) {
	  string::size_type a = find_amino1_insert(*isa);
	  for(string::iterator isb = insb.begin(); isb != insb.end(); ++isb) {
	    string::size_type b = find_amino1_insert(*isb);
	    nl_II[i][j][a][b] += w;
	    nl_II[j][i][b][a] = nl_II[i][j][a][b];
	  }
	}
      }

      if(i == 0) continue;
      string::size_type a = find_amino1_core(mseq[i].O);
      nl_OO[i][i][a][a] += w;
      for(int j = i+1; j <= model_length; ++j) {
	string::size_type b = find_amino1_core(mseq[j].O);
	nl_OO[i][j][a][b] += w;
      }
      for(int j = 0; j <= model_length; ++j) {
	string ins = mseq[j].I;
	for(string::iterator is = ins.begin(); is != ins.end(); ++is) {
	  string::size_type b = find_amino1_insert(*is);
	  nl_OI[i][j][a][b]  += w;
	}
      }
    }
  }

  return;
}


void Model::read_ipairs(const string filename)
{
  ifstream fin;
  cerr << "Read_ipairs: " << filename << endl;
  fin.open(filename.c_str());
  if (!fin.is_open()) {
    cerr << "Read_ipairs: could not open: " << filename << endl;
    exit(1);
  }

  string pdbid, asym_id;
  int i, j;
  while(fin >> pdbid >> asym_id >> i >> j) {
    IPair ip(i,j);
    intpairs_ref.push_back(ip);
    //    cerr << "Model::read_ipairs: " << i << "\t" << j << endl;
  }

  fin.close();

}

void Model::set_default_ipairs(void)
{
  for(int i = 1; i <= model_length - min_seq_sep; ++i) {
    for(int j = i + min_seq_sep; j <= model_length; ++j) {
      IPair ip(i,j);
      intpairs.push_back(ip);
    }      
  }
  return;
}

void Model::read_reference_3Dalignment(const string filename)
{
  // This reads an alignment produced by HMMer3 in the Stockholm format.
  ifstream fin;
  cerr << "read_reference_3Dalignment: " << filename << endl;
  fin.open(filename.c_str());
  if (!fin.is_open()) {
    cerr << "read_reference_3Dalignment: could not open: " << filename << endl;
    exit(1);
  }

  string line,head,subhead;
  string refseq=(""), modseq(""), subseq;
  while(getline(fin,line)){
    if(line=="") continue;
    
    istringstream iss(line);
    iss >> head;
    if(head == "#" || head == "//") continue;
    if(head == "#=GC") {
      iss >> subhead;
      if(subhead == "RF") {
	iss >> subseq;
	modseq += subseq;
      }
      else continue;
    }
    else if (head == "#=GR") continue;
    else {
      iss >> subseq;
      refseq += subseq;
    }
  }
  fin.close();

  cerr << refseq.length() << "\t" << refseq << endl;
  cerr << modseq.length() << "\t" << modseq << endl;

  assert(refseq.length() == modseq.length());
  
  int i = 0, j = 0;
  const int len = refseq.length();
  for(int k = 0; k < len; ++k) {
    if(refseq[k] != '-') ++i;
    if(modseq[k] == 'x') ++j;
    if(refseq[k] != '-' && modseq[k] == 'x') {
      intref2mod[i] = j;
    }
  }

  set_intpairs();
}

void Model::align_ipairs_1st(Msa& msa)
{
  vector<double> count(msa_length,0);

  string seq = msa.sequence(0);
  int i = 0, j = 0;
  for (int k = 0; k < msa_length; ++k) {
    if(seq[k] != '-') i++;
    if(msa_sites[k]) j++;
    if(seq[k] != '-' && msa_sites[k]) {
      intref2mod[i] = j;
    }
  }
  set_intpairs();
}

void Model::set_intpairs(void)
{
  cerr << "Model::set_intpairs: min_seq_sep= " << min_seq_sep << endl;

  darray1 cn(boost::extents[model_length+1]);
  std::fill(cn.data(), cn.data() + cn.num_elements(), 0.0);
  
  const int nint = intpairs_ref.size();
  intpairs.clear();
  int n=0;
  for(int k = 0; k < nint; ++k) {
    const int i = intpairs_ref[k].i;
    const int j = intpairs_ref[k].j;
    if(intref2mod.count(i) > 0 && intref2mod.count(j) > 0) {
      IPair ip(intref2mod[i],intref2mod[j]);
      const bool flg = (abs(ip.j - ip.i) >= min_seq_sep);
      cerr << "ipair\t" << ++n << "\t" << ip.i << "\t" << ip.j
	   << "\t" << flg << endl;
      if(flg) {
	intpairs.push_back(ip);
	cn[ip.i] += 1.0;
	cn[ip.j] += 1.0;
      }
    }
  }

  for(int i = 1; i <= model_length; ++i) {
    cerr << "contact_number: " << i << "\t" << cn[i] << endl;
  }
}

void Model::initialize_expected_densities(void)
{
  std::fill(en_O.data(), en_O.data() + en_O.num_elements(), 0.0);
  std::fill(en_I.data(), en_I.data() + en_I.num_elements(), 0.0);

  std::fill(nE_O.data(), nE_O.data() + nE_O.num_elements(), 0.0);
  std::fill(nE_I.data(), nE_I.data() + nE_I.num_elements(), 0.0);
  std::fill(nE_OO.data(), nE_OO.data() + nE_OO.num_elements(), 0.0);
  
  std::fill(ens_OO.data(), ens_OO.data() + ens_OO.num_elements(), 0.0);
  std::fill(ens_OI.data(), ens_OI.data() + ens_OI.num_elements(), 0.0);
  std::fill(ens_IO.data(), ens_IO.data() + ens_IO.num_elements(), 0.0);
  std::fill(ens_II.data(), ens_II.data() + ens_II.num_elements(), 0.0);

  std::fill(enl_OO.data(), enl_OO.data() + enl_OO.num_elements(), 0.0);
  std::fill(enl_OI.data(), enl_OI.data() + enl_OI.num_elements(), 0.0);
  std::fill(enl_IO.data(), enl_IO.data() + enl_IO.num_elements(), 0.0);
  std::fill(enl_II.data(), enl_II.data() + enl_II.num_elements(), 0.0);
}

modelSeq& Model::get_PreModel(const int i)
{
  return PreModel[i];
}

void Model::add_pseudo_counts(darray2& tn_O, darray2& tn_I,
			      darray4& tns_OO, darray4& tns_OI,
			      darray4& tns_IO, darray4& tns_II,
			      darray4& tnl_OO, darray4& tnl_OI,
			      darray4& tnl_IO, darray4& tnl_II,
			      const bool full_nonbonded)
{
  const double C = 1.0/(GAMMA + 1.0);
  double lambda;
  for(int i = 0; i <= model_length; ++i) {
    for(int a = 0; a < nstate; ++a) {
      tn_I[i][a] = C*(tn_I[i][a] + lambdaH_I);
      for(int b = 0; b < nstate; ++b) {
	if(i == model_length && b != terminal_state) continue;
	lambda = (i == model_length) ? (0.5*GAMMA/nstate) : lambdaJ_IO;
	tns_IO[i][a][b] = C*(tns_IO[i][a][b] + lambda);
      }
      for(int b = 0; b < nstate; ++b) {
	tns_II[i][a][b] = C*(tns_II[i][a][b] + lambdaJ_II);
      }
    }
    for(int a = 0; a < nstate; ++a) {
      if(i == 0 && a != terminal_state) continue;
      tn_O[i][a] = C*(tn_O[i][a] + lambdaH_O);
      
      for(int b = 0; b < nstate; ++b) {
	if(i == model_length && b != terminal_state) continue;
	lambda = ((i == 0 || i == model_length)
		  ? (0.5*GAMMA/nstate) : lambdaJ_OO);
	tns_OO[i][a][b] = C*(tns_OO[i][a][b] + lambda);
      }
      for(int b = 0; b < nstate; ++b) {
	lambda = ((i == 0) ? (0.5*GAMMA/nstate) : lambdaJ_OI);
	tns_OI[i][a][b] = C*(tns_OI[i][a][b] + lambda);
      }
    }

    if(!full_nonbonded) continue;
    
    for(int j = i; j <= model_length; ++j) {
      for(int a = 0; a < nstate; ++a) {
	for(int b = 0; b < nstate; ++b) {
	  tnl_II[i][j][a][b] = C*(tnl_II[i][j][a][b] + lambdaK_II);
	  tnl_II[j][i][b][a] = tnl_II[i][j][a][b];
	}
      }
    }
    if(i>0) {
      for(int a = 0; a < nstate; ++a) {
	tnl_OO[i][i][a][a] = C*(tnl_OO[i][i][a][a] + GAMMA/nstate);
	for(int j = i+1; j <= model_length; ++j) {
	  for(int b = 0; b < nstate; ++b) {
	    tnl_OO[i][j][a][b] = C*(tnl_OO[i][j][a][b] + lambdaK_OO);
	    tnl_OO[j][i][b][a] = tnl_OO[i][j][a][b];
	  }
	}
      }

      for(int j = 0; j <= model_length; ++j) {
      	for(int a = 0; a < nstate; ++a) {
	  for(int b = 0; b < nstate; ++b) {
	    tnl_OI[i][j][a][b] = C*(tnl_OI[i][j][a][b] + lambdaK_OI);
	    tnl_IO[j][i][b][a] = tnl_OI[i][j][a][b];
	  }
	}
      }
    }
  }

  if(!full_nonbonded) {
    for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
      const int i = it->i, j = it->j;
      for(int a = 0; a < nstate; ++a) {
	for(int b = 0; b < nstate; ++b) {
	  tnl_OO[i][j][a][b] = C*(tnl_OO[i][j][a][b] + lambdaK_OO);
	  tnl_OO[j][i][b][a] = tnl_OO[i][j][a][b];
	}
      }
    }
  }

  return;
}

void Model::normalize_expected_densities(const double nsamples,
					 const bool full_nonbonded)
{
  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate; ++a) {
      en_I[i][a] /= nsamples;
      nE_I[i][a] /= nsamples;
      for(int b = 0; b < nstate; ++b) {
	if(i == model_length && b != terminal_state) continue;
	ens_IO[i][a][b] /= nsamples;
      }
      for(int b = 0; b < nstate; ++b) {
	ens_II[i][a][b] /= nsamples;
      }
    }
    for(int a = 0; a < nstate; ++a) {
      en_O[i][a] /= nsamples;
      nE_O[i][a] /= nsamples;
      for(int b = 0; b < nstate; ++b) {
	ens_OO[i][a][b] /= nsamples;
      }
      for(int b = 0; b < nstate; ++b) {
	ens_OI[i][a][b] /= nsamples;
      }
    }

    if(!full_nonbonded) continue;
    
    for(int j = i; j <= model_length; ++j) {
      for(int a = 0; a < nstate; ++a) {
	for(int b = 0; b < nstate; ++b) {
	  enl_II[i][j][a][b] /= nsamples;
	}
      }
    }
    for(int a = 0; a < nstate; ++a) {
      enl_OO[i][i][a][a] /= nsamples;
      for(int j = i+1; j <= model_length; ++j) {
	for(int b = 0; b < nstate; ++b) {
	  enl_OO[i][j][a][b] /= nsamples;
	  nE_OO[i][j][a][b] /= nsamples;
	}
      }
    }

    for(int j = 1; j <= model_length; ++j) {
      for(int a = 0; a < nstate; ++a) {
	for(int b = 0; b < nstate; ++b) {
	  enl_OI[i][j][a][b] /= nsamples;
	}
      }
    }
  }

  if(!full_nonbonded) {
    for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
      const int i = it->i, j = it->j;
      for(int a = 0; a < nstate; ++a) {
	for(int b = 0; b < nstate; ++b) {
	  enl_OO[i][j][a][b] /= nsamples;
	  nE_OO[i][j][a][b] /= nsamples;
	}
      }
    }
  }
  
  en_O[0][terminal_state] = 1.0;

  /***** Should this be done???? *****/
  /*
  add_pseudo_counts(en_O, en_I, ens_OO, ens_OI, ens_IO, ens_II,
  		    enl_OO, enl_OI, enl_IO, enl_II, full_nonbonded);
  */
  return;
}

void Model::init_KOO_corr(void)
{
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    const int i = it->i, j = it->j;
    for (int a = 0; a < nstate; ++a) {
      const double na = n_O[i][a];
      for (int b = 0; b < nstate; ++b) {
	const double n00 = na*n_O[j][b];
	const double nOO =  nl_OO[i][j][a][b];
	if(nOO > n00 && nOO > lambdaK_OO) {
	  K_OO[j][i][b][a] = K_OO[i][j][a][b] = 1.0;
	}
      }
    }
  }
  return;
}


