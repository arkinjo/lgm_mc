
#include "model.hpp"

using namespace std;

void Model::dump_premodel(const Msa& msa, const string& outfile)
{
  ofstream ofs(outfile.c_str());
  
  for(int m = 0; m < num_seq; ++m) {
    ofs << msa.header(m) << endl;
    for(int i = 0; i <= model_length+1; ++i) {
      ofs << "Site\t" << i << "\t" << PreModel[m][i].O
	  << "\t";
      if(PreModel[m][i].I.length() == 0) {
	ofs << "-" << endl;
      }
      else {
	ofs << PreModel[m][i].I << endl;
      }
    }
    ofs << "//" << endl;
  }

  ofs.close();
  return;
}

void Model::read_premodel1(const string& filename)
{
  ifstream ifs(filename.c_str());
  string header;
  int isite;
  char aa, bb;
  string insert;
  modelSeq ms(model_length+2);

  num_seq = 1;
  
  for(int i = 0; i <= model_length + 1; ++i) {
    ms[i].O = '-';
    ms[i].I = "";
  }
  
  ifs >> header;

  for(int i = 0; i <= model_length + 1; ++i) {
    ifs >> header >> isite >> aa >> insert;
    ms[i].O = aa;
    if (insert.compare("-") != 0) {
      ms[i].I = insert;
    }
  }

  PreModel.resize(1);
  PreModel[0].resize(model_length+2);
  PreModel[0] = ms;

  return;
}

void Model::write_parameters(const string filename, const bool long_range_full)
{
  ofstream ofs(filename.c_str());

  ofs.setf(std::ios::scientific);
  ofs.precision(15);
  
  ofs << "ModelLength\t" << model_length << endl;
  ofs << "NumberOfSequences,SeqLength\t" << num_seq << "\t" << SeqLength << endl;
  ofs << "gamma,NS_CUTOFF,NL_CUTOFF\t" << GAMMA
      << "\t" << NS_CUTOFF
      << "\t" << NL_CUTOFF
      << endl;
  set_lambda();
  
  ofs << "#SingleSite:observed_density:expected_density:mu." << endl;
  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      ofs << "O\t" << i << "\t" << amino1_core[a] << "\t"
	  << n_O[i][a] << "\t" << en_O[i][a] << "\t"  << H_O[i][a] << endl;
    }
    if(i == model_length) continue;
    for(int a = 0; a < nstate_insert; ++a) {
      ofs << "I\t" << i << "\t" << amino1_insert[a] << "\t"
	  << n_I[i][a] << "\t" << en_I[i][a] << "\t" << H_I[i][a] << endl;
    }
  }

  ofs << "#Short-RangePair:observed_density:expected_density:J." << endl;
  for(int i = 1; i < model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ofs << "OI_s\t" << i << "\t" 
	    << amino1_core[a] << "\t" << amino1_insert[b] << "\t"
	    << ns_OI[i][a][b] << "\t" << ens_OI[i][a][b] << "\t"
	    << J_OI[i][a][b] << "\t"
	    << ns_OI[i][a][b] - n_O[i][a]*n_Int[i][b] - cs0_OI << endl;
      }

      for(int b = 0; b < nstate_core; ++b) {
	  ofs << "OO_s\t" << i << "\t"
	      << amino1_core[a] << "\t" << amino1_core[b] << "\t"
	      << ns_OO[i][a][b] << "\t" << ens_OO[i][a][b] << "\t"
	      << J_OO[i][a][b] << "\t"
	      << ns_OO[i][a][b] - n_O[i][a]*n_O[i+1][b] - cs0_OO << endl;
      }
    }

    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ofs << "II_s\t" << i << "\t" 
	    << amino1_insert[a] << "\t" << amino1_insert[b] << "\t"
	    << ns_II[i][a][b] << "\t" << ens_II[i][a][b] << "\t"
	    << J_II[i][a][b] << "\t"
	    << ns_II[i][a][b] - n_I[i][a]*n_I[i][b] - cs0_II << endl;
      }
      for(int b = 0; b < nstate_core; ++b) {
	ofs << "IO_s\t" << i << "\t" 
	    << amino1_insert[a] << "\t" << amino1_core[b] << "\t"
	    << ns_IO[i][a][b] << "\t" << ens_IO[i][a][b] << "\t"
	    << J_IO[i][a][b] << "\t"
	    << ns_IO[i][a][b] - n_Ict[i][a]*n_O[i+1][b] - cs0_IO << endl;
      }
    }
  }
 
  ofs << "#Long-RangePair:observed_density:expected_density:K:contact. "
      << intpairs.size() << endl;
  if(!long_range_full) {
    for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
      const int i = it->i, j = it->j;
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  ofs << "OO_l\t" << i << "\t" << j << "\t" 
	      << amino1_core[a] << "\t" << amino1_core[b] << "\t"
	      << nl_OO[i][j][a][b] << "\t"
	      << enl_OO[i][j][a][b] << "\t"
	      << K_OO[i][j][a][b] << "\t"
	      << nl_OO[i][j][a][b] - n_O[i][a]*n_O[j][b] - cl0_OO << endl;
	}
      }
    }
  } else {
    for(int i = 1; i <= model_length; ++i) {
      for(int j = i; j <= model_length; ++j) {
	for(int a = 0; a < nstate_core; ++a) {
	  for(int b = 0; b < nstate_core; ++b) {
	    ofs << "OO_l\t" << i << "\t" << j << "\t" 
		<< amino1_core[a] << "\t" << amino1_core[b] << "\t"
		<< nl_OO[i][j][a][b] << "\t"
		<< enl_OO[i][j][a][b] << "\t"
		<< K_OO[i][j][a][b] << "\t"
		<< nl_OO[i][j][a][b] - n_O[i][a]*n_O[j][b] - cl0_OO << endl;
	  }
	}
      }
      for(int j = 1; j < model_length; ++j) {
	for(int a = 0; a < nstate_core; ++a) {
	  for(int b = 0; b < nstate_insert; ++b) {
	    ofs << "OI_l\t" << i << "\t" << j << "\t" 
		<< amino1_core[a] << "\t" << amino1_insert[b] << "\t"
		<< nl_OI[i][j][a][b] << "\t"
		<< enl_OI[i][j][a][b] << "\t"
		<< 0 << "\t"
		<< nl_OI[i][j][a][b] - n_O[i][a]*n_I[j][b] << endl;
	  }
	}
      }
      if(i == model_length) continue;
      for(int j = i; j < model_length; ++j) {
	for(int a = 0; a < nstate_insert; ++a) {
	  for(int b = 0; b < nstate_insert; ++b) {
	    ofs << "II_l\t" << i << "\t" << j << "\t" 
		<< amino1_insert[a] << "\t" << amino1_insert[b] << "\t"
		<< nl_II[i][j][a][b] << "\t"
		<< enl_II[i][j][a][b] << "\t"
		<< 0 << "\t"
		<< nl_II[i][j][a][b] - n_I[i][a]*n_I[j][b] << endl;
	  }
	}
      }
    }
  }

  ofs.close();
  return;
}

void Model::read_parameters(const string filename)
{
  cerr << "# Reading parameters from: " << filename << endl;
  ifstream ifs(filename.c_str());
  string header;
  int ii,jj;
  char aa, bb;
  ifs >> header >> model_length;
  ifs >> header >> num_seq >> SeqLength;
  ifs >> header >> GAMMA >> NS_CUTOFF >> NL_CUTOFF;
  set_lambda();
  allocate(model_length);

  ifs >> header;
  cerr << "#processing " << header << endl;
  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      ifs >> header >> ii >> aa >> n_O[i][a] >> en_O[i][a] >>  H_O[i][a];
      n0_O[i][a] = en_O[i][a]; // reference density
    }
    if(i == model_length) continue;
    for(int a = 0; a < nstate_insert; ++a) {
      ifs >> header >> ii >> aa >> n_I[i][a] >> en_I[i][a] >> H_I[i][a];
      n0_I[i][a] = en_I[i][a]; // reference density
    }
  }

  // set len_I, n_Int and n_Ict.

  
  ifs >> header;
  cerr << "#processing " << header << endl;
  double corr;
  for(int i = 1; i < model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ifs >> header >> ii >> aa >> bb
	    >> ns_OI[i][a][b] >> ens_OI[i][a][b] >> J_OI[i][a][b] >> corr;
	ns0_OI[i][a][b] = ens_OI[i][a][b];
      }
      for(int b = 0; b < nstate_core; ++b) {
	ifs >> header >> ii >> aa >> bb 
	    >> ns_OO[i][a][b] >> ens_OO[i][a][b] >> J_OO[i][a][b] >> corr;
	ns0_OO[i][a][b] = ens_OO[i][a][b];
      }
    }

    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ifs >> header >> ii >> aa >> bb
	    >> ns_II[i][a][b] >> ens_II[i][a][b] >> J_II[i][a][b] >> corr;
	ns0_II[i][a][b] = ens_II[i][a][b];
      }
      for(int b = 0; b < nstate_core; ++b) {
	ifs >> header >> ii >> aa >> bb
	    >> ns_IO[i][a][b] >> ens_IO[i][a][b] >> J_IO[i][a][b] >> corr;
	ns0_IO[i][a][b] = ens_IO[i][a][b];
      }
    }
  }

  int nint;
  ifs >> header >> nint;
  cerr << "#processing " << header << "\t" << nint << endl;
  intpairs.clear();
  for(int k = 0; k < nint; ++k) {
    int i,j;
    double x,y,z;
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	ifs >> header >> i >> j >> aa >> bb >> x >> y >> z >> corr;
	nl_OO[i][j][a][b] = nl_OO[j][i][b][a] = x;
	enl_OO[i][j][a][b] = enl_OO[j][i][b][a] = y;
	K_OO[i][j][a][b] = K_OO[j][i][b][a] = z;
	nl0_OO[i][j][a][b] = nl0_OO[j][i][b][a] = y;
	if(a == 0 && b == 0) {
	  IPair ip(i,j);
	  intpairs.push_back(ip);
	}
      }
    }
  }

  ifs.close();
  SeqLength = get_seqlength();
  sum_I();
  return;
}

void Model::write_expected_densities(const string filename)
{
  ofstream ofs(filename.c_str());
  ofs.setf(std::ios::scientific);

  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      ofs << "O\t" << i << "\t" << amino1_core[a] << "\t"
	  << en_O[i][a] << "\t"
	  << en_O[i][a] - n_O[i][a] << "\t" 
	  << H_O[i][a] << "\t"
	  << endl;
    }
    if(i == model_length) continue;
    for(int a = 0; a < nstate_insert; ++a) {
      ofs << "I\t" << i << "\t" << amino1_insert[a] << "\t"
	  << en_I[i][a]  << "\t"
	  << en_I[i][a] - n_I[i][a] << "\t"
	  << H_I[i][a] << "\t"
	  << endl;
    }
  }

  ofs << "#Short-RangePair:observed_density:J." << endl;
  for(int i = 1; i < model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	ofs << "OO_s\t" << i << "\t" 
	    << amino1_core[a] << "\t" << amino1_core[b] << "\t"
	    << ens_OO[i][a][b] << "\t"
	    << ens_OO[i][a][b] - ns_OO[i][a][b] << endl;
      }
    }

    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ofs << "OI_s\t" << i << "\t" 
	    << amino1_core[a] << "\t" << amino1_insert[b] << "\t"
	    << ens_OI[i][a][b] << "\t"
	    << ens_OI[i][a][b] - ns_OI[i][a][b] << endl;
      }
    }

    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	ofs << "II_s\t" << i << "\t" 
	    << amino1_insert[a] << "\t" << amino1_insert[b] << "\t"
	    << ens_II[i][a][b] << "\t"
	    << ens_II[i][a][b] - ns_II[i][a][b] << endl;
      }
    }

    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	ofs << "IO_s\t" << i << "\t" 
	    << amino1_insert[a] << "\t" << amino1_core[b] << "\t"
	    << ens_IO[i][a][b] << "\t"
	    << ens_IO[i][a][b] - ns_IO[i][a][b] << endl;
      }
    }
  }

  ofs << "#Entropy,Divergence:i:ent_O:div_O:dev_O:dev_I." << endl;
  for(int i = 0; i <= model_length; ++i) {
    ofs << "Divergence" << "\t"
	<< i << "\t"
      	<< ent_O[i] << "\t"
	<< div_O[i] << "\t"
      	<< dev_O[i] << "\t"
	<< dev_I[i] << "\t"
	<< endl;
  }

  ofs.close();
  return;
}

void Model::write_covariance(const string filename)
{
  ofstream ofs(filename.c_str());
  ofs.setf(std::ios::scientific);

  for(int i = 1; i <= model_length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      for(int j = 0; j <= model_length; ++j) {
	for(int b = 0; b < nstate_core; ++b) {
	  ofs << "OO" << "\t" << i << "\t" << j << "\t"
	      << amino1_core[a] << "\t" << amino1_core[b] << "\t"
	      << nl_OO[i][j][a][b] - n_O[i][a]*n_O[j][b] << "\t"
	      << enl_OO[i][j][a][b] - en_O[i][a]*en_O[j][b] << endl;
	}
	if(j == model_length) continue;
	for(int b = 0; b < nstate_insert; ++b) {
	  ofs << "OI" << "\t" << i << "\t" << j << "\t"
	      << amino1_core[a] << "\t" << amino1_insert[b] << "\t"
	      << nl_OI[i][j][a][b] - n_O[i][a]*n_I[j][b] << "\t"
	      << enl_OI[i][j][a][b] - en_O[i][a]*en_I[j][b] << endl;
	}
      }
    }
    if(i == model_length) continue;
    for(int a = 0; a < nstate_insert; ++a) {
      for(int j = 1; j <= model_length; ++j) {
	if(j != model_length) {
	  for(int b = 0; b < nstate_insert; ++b) {
	    ofs << "II" << "\t" << i << "\t" << j << "\t"
		<< amino1_insert[a] << "\t" << amino1_insert[b] << "\t"
		<< nl_II[i][j][a][b] - n_I[i][a]*n_I[j][b] << "\t"
		<< enl_II[i][j][a][b] - en_I[i][a]*en_I[j][b] << endl;
	  }
	}
	for(int b = 0; b < nstate_core; ++b) {
	  ofs << "IO" << "\t" << j << "\t" << i << "\t"
	      << amino1_insert[a] << "\t" << amino1_core[b] << "\t"
	      << nl_OI[j][i][b][a] - n_I[i][a]*n_O[j][b] << "\t"
	      << enl_OI[j][i][b][a] - en_I[i][a]*en_O[j][b] << endl;
	}
      }
    }
  }
}

void Model::dump_covariance(darray2& tn_O, darray2& tn_I,
			    darray4& tnl_OO, darray4& tnl_OI,
			    darray4& tnl_IO, darray4& tnl_II,
			    const string filename)
{
  ofstream ofs(filename.c_str());
  ofs.setf(std::ios::scientific);

  for(int i = 0; i <= model_length; ++i) {
    if(i > 0) {
      for(int a = 0; a < nstate_core; ++a) {
	for(int j = 0; j <= model_length; ++j) {
	  if(j > 0) {
	    for(int b = 0; b < nstate_core; ++b) {
	      ofs << tnl_OO[i][j][a][b] - tn_O[i][a]*tn_O[j][b] << "\t";
	    }
	  }
	  for(int b = 0; b < nstate_insert; ++b) {
	    ofs << tnl_OI[i][j][a][b] - tn_O[i][a]*tn_I[j][b] << "\t";
	  }
	}
	ofs << endl;
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      for(int j = 0; j <= model_length; ++j) {
	if(j > 0) {
	  for(int b = 0; b < nstate_core; ++b) {
	    ofs << tnl_IO[i][j][a][b] - tn_I[i][a]*tn_O[j][b] << "\t";
	  }
	}
	for(int b = 0; b < nstate_insert; ++b) {
	  ofs << tnl_II[i][j][a][b] - tn_I[i][a]*tn_I[j][b] << "\t";
	}
      }
      ofs << endl;
    }
  }
}
