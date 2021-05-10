#include "AnaTraj.hpp"

AnaTraj::AnaTraj(Model& model)
{
  const int mlen = model.length();
  vec_len = mlen*nstate_core + (mlen-1)*nstate_insert;
}

bool AnaTraj::pick_seq(Model& model, MCseq& mcs)
{
  string line,seq;
  double energy;
  
  if(getline(fin,line)) {
    istringstream iss(line);
    iss >> energy >> seq;
    //    cerr << "AnaTraj::pick_seq " << energy << "\t" << seq << endl;
  }
  else return false;

  int i=0,length=0;
  mcs.aaseq[0].I.clear();
  for(string::iterator it = seq.begin(); it != seq.end(); ++it) {
    if(amino1_core_p(*it)) {
      i++;
      const int a = find_amino1_core(*it);
      mcs.aaseq[i].O = a;
      mcs.aaseq[i].I.clear();
      if(a != terminal_state) length++;
    } else {
      mcs.aaseq[i].I.push_back(find_amino1_insert(*it));
      length++;
    }
  }
  mcs.Energy = energy;
  //  mcs.get_energy(model);
  //cerr << "Energy(read,comp) " << energy << "\t" << mcs.Energy << endl;

  mcs.aaseq[0].O = mcs.aaseq[mcs.Length+1].O = terminal_state;

  mcs.SeqLength = length;

  return true;
}

void AnaTraj::read_eigen_components(const string filename)
{

  ifstream fi;
  fi.open(filename.c_str());
  
  fi >> Neigen;
  cerr << "read_eigen_components: " << Neigen << "\t" << vec_len << endl;
  darray1::extent_gen ext1;
  darray2::extent_gen ext2;
  eigen_vals.resize(ext1[Neigen]);
  eigen_vecs.resize(ext2[Neigen][vec_len]);

  double a,b,c;
  fi >> a >> b >> c;
  eigen_vals[0] = a;
  eigen_vals[1] = b;
  eigen_vals[2] = c;
  for(int i = 0; i < vec_len; ++i) {
    fi >> a >> b >> c;
    eigen_vecs[0][i] = a;
    eigen_vecs[1][i] = b;
    eigen_vecs[2][i] = c;
  }
  fi.close();
}

darray1 AnaTraj::pca_project(ModSeq& mseq)
{
  vector<double> v = mseq.seq2vec();
  darray1 p(boost::extents[Neigen]);
  for(int k = 0; k < Neigen; ++k) {
    p[k] = 0.0;
    const double s = sqrt(eigen_vals[k]);
    for(int i = 0; i < vec_len; ++i) {
      p[k] += s*eigen_vecs[k][i]*v[i];
    }
  }

  return p;
}

vector<int> PCA_project_seqs(AnaTraj& atraj, ostream& fo, const int nseqs, vector<double>& weights)
{
  boost::random::mt19937 rng;
  boost::random::discrete_distribution<> dist(weights);
  vector<int> inds;

  for(int ii = 0; ii != nseqs; ++ii) {
    const int i = dist(rng);
    inds.push_back(i);
  }
  std::sort(inds.begin(), inds.end());
  return inds;
}
