#ifndef _ANATRAJ_H_
#define _ANATRAJ_H_

#include <algorithm>
#include "model.hpp"
#include "MCseq.hpp"

struct AnaTraj {

  string ftraj;
  ifstream fin;
  int vec_len;
  int Neigen;
  darray1 eigen_vals;
  darray2 eigen_vecs;
  AnaTraj(Model& model);
  
  inline void set_trajectory_file(const string filename) {
    ftraj = filename;
    fin.open(filename.c_str());
  };
  inline void rewind() { fin.clear(); fin.seekg(0,fin.beg); };
  inline void unset_trajectory_file() { fin.close(); };
  
  bool pick_seq(Model& model, MCseq& mcs);
  void read_eigen_components(const string filename);
  darray1 pca_project(ModSeq& mseq);
  
};

vector<int> PCA_project_seqs(AnaTraj& atraj, ostream& fo, const int nseqs,
			     vector<double>& weights);
#endif // _ANATRAJ_H_
