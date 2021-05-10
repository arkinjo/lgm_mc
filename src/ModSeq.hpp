#ifndef _MODSEQ_H_
#define _MODSEQ_H_
#include "base.hpp"
#include "model.hpp"

struct ISite {
  // c.f., SSite in model.hpp 
  int O; // core
  vector<int> I; // insert
  ISite();
};

typedef vector<ISite> ISeq;

// representation of a sequence.
struct ModSeq {
  ModSeq() {};
  ModSeq(const int len) { set_model_length(len); };
  ISeq aaseq;
  double Energy, Ebonded, Enonbonded, Egibbs;
  int Length; // model length
  int SeqLength; // Number of "real" residues in aaseq.
  boost::random::mt19937 rng;
  
  inline void set_seed(const int s) { rng.seed(s); };
  void set_model_length(const int len);
  double get_energy(Model& model);
  double get_divergence(Model& model);
  int get_seqlength(void);
  void set_random_sequence(Model& model, const int ilen=0);
  void set_natural_sequence(Model &model, const int imodel);
  void set_natural_sequence(Model &model);
  void dump_aaseq(void);
  vector<double> seq2vec(void);
};

#endif // _MODSEQ_H_
