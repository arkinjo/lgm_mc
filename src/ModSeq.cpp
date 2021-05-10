#include "ModSeq.hpp"

void ModSeq::set_model_length(const int len)
{
  Length = len;
  for(int i = 0; i <= len + 1; ++i) {
    ISite isite;
    aaseq.push_back(isite);
  }
}

int ModSeq::get_seqlength(void)
{
  SeqLength=0;
  for(int i = 1; i <= Length; ++i) {
    if(aaseq[i].O != terminal_state) SeqLength++;
    if(i < Length) SeqLength += aaseq[i].I.size();
  }
  return SeqLength;
}

double ModSeq::get_energy(Model& model)
{
  Energy = Ebonded = Enonbonded = Egibbs = 0.0;
  SeqLength = 0;
  for(int i = 1; i <= Length; ++i) {
    Egibbs -= model.H_O[i][aaseq[i].O];
    if(aaseq[i].O != terminal_state) SeqLength++;

    if(i==Length) continue;
    const int ilen = aaseq[i].I.size();
    SeqLength += ilen;
    if(ilen > 0) {
      for(int j = 0; j < ilen; ++j) {
	Egibbs -= model.H_I[i][aaseq[i].I[j]];
      }
      Ebonded -= model.J_OI[i][aaseq[i].O][aaseq[i].I[0]]
	+ model.J_IO[i][aaseq[i].I[ilen-1]][aaseq[i+1].O];
      for(int j = 0; j < ilen-1; ++j) {
	Ebonded -= model.J_II[i][aaseq[i].I[j]][aaseq[i].I[j+1]];
      }
    } else {
      Ebonded -= model.J_OO[i][aaseq[i].O][aaseq[i+1].O];
    }
  }

  IntPairs &intpairs = model.get_intpairs();
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    int i = it->i, j = it->j;
    Enonbonded -= model.K_OO[i][j][aaseq[i].O][aaseq[j].O];
  }

  Energy = Ebonded + Enonbonded + Egibbs;
  return Energy;
}

double ModSeq::get_divergence(Model& model)
{
  double div=0.0;
  
  for(int i = 1; i <= Length; ++i) {
    div -= log(model.n_O[i][aaseq[i].O]);
  }

  return (div/Length);
}

void ModSeq::set_random_sequence(Model& model, const int ilen)
{
  boost::random::uniform_int_distribution<> rand_core(0,nstate_core-1),
    rand_insert(0,nstate_insert-1), rand_insert_len(0,ilen);

  for(int i = 1; i <= Length; ++i) {
    aaseq[i].O = rand_core(rng);
    aaseq[i].I.clear();
    if(i == Length) continue;
    while(rand_insert_len(rng)) {
      aaseq[i].I.push_back(rand_insert(rng));
    }

  }

  get_energy(model);
  return;
}

void ModSeq::set_natural_sequence(Model& model, const int imodel)
{
  modelSeq &aseq = model.get_PreModel(imodel);

  for(int i = 1; i <= Length; ++i) {
    aaseq[i].O = find_amino1_core(aseq[i].O);
    aaseq[i].I.clear();
    if(i == Length) continue;
    for(string::iterator it = aseq[i].I.begin(); it != aseq[i].I.end(); ++it) {
      aaseq[i].I.push_back(find_amino1_insert(*it));
    }
  }
  get_energy(model);
  return;
}

void ModSeq::set_natural_sequence(Model& model)
{
  int nseq = model.get_num_seq();
  boost::random::uniform_int_distribution<> rand(0,nseq - 1);
  
  const int imodel = rand(rng);

  set_natural_sequence(model, imodel);
  
  return;
}

void ModSeq::dump_aaseq(void)
{
  int len=0;
  for(int i = 1; i <= Length; ++i) {
    cerr << amino1_core[aaseq[i].O];
    if(aaseq[i].O != terminal_state) len++;
    for(int k = 0; k < aaseq[i].I.size(); ++k) {
      cerr << amino1_insert[aaseq[i].I[k]];
      len++;
    }
  }
  cerr << endl;
  cerr << "(length= " << len << ")" << endl;
}

vector<double> ModSeq::seq2vec(void)
{
  const int vlen = Length*nstate_core + (Length-1)*nstate_insert;
  vector<double> v(vlen,0.0);
  int k=0;
  for(int i = 1; i <= Length; ++i) {
    for(int a = 0; a < nstate_core; ++a) {
      v[k] = (aaseq[i].O == a) ? 1.0 : 0.0;
      k++;
    }
    if(i==Length) continue;
    vector<double> ni(nstate_insert,0.0);
    for(int l = 0; l < aaseq[i].I.size(); ++l) {
      ni[aaseq[i].I[l]]++;
    }
    for(int a = 0; a < nstate_insert; ++a) {
      v[k] = ni[a];
      k++;
    }
  }
  assert(k == vlen);
  return v;
}
