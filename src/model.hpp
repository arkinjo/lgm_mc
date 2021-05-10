#ifndef MODEL_H_
#define MODEL_H_

#include "base.hpp"
#include "msa.hpp"

using namespace std;

const double DEFAULT_GAMMA = 0.1;

struct SSite 
{
  char O; // matching or deletion
  string I; // insertion
  SSite();
};

typedef vector<SSite> modelSeq;

struct IPair
{
  int i;
  int j;
  IPair(int i, int j);
};

typedef vector<IPair> IntPairs; // list of interacting pairs of core sites.

struct Model 
{
  static const bool core_only = true;
  int model_length;
  int msa_length;
  int num_seq;
  double temperature,beta;
  vector<bool> msa_sites;
  boost::mt19937 rng;
  boost::normal_distribution<> nd;

  void assign_sites_1st(const Msa&);
  void assign_sites(const Msa&);
  void assign_sites_hmm(const Msa&);
  void assign_sites_blast(const Msa&);
  vector<double> pinsert;
  modelSeq seq_to_modelSeq(const string &seq);
  modelSeq seq_to_modelSeq_blast(const Msa&, const int);
  
  vector<modelSeq> PreModel;
  vector<double> weights;
  double SeqLength;
  IntPairs intpairs_ref; // interacting pairs of residues in the reference sequence/structure.
  map<int,int> intref2mod; // mapping from residue positions in the reference to core sites.
  double GAMMA;  

  Model();
  Model(const Msa&, const bool hmm=true, const bool blast=false);
  ~Model();

  modelSeq &get_PreModel(const int i);
  inline int get_num_seq(void) { return num_seq; };
  void set_weights(void); //sets weights to sequences in MSA (weights[num_seq])
  inline double get_weight(const int i) const  { return weights[i]; };
  void set_temperature(const double t);
  inline double get_temperature(void) { return temperature; };
  void stat_site(const vector<double>& weights);
  void stat_site_bonded(const vector<double>& weights);
  void stat_site_nonbonded(const vector<double>& weights);
  void stat_all(void);

  void add_pseudo_counts(darray2& tnO, darray2& tnI,
			 darray4& tnsOO, darray4& tnsOI,
			 darray4& tnsIO, darray4& tnsII,
			 darray4& tnlOO, darray4& tnlOI,
			 darray4& tnlIO, darray4& tnlII,
			 const bool full_nonbonded);
  void allocate(const int mlen);
  void dump_premodel(const Msa& msa, const string& outfile);
  void read_premodel1(const string& filename);

  inline int length() const { return model_length; };
  inline IntPairs& get_intpairs() { return intpairs; };
  IntPairs intpairs; // interacting pairs in terms of core sites.

  double get_seqlength_generic(darray2& no, darray2& ni,
			       double& lO, double& lI);
  double get_seqlength(double& lO, double& lI);
  double get_seqlength(void);
  double get_seqlength_ref(double& lO, double& lI);
  double get_seqlength_ref(void);
  double get_seqlength_exp(double& lO, double& lI);
  double get_seqlength_exp(void);
  
  // regularization parameters

  double lambdaH_O,lambdaH_I;
  double lambdaJ_OO, lambdaJ_II, lambdaJ_OI, lambdaJ_IO;
  double lambdaK_OO, lambdaK_OI, lambdaK_II;
  //default correlation coefficient due to pseudo-counts.
  double cs0_OO,cs0_OI,cs0_IO,cs0_II;
  double cl0_OO,cl0_OI,cl0_IO,cl0_II;
  
  void set_lambda(void);
  inline double get_gamma(void) { return GAMMA;};
  inline void set_gamma(const double g) {
    GAMMA = g;
    set_lambda();
  };

  double NS_CUTOFF;
  double NL_CUTOFF;
  void set_nl_cutoff(const int ninc);
  //  darray1 ilength;
  // observed densities
  darray1 site_energy, site_energy_exp, site_energy_ref;
  darray1 len_I;
  darray2 n_O, n_I;  // site number densities
  darray2 n_Int, n_Ict;  // N- and C-terminii of each insert.
  darray4 ns_OO, ns_OI, ns_IO, ns_II;  // short-range pair number densities 
  darray4 nl_OO,nl_OI,nl_IO,nl_II;  // long-range pair number densities


  // reference densities
  darray2 n0_O, n0_I;  // reference site number densities
  darray4 ns0_OO, ns0_OI, ns0_IO, ns0_II;  // short-range pair number densities 
  darray4 nl0_OO,nl0_OI,nl0_IO,nl0_II;  // long-range pair number densities
  
  // expected densities
  darray2 en_O, en_I;  // site number densities
  darray4 ens_OO, ens_OI, ens_IO, ens_II;  // short-range pair number densities 
  darray4 enl_OO,enl_OI, enl_IO, enl_II;  // long-range pair number densities

  // correlation between single-site density and energy
  darray2 nE_O, nE_I;

  // correlation between non-bonded pair density and energy
  darray4 nE_OO;
  
  void sum_I(void);
  void init_params(const double ep=-2.0);
  void init_KOO_corr(void);
  void initialize_expected_densities(void);
  void normalize_expected_densities(const double nsamples,
				    const bool full_nonbonded=false);
  
  double Mu_len;
  void set_mu_len(const double mu);
  void set_mu_lenI(const double mu);
  darray2 H_O, H_I; // "chemical potential"
  darray4 J_OO, J_OI, J_IO, J_II; // short-range coupling
  darray4 K_OO; // long-range coupling

  darray1 div_O,div_R; // KL-divergence between en_S and n_S.
  darray1 ddiv_O, mdiv_O; // d(div_O)/dT , d(diV_O)/dmu

  darray1 ent_O; // KL-divergence between en_S and n_S.
  darray1 dev_O,dev_I; // deviation between en_S and n_S.
  darray2 div_OO;
  darray1 div_OOs, ddiv_OOs; // sum_j d(div_OiOj)/dT

  // Gauge fixing
  void set_Ising_gauge_H(void);
  void set_Ising_gauge_J(void);
  void set_Ising_gauge_K(void);
  inline void set_Ising_gauge(void) {
    set_Ising_gauge_K();
    set_Ising_gauge_J();
    set_Ising_gauge_H();
  }
  double zero_Jp(const int i);
  double zero_Jm(const int i);

  void eliminate_H(void);
  void set_JOp(const int i, darray1& jv);
  void set_JIp(const int i, darray1& jv);
  void set_JOm(const int i, darray1& jv);
  void set_JIm(const int i, darray1& jv);
  double symmetrize_JI(const int i);
  double symmetrize_JO(const int i);
  double normalize_JOp(const int i);
  double normalize_JOm(const int i);
  double update_J_gauge(void);
  double fix_B_gauge_J(void);
  double fix_B_gauge_K(void);
  void fix_B_gauge(void);

  // "Sequence" computation
  double get_divergence(darray2& tn_O, darray1& tdiv);
  double get_divergence_from_ref(void);
  double get_divergence(void);
  double get_deviation(double& dev_O, double& dev_I);
  double get_site_entropy_generic(darray2& tn_O);
  double get_site_entropy(void);
  double get_site_entropy_ref(void);
  double get_site_entropy_exp(void);
  
  double get_energy_generic(darray2& tnO, darray2& tnI,
			    darray4& tns_OO, darray4& tns_OI,
			    darray4& tns_IO, darray4& tns_II,
			    darray4& tnl_OO,
			    double &eshort, double &elong, double &egibbs,
			    darray1& site_energy);

  double get_energy(double &eshort, double &elong, double &egibbs);
  double get_energy(void);
  double get_energy_ref(double &eshort, double &elong, double &egibbs);
  double get_energy_ref(void);
  double get_energy_exp(double &eshort, double &elong, double &egibbs);
  double get_energy_exp(void);
  double chem_pot_work(void);
  double get_gibbs_energy(void);
  
  // IO stuffs
  void read_ipairs(const string filename);
  void read_reference_3Dalignment(const string filename);
  void align_ipairs_1st(Msa& msa);
  void set_intpairs(void);
  void set_default_ipairs(void);
  
  void write_parameters(const string filename,
			const bool long_range_full = false);
  void read_parameters(const string filename);
  void write_covariance(const string filename);
  void dump_covariance(darray2& tn_O, darray2& tn_I,
		       darray4& tnl_OO, darray4& tnl_OI,
		       darray4& tnl_IO, darray4& tnl_II,
		       const string filename);
  inline void dump_covariance_obs(const string filename) {
    dump_covariance(n_O, n_I, nl_OO, nl_OI, nl_IO, nl_II, filename);
  };
  inline void dump_covariance_exp(const string filename) {
    dump_covariance(en_O, en_I, enl_OO, enl_OI, enl_IO, enl_II, filename);
  };
  
  void write_expected_densities(const string filename);
  
  void check_const_site(void);
  void check_const_site_exp(void);
  void check_const_bonded(void);
  void check_const_bonded_exp(void);
  void check_const_nonbonded(void);
  void check_const_nonbonded_exp(void);

  double check_site(const bool show);
  double check_bonded(const bool show);
  double check_nonbonded(const bool show);

  void set_null_sequence(void);//for testing
  void set_uniform_sequence(void);//for testing

};

#endif // MODEL_H_
