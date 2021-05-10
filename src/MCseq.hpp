#ifndef MCSEQ_H_
#define MCSEQ_H_

#include "ModSeq.hpp"
#include "model.hpp"

using namespace std;


enum SelectionType { CANONICAL, MULTI_CANONICAL };

struct MCseq: ModSeq {
  MCseq(Model& model);
  int imin_natural,imax_natural;
  double emin_natural,emax_natural;
  double eave_natural,evar_natural;
  void scan_natural_sequences(Model& model);
  
  double WeightTotal,Evar,Eave;
  double Lvar,Lave;
  inline void init_stats(Model& model) {
    WeightTotal = Eave = Evar = 0.0;
    Lave = Lvar = 0.0;
    model.initialize_expected_densities();
  };
  inline void finish_stats(Model& model) {
    model.normalize_expected_densities(WeightTotal, stats_full_nonbonded);
    Evar /= WeightTotal;
    Eave /= WeightTotal;
    Evar -= (Eave*Eave);
    Lvar /= WeightTotal;
    Lave /= WeightTotal;
    Lvar -= (Lave*Lave);
  };
  bool trajectory_p;
  bool gibbs_sampler_p;
  void save_DOS(const string filename);
  int read_DOS(const string filename, const bool reset_max);
  void save_trajectory();
  ofstream ftraj;
  void set_trajectory_file(string filename);
  inline void unset_trajectory_file(void) {
    ftraj.close();
    trajectory_p = false;
  };

  SelectionType select_type;
  bool select_criteria(const double dene, const double coef);
  bool Metropolis(const double dene, const double coef);
  bool WL(const double dene, const double coef);
  double WLprec; // WL DOS precision limit
  //// for Wang-Landau simulations
  void setup_WL(const double min_energy, const double max_energy,
		const double delta_energy=0.1);
  double DELTA_ENERGY, WLfactor, MAX_ENERGY, MIN_ENERGY, MAX_DOS, MIN0_ENERGY;
  map<int,double> EHist, DOS, EHistTraj;
  bool check_EHist(void);
  int traj_min_count(void);
  void clear_DOS(void);
  void clear_EHist(void);
  //  int get_DOS_bindex(const double e);
  inline int get_DOS_bindex(const double e)
  {
    return (int) floor(e/DELTA_ENERGY);
  };


  ////
  map< int, vector<int> > IntList; // interactions list for fast look-up.
  double Temperature, Beta;
  inline void set_temperature(const double temp) {
    Temperature = temp;
    Beta = 1.0/temp;
  };
  
  double dist_core[nstate_core], dist_insert[nstate_insert];
  double ep_core[nstate_core], ep_insert[nstate_insert];

  double mutate_core_gs(const Model& model, const int i, int& anew);
  double extend_insert_gs(const Model& model, const int i, const int j,
			  int& anew);
  double mutate_insert_gs(const Model& model, const int i, const int j,
			  int& anew);

  double mutate_core(const Model& model, const int i, const int anew);
  double mutate_insert(const Model& model, const int i, const int j, const int anew);
  double extend_insert(const Model& model, const int i, const int j, const int anew);
  double shorten_insert(const Model& model, const int i, const int j);
  
  int NMONITOR;
  bool MONITOR_p;
  int NTRAJ;
  double sweep(const Model& model);
  // Gibbs sampling
  double sweep_gs(const Model& model);

  int anneal(Model& model, const int nsweeps, const int nsteps,
	     const double temp_high, const double temp_low);
  
  int run(Model& model, const int nsweeps, const bool samplep);
  bool stats_full_nonbonded;
  bool STATS_MULTICANONICAL;
  double count_stats(Model& model, const bool just_weight=false);
};

#endif //MCSEQ_H_
