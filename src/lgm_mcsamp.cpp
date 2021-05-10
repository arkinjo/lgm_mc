#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>

#include "msa.hpp"
#include "model.hpp"
#include "MCseq.hpp"
#include "OptParam.hpp"
using namespace std;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h","produce help message")
    ("param,p", po::value<string>(), "LGM parameter file")
    ("logbase,l", po::value<string>(), "basename for log files")
    ("traj,t", po::value<bool>(), "save trajectory in 'LOGBASE.traj' or not")
    ("param_out", po::value<bool>(), "save parameters in 'LOGBASE.param_out' or not")
    ("gibbs", po::value<bool>(), "use Gibbs sampler or not")
    ("nsweep", po::value<int>(), "number of sweeps")
    ("nsweep_eq", po::value<int>(), "number of sweeps for equilibration")
    ("seed", po::value<int>(), "random seed")
    ("temp", po::value<double>(), "temperature")
    ("temp_eq", po::value<double>(), "temperature in the equilibration stage")
    ("mu_len", po::value<double>(), "mu_len")
    ("mu_I", po::value<double>(), "mu_I")
    ("pm_site", po::value<int>(), "site for point mutation")
    ("pm_res", po::value<char>(), "residue for point mutation")
    ("pm_mu", po::value<double>(), "mu_O(a) for  point mutation")
    ("nmon", po::value<int>(), "monitoring interval")
    ("ntraj", po::value<int>(), "interval for saving trajectory")
    ;

  po::positional_options_description p;
  po::variables_map vm;
  po::store(po::command_line_parser(argc,argv).
	    options(desc).positional(p).run(), vm);
  po::notify(vm);

  if(vm.count("help")) {
    cout << desc << endl;
    return 1;
  }

  string logbase = "mcsamp";
  if(vm.count("logbase")) {
    logbase = vm["logbase"].as<string>();
  }
  cerr << "# LOGBASE= " << logbase << endl;

  if(!vm.count("param")) {
    cerr << "Specify the parameter file!!" << endl;
    exit(1);
  }

  Model model;
  model.read_parameters(vm["param"].as<string>());

  double temp_eq = 1;
  if(vm.count("temp_eq")) {
    temp_eq  = vm["temp_eq"].as<double>();

  }
  cerr << "temperature(eq.) = " << temp_eq << endl;
  model.set_temperature(temp_eq);

  double Mu_len=0;
  if(vm.count("mu_len")) {
    Mu_len  = vm["mu_len"].as<double>();
  }
  cerr << "Mu_len = " << Mu_len << endl;
  model.set_mu_len(Mu_len);
  double Mu_I = 0;
  if(vm.count("mu_I")) {
    Mu_I  = vm["mu_I"].as<double>();
  }
  cerr << "Mu_I = " << Mu_I << endl;
  model.set_mu_lenI(Mu_I);
  
  int pm_site = -1;
  char pm_res = 'A';
  double pm_mu = 10.0;
  if(vm.count("pm_site")) {
    pm_site  = vm["pm_site"].as<int>();
  }
  if(vm.count("pm_res")) {
    pm_res  = vm["pm_res"].as<char>();
  }
  if(vm.count("pm_mu")) {
    pm_mu  = vm["pm_mu"].as<double>();
  }
  if(pm_site > 0) {
    cerr << "point_mutation(site,res,mu): " << pm_site << "\t" << pm_res << "\t" << pm_mu
	 << endl;
    model.H_O[pm_site][find_amino1_core(pm_res)] += pm_mu;
  }
  
  MCseq mcs(model);
  if(vm.count("seed")) {
    const int seed = vm["seed"].as<int>();
    cerr << "random seed = " << seed << endl;
    mcs.rng.seed(seed);
  }

  mcs.select_type = CANONICAL;

  bool gibbs = false;
  if(vm.count("gibbs")) {
    gibbs = vm["gibbs"].as<bool>();
  }
  mcs.gibbs_sampler_p = gibbs;
  
  int nsweeps_eq = 1000;
  if(vm.count("nsweep_eq")) {
    nsweeps_eq = vm["nsweep_eq"].as<int>();
  }
    
  int nsweeps = 100000;
  if(vm.count("nsweep")) {
    nsweeps = vm["nsweep"].as<int>();
  }
    
  mcs.NMONITOR = 100; //nsweeps / 10 ;
  if(vm.count("nmon")) {
    mcs.NMONITOR = vm["nmon"].as<int>();
  }
  if(vm.count("ntraj")) {
    mcs.NTRAJ = vm["ntraj"].as<int>();
  }
  cerr << "Saving trajectory for every " << mcs.NTRAJ
       << " sweeps" << endl;

  double temp = 1;
  if(vm.count("temp")) {
    temp  = vm["temp"].as<double>();

  }
  cerr << "temperature = " << temp << endl;

  
  mcs.set_random_sequence(model,0);
  //  mcs.gibbs_sampler_p = true;
  if(nsweeps_eq > 0) {
    cerr << "# Annealing for " << nsweeps_eq << " sweeps" << endl;
    mcs.MONITOR_p = false;
    mcs.get_energy(model);
    //    mcs.anneal(model, nsweeps_eq/10, 10, 1.5, temp);
    mcs.set_temperature(temp);
    mcs.run(model, nsweeps_eq, false);
  }
  mcs.set_temperature(temp);
  //  mcs.gibbs_sampler_p = gibbs;
  
  if(vm.count("traj")) {
    const string ftraj = logbase + ".traj";
    cerr << "# Trajectory saved to " <<  ftraj << endl;
    mcs.set_trajectory_file(ftraj);
  }
  const clock_t begin_time = clock();
  cerr << "# Production run for " << nsweeps << " sweeps" << endl;
  mcs.MONITOR_p = true;
  mcs.stats_full_nonbonded = true;// ******** 
  mcs.init_stats(model);
  int nsamples = mcs.run(model, nsweeps, true);
  mcs.finish_stats(model);
  assert(nsamples == nsweeps);
  const clock_t end_time = clock();
  cerr << "Time " << float(end_time - begin_time)/CLOCKS_PER_SEC << endl;

  if(mcs.trajectory_p) {
    mcs.unset_trajectory_file();
  }
  cerr << "#  run for " << nsweeps << " sweeps" << endl;

  double slenO0,slenI0;
  const double slen0 = model.get_seqlength(slenO0,slenI0);
  double egibbs0,ebonded0,enonbonded0;
  const double energy0 = model.get_energy(ebonded0,enonbonded0,egibbs0);

  cerr << "slen0 " << slen0 << "\t"
       << slenO0 << "\t"
       << slenI0 << "\t" << endl;
  const double SpecHeat = mcs.Evar/(temp*temp);
  double slenO,slenI;
  const double slen = model.get_seqlength_exp(slenO,slenI);
  cerr << "slen " << slen << "\t"
       << slenO << "\t"
       << slenI << "\t" << endl;

  double devO,devI;
  const double dev = model.get_deviation(devO,devI);
  double egibbs,ebonded,enonbonded;
  const double energy = model.get_energy_exp(ebonded,enonbonded,egibbs);
  const double ent = model.get_site_entropy_exp();
  const double divO = model.get_divergence();
  const double divR = model.get_divergence_from_ref();
  cout << "Summary:\t"
       << nsweeps << "\t"
       << Mu_I << "\t"
       << mcs.Temperature << "\t"
       << SpecHeat << "\t"
       << slen - slen0 << "\t"
       << slenO - slenO0 << "\t"
       << slenI - slenI0 << "\t"
       << dev  << "\t"
       << devO << "\t"
       << devI << "\t"
       << energy - energy0 << "\t"
       << ebonded - ebonded0 << "\t"
       << enonbonded - enonbonded0 << "\t"
       << egibbs - egibbs0 << "\t"
       << ent << "\t"
       << divO << "\t"
       << divR << endl;


  for(int i = 1; i <= model.length(); ++i) {
    cout << "Site\t" << i << "\t"
	 << model.dev_O[i] << "\t"
	 << model.dev_I[i] << "\t"
	 << model.ent_O[i] << "\t" 
	 << model.div_O[i] << "\t"
      	 << model.div_R[i] << "\t" 
	 << model.site_energy_exp[i] << endl;
  }

  if(vm.count("param_out")) {
    model.write_parameters(logbase + ".param_out", mcs.stats_full_nonbonded);
  }

  return 0;
}
