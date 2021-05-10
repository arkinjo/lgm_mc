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
    ("dos", po::value<string>(), "WL DOS file")
    ("param,p", po::value<string>(), "LGM parameter file")
    ("logbase,l", po::value<string>(), "basename for log files")
    ("nsweep", po::value<int>(), "number of sweeps")
    ("nloop", po::value<int>(), "number of loops")
    ("mu_len", po::value<double>(), "Mu_len")
    ("mu_I", po::value<double>(), "Mu_I")
    ("traj,t", po::value<bool>(), "save trajectory in 'LOGBASE.traj' or not")
    ("seed", po::value<int>(), "random seed")
    ("ntraj", po::value<int>(), "interval for saving trajectory")
    ("temp_max", po::value<double>(), "highest temperature")
    ("min_cnt", po::value<int>(), "minimum count for each bin")
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

  string logbase = "wldos";
  if(vm.count("logbase")) {
    logbase = vm["logbase"].as<string>();
  }
  cerr << "# LOGBASE= " << logbase << endl;

  Model model;
  if(vm.count("param")) {
    model.read_parameters(vm["param"].as<string>());
  } else {
    cerr << "Specify the parameter file!" << endl;
    exit(1);
  }
  if(vm.count("mu_len")) {
    const double mu_len  = vm["mu_len"].as<double>();
    cerr << "Mu_len = " << mu_len << endl;
    model.set_mu_len(mu_len);
  }
  if(vm.count("mu_I")) {
    const double mu_I  = vm["mu_I"].as<double>();
    cerr << "Mu_I = " << mu_I << endl;
    model.set_mu_lenI(mu_I);
  }

  double tmax=1.2;
  if(vm.count("temp_max")) {
    tmax = (vm["temp_max"].as<double>());
  }

  int min_cnt = 1000;
  if(vm.count("min_cnt")) {
    min_cnt = (vm["min_cnt"].as<int>());
  }
  cerr << "Minimum count for convergence: " << min_cnt << endl;
  
  MCseq mcs(model);
  if(vm.count("seed")) {
    mcs.set_seed(vm["seed"].as<int>());
    cerr << "random seed= " << (vm["seed"].as<int>()) << endl;;
  }
    
  int nsweeps = 100000;
  if(vm.count("nsweep")) {
    nsweeps = vm["nsweep"].as<int>();
  }
  int nloops = 100;
  if(vm.count("nloop")) {
    nloops = vm["nloop"].as<int>();
  }

  mcs.NTRAJ = 1;
  if(vm.count("ntraj")) {
    mcs.NTRAJ = vm["ntraj"].as<int>();
  }
  cerr << "Saving trajectory for every " << mcs.NTRAJ
       << " sweeps" << endl;

  mcs.NMONITOR = 10000; //nsweeps / 10 ;
  mcs.setup_WL(0,0,0.1);
  if(vm.count("dos")) {
    mcs.read_DOS(vm["dos"].as<string>(), true);
    cerr << "Reading DOS file: " << vm["dos"].as<string>() << endl;;
  } else {
    cerr << "Specify the DOS file with the '--dos' option!!" << endl;
    exit(1);
  }

  if(vm.count("traj")) {
    const string ftraj = logbase + ".traj";
    cerr << "# Trajectory saved to " <<  ftraj << endl;
    mcs.set_trajectory_file(ftraj);
  }
  const clock_t begin_time = clock();
  cerr << "# WL sampling run for "
       << nsweeps << " sweeps." << endl;

  int iloop,isweeps;
  mcs.gibbs_sampler_p = false;

  mcs.init_stats(model);
  mcs.stats_full_nonbonded = false;
  for(int iloop = 0; iloop != nloops; ++iloop) {
    cerr << "#Loop: " << iloop << endl;
    mcs.trajectory_p = false;
    mcs.select_type = CANONICAL;
    mcs.set_random_sequence(model);
    mcs.set_temperature(tmax);
    mcs.run(model, 1000, false);
    mcs.set_temperature(1.0);
    while(1) {
      mcs.run(model, 10, false);
      if(mcs.Energy > mcs.MIN_ENERGY && mcs.Energy < mcs.MAX_ENERGY) break;
    }
    mcs.select_type = MULTI_CANONICAL;
    isweeps = mcs.run(model, 1000, false);
    mcs.trajectory_p = true;
    isweeps = mcs.run(model, nsweeps, true);
    const int traj_min_cnt = mcs.traj_min_count();
    if(traj_min_cnt > min_cnt) {
      cerr << "Multicanonical sampling converged in loop: " << iloop << endl;
      break;
    }
  }
  mcs.finish_stats(model);
  const clock_t end_time = clock();
  cerr << "Time " << float(end_time - begin_time)/CLOCKS_PER_SEC << endl;

  if(vm.count("traj")) {
    mcs.unset_trajectory_file();
  }

  return 0;
}
