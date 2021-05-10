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
    ("nloop", po::value<int>(), "number of optimization loops")
    ("nsweep", po::value<int>(), "number of sweeps")
    ("temp_min", po::value<double>(), "lowest temperature")
    ("temp_max", po::value<double>(), "highest temperature")
    ("min_energy", po::value<double>(), "minimun energy")
    ("max_energy", po::value<double>(), "maximum energy")
    ("delta_energy", po::value<double>(), "width of energy bin")
    ("mu_len", po::value<double>(), "Mu_len")
    ("mu_I", po::value<double>(), "Mu_I")
    ("eps,e", po::value<double>(), "precision of DOS")
    ("dos", po::value<string>(), "WL DOS file for continuation!")
    ("wlfactor,w", po::value<double>(), "WLfactor for (re)starting")
    ("seed", po::value<int>(), "random seed")
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


  
  if(!vm.count("param")) {
    cerr << "Specify the parameter file!" << endl;
    exit(1);
  }
  
  Model model;
  model.read_parameters(vm["param"].as<string>());
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

    
  MCseq mcs(model);
  if(vm.count("seed")) {
    mcs.set_seed(vm["seed"].as<int>());
    cerr << "random seed= " << (vm["seed"].as<int>()) << endl;;
  }
  double tmax=1.2,tmin=0.5;
  if(vm.count("temp_max")) {
    tmax = (vm["temp_max"].as<double>());
  }
  if(vm.count("temp_min")) {
    tmin = (vm["temp_min"].as<double>());
  }
  cerr << "Temperature range: " << tmin << " ... " << tmax << endl;

  mcs.select_type = CANONICAL;

  mcs.NMONITOR = 10000; 
  mcs.gibbs_sampler_p = true;
  const int nsweeps_eq = 10000;
  double min_energy=0;
  if(vm.count("min_energy")) {
    min_energy = vm["min_energy"].as<double>();
  } else if(!vm.count("dos")) {
    mcs.init_stats(model);
    for(int k = 0; k < 10; ++k) {
      mcs.set_random_sequence(model);
      mcs.set_temperature(1.0);
      mcs.run(model, 10000, false);
      mcs.anneal(model,1000, 100, 1.0, tmin);
      mcs.set_temperature(tmin);
      mcs.run(model, 10000, true);
    }
    mcs.finish_stats(model);
    cerr << "Tmin " << tmin
	 << "\t" << mcs.Eave << "\t" << sqrt(mcs.Evar) << endl;
    min_energy = floor(mcs.Eave);
  }

  double max_energy=100;
  if(vm.count("max_energy")) {
    max_energy = vm["max_energy"].as<double>();
  } else if(!vm.count("dos")) {
    mcs.init_stats(model);
    mcs.set_temperature(tmax);
    for(int k = 0; k < 100; ++k) {
      mcs.set_random_sequence(model);
      mcs.run(model, 1000, false);
      mcs.run(model, 1000, true);
    }
    mcs.finish_stats(model);
    cerr << "Tmax " << tmax
	 << "\t" << mcs.Eave << "\t" << sqrt(mcs.Evar) << endl;
    max_energy = ceil(mcs.Eave);
  }

  double delta_energy = 0.2;
  if(vm.count("delta_energy")) {
    delta_energy = vm["delta_energy"].as<double>();
  }
  mcs.setup_WL(min_energy, max_energy, delta_energy);

  if(vm.count("dos")) {
    mcs.read_DOS(vm["dos"].as<string>(), true);
  } 
  if(vm.count("wlfactor")) {
    mcs.WLfactor = vm["wlfactor"].as<double>();
  } 
  min_energy = mcs.MIN_ENERGY;
  max_energy = mcs.MAX_ENERGY;

    
  double eps = 1e-5;
  if(vm.count("eps")) {
    eps = vm["eps"].as<double>();
  }
  mcs.WLprec = eps;
    

  int nsweeps = 100000;
  if(vm.count("nsweep")) {
    nsweeps = vm["nsweep"].as<int>();
  }


  mcs.NMONITOR = 10000; //nsweeps / 10 ;


  int nloops = 10;
  if(vm.count("nloop")) {
    nloops = vm["nloop"].as<int>();
  }
    
  const clock_t begin_time = clock();
  cerr << "# WL sampling run for "
       << nloops << " loops of " << nsweeps << " sweeps." << endl;
  int iloop,isweeps;

  for(iloop = 1; iloop <= nloops; ++iloop) {
    cerr << "#Loop: " << iloop << endl;
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
    isweeps = mcs.run(model, nsweeps, false);
    mcs.save_DOS(logbase + ".dos");
    if(isweeps < nsweeps) break;
  }
  const clock_t end_time = clock();
  cerr << "Time " << float(end_time - begin_time)/CLOCKS_PER_SEC << endl;
  if(iloop == nloops + 1) {
    cerr << "### Wang-Landau sampling did NOT converge!" << endl;
  } else {
    cerr << "### Wang-Landau sampling CONVERGED in "
	 << iloop << " loops and "
	 << isweeps << " sweeps!" << endl;
  }
return 0;
}
