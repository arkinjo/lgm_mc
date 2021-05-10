#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>

#include "msa.hpp"
#include "model.hpp"
#include "MCseq.hpp"
#include "OptParam.hpp"
#include "AnaTraj.hpp"
using namespace std;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h","produce help message")
    ("dos", po::value<string>(), "WL DOS file")
    ("param,p", po::value<string>(), "input LGM parameter file")
    ("param_out", po::value<string>(), "output LGM parameter file")
    ("logbase,l", po::value<string>(), "basename for log files")
    ("temp", po::value<double>(), "Temperature")
    ("mu_len", po::value<double>(), "Mu_len")
    ("mu_I", po::value<double>(), "Mu_I")
    ("mulcano", po::value<bool>(), "Multicanonical?")
    ("ftraj", po::value<string>(), "trajectory file to be analyzed")
    ("seed", po::value<int>(), "random seed")
    ("min_energy",po::value<double>(), "lower bound of energy range")
    ("max_energy",po::value<double>(), "upper bound of energy range")
    ("delta_energy",po::value<double>(), "max_energy = min_energy + delta_energy")
    ("min_cnt",po::value<int>(), "minimum count of each bin")
    ("ebin",po::value<int>(), "target energy bin")
    ;

  po::positional_options_description p;
  po::variables_map vm;
  po::store(po::command_line_parser(argc,argv).
	    options(desc).positional(p).run(), vm);
  po::notify(vm);

  if(vm.count("help")) {
    cerr << desc << endl;
    return 1;
  }
  string logbase = "ana";
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

  int min_cnt = 100000000;
  if(vm.count("min_cnt")) {
    min_cnt = vm["min_cnt"].as<int>();
  }
  cerr << "minimum count of each bin= " << min_cnt << endl;
  
  bool erange = false;
  double min_energy=-1000, delta_energy=1.0;
  if(vm.count("min_energy")) {
    min_energy = vm["min_energy"].as<double>();
    erange = true;
  }
  if(vm.count("delta_energy")) {
    delta_energy = vm["delta_energy"].as<double>();
  } else if(vm.count("max_energy")) {
    delta_energy = vm["max_energy"].as<double>() - min_energy;
  }
  const double max_energy = min_energy + delta_energy;

  MCseq mcs(model);
  
  if(vm.count("seed")) {
    mcs.set_seed(vm["seed"].as<int>());
    cerr << "random seed= " << (vm["seed"].as<int>()) << endl;;
  }

  mcs.STATS_MULTICANONICAL = true;
  if(vm.count("mulcano")) {
    mcs.STATS_MULTICANONICAL = vm["mulcano"].as<bool>();
  }

  if(vm.count("dos")) {
    mcs.read_DOS(vm["dos"].as<string>(), true);
  } else {
    cerr << "Specify the DOS file with the '--dos' option!!" << endl;
    exit(1);
  }
  double temp = 1.0;
  if(vm.count("temp")) {
    temp = vm["temp"].as<double>();
  }
  model.set_temperature(temp);
  mcs.set_temperature(temp);
  cerr << "Temperature= " << temp << endl;
    
    
  AnaTraj atraj(model);
  string ftraj;
  if(vm.count("ftraj")) {
    ftraj = vm["ftraj"].as<string>();
    cerr << "# Trajectory file: " <<  ftraj << endl;
    atraj.set_trajectory_file(ftraj);
  }

  double Mu_len = 0;
  if(vm.count("mu_len")) {
    Mu_len = vm["mu_len"].as<double>();
  }
  cerr << "Mu_len= " << Mu_len << endl;
  model.set_mu_len(Mu_len);

  double Mu_I = 0;
  if(vm.count("mu_I")) {
    Mu_I = vm["mu_I"].as<double>();
  }
  cerr << "Mu_I= " << Mu_I << endl;
  model.set_mu_lenI(Mu_I);

  int Nseqs=0;
  mcs.init_stats(model);
  mcs.stats_full_nonbonded = false;
  const clock_t begin_time = clock();
  cerr << "# Sampling from Trajectory file " << endl;

  mcs.clear_EHist();
  
  while(atraj.pick_seq(model, mcs)) {
    if(mcs.Energy < mcs.MIN_ENERGY || mcs.Energy >= mcs.MAX_ENERGY) continue;
    
    const int bin = mcs.get_DOS_bindex(mcs.Energy);
    if(mcs.EHist[bin]>min_cnt) continue;
    mcs.EHist[bin]++;
    const double div = mcs.get_divergence(model);
    if(!erange
       || (erange && mcs.Energy >= min_energy && mcs.Energy < max_energy)) {
      const double e = mcs.Energy;
      const double w = mcs.count_stats(model);
      Nseqs++;
      if(Nseqs % 10000 == 0) cerr << Nseqs << " sequences processed." << endl;
    }
  }
  atraj.unset_trajectory_file();
  cerr << "total weight= " << mcs.WeightTotal << endl;
  mcs.finish_stats(model);
  const clock_t end_time = clock();

  int imin = -1;
  double ehist_min = 1e+100;
  for(map<int,double>::iterator it = mcs.DOS.begin();
      it != mcs.DOS.end(); ++it) {
    if(mcs.EHist[it->first] < ehist_min) {
      imin = it->first;
      ehist_min = mcs.EHist[it->first];
    }
  }
  cerr << "min_EHist/min_cnt= "
       << imin
       << "\t" << imin*mcs.DELTA_ENERGY
       << "\t" << ehist_min << "\t" << min_cnt << endl;
  
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
  const double ent0 = model.get_site_entropy();
  const double ent = model.get_site_entropy_exp();
  const double divO = model.get_divergence();
  const double divR = model.get_divergence_from_ref();
  cout << "Summary:\t"
       << Nseqs << "\t"
       << Mu_I << "\t"
       << mcs.Temperature << "\t"
       << SpecHeat << "\t"
       << slen - slen0 << "\t"
       << slenO - slenO0 << "\t"
       << slenI - slenI0 << "\t"
       << dev << "\t"
       << devO << "\t"
       << devI << "\t"
       << energy - energy0 << "\t"
       << ebonded - ebonded0 << "\t"
       << enonbonded - enonbonded0 << "\t"
       << egibbs << "\t"
       << ent - ent0 << "\t"
       << divO << "\t"
       << divR << endl;

  for(int i = 1; i <= model.length(); ++i) {
    cout << "Site\t" << i << "\t"
	 << model.ent_O[i] << "\t" 
	 << model.div_O[i] << "\t"
      	 << model.ddiv_O[i] << "\t"
	 << model.div_OOs[i] << "\t"
      	 << model.ddiv_OOs[i] << "\t"
	 << model.mdiv_O[i] << "\t"
	 << model.site_energy_exp[i] << "\t"
	 << model.dev_I[i] << endl;
  }

  if(vm.count("param_out")) {
    const string pout = vm["param_out"].as<string>();
    cerr << "#  Writing parameters to " << pout << endl;
    model.write_parameters(pout);
    const string covout = logbase + ".cov";
    cerr << "#  Writing covariance matrix to " << covout << endl;
    model.dump_covariance_exp(covout);
  }

return 0;
}
