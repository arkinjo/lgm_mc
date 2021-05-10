#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>

#include "model.hpp"
#include "PartFunc.hpp"
using namespace std;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h","produce help message")
    ("param,p", po::value<string>(), "input LGM parameter file")
    ("param_out", po::value<string>(), "output LGM parameter file")
    ("logbase,l", po::value<string>(), "basename for log files")
    ("temp", po::value<double>(), "Temperature")
    ("mu_len", po::value<double>(), "Mu_len")
    ("mu_I", po::value<double>(), "Mu_I")
    ("eps", po::value<double>(), "epsilon for convergence")
    ("alpha", po::value<double>(), "for updating mean-fields")
    ("niter,n", po::value<int>(), "max number of iterations")
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
  
  double Temperature = 1.0;
  if(vm.count("temp")) {
    Temperature = vm["temp"].as<double>();
  }
  model.set_temperature(Temperature);
  cerr << "Temperature= " << Temperature << endl;

  double Mu_len = 0;
  if(vm.count("mu_len")) {
    Mu_len = vm["mu_len"].as<double>();
  }
  double Mu_I = 0;
  if(vm.count("mu_I")) {
    Mu_I = vm["mu_I"].as<double>();
  }
  double eps = 1e-12;
  if(vm.count("eps")) {
    eps = vm["eps"].as<double>();
  }
  double alpha = 1.0;
  if(vm.count("alpha")) {
    alpha = vm["alpha"].as<double>();
  }
  int niter = 30000;
  if(vm.count("niter")) {
    niter = vm["niter"].as<int>();
  }
  cerr << "eps,alpha= " << eps << "\t" << alpha << endl;
  cerr << "Mu_len,Mu_I= " << Mu_len << "\t" << Mu_I << endl;
  model.set_mu_len(Mu_len);
  model.set_mu_lenI(Mu_I);

  PartFunc pfunc(model);
    
  const clock_t begin_time = clock();
  pfunc.runSCF(model, niter, eps, alpha);
  const clock_t end_time = clock();

  const double SpecHeat = 0.0;
  double slenO0,slenI0;
  const double slen0 = model.get_seqlength(slenO0,slenI0);
  double slenO,slenI;
  const double slen = model.get_seqlength_exp(slenO,slenI);
  double devO,devI;
  const double dev = model.get_deviation(devO,devI);
  double egibbs0,ebonded0,enonbonded0;
  model.get_energy(ebonded0,enonbonded0,egibbs0);
  const double energy0 = ebonded0 + enonbonded0;
  double egibbs,ebonded,enonbonded;
  model.get_energy_exp(ebonded,enonbonded,egibbs);
  const double energy = ebonded + enonbonded;
  const double ent0 = model.get_site_entropy();
  const double ent = model.get_site_entropy_exp();
  const double divO = model.get_divergence();
  const double divR = model.get_divergence_from_ref();
  cout << "Summary:\t"
       << "inf" << "\t"
       << Mu_I << "\t"
       << Temperature << "\t"
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
	 << model.dev_O[i] << "\t" << model.dev_I[i] << "\t"
	 << model.ent_O[i] << "\t" 
	 << model.div_O[i] << "\t"
      	 << model.div_R[i] << "\t"
	 << model.site_energy_exp[i] << endl;
  }

  if(vm.count("param_out")) {
    const string pout = vm["param_out"].as<string>();
    cerr << "#  Writing parameters to " << pout << endl;
    model.write_parameters(pout);
  }

  return 0;
}
