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
    ("in,i", po::value<string>(), "set input file (multiple FASTA)")
    ("ipair", po::value<string>(), "interacting pairs of reference sequence")
    ("rmap", po::value<string>(), "HMMer3 alignment of reference sequence")
    ("param,p", po::value<string>(), "LGM parameter file (for restart)")
    ("logbase,l", po::value<string>(), "basename for log files")
    ("alphaH", po::value<double>(), "alphaH")
    ("alphaJ", po::value<double>(), "alphaJ")
    ("alphaK", po::value<double>(), "alphaK")
    ("betaH", po::value<double>(), "betaH")
    ("betaJ", po::value<double>(), "betaJ")
    ("betaK", po::value<double>(), "betaK")
    ("nloop", po::value<int>(), "number of optimization loops")
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
  string logbase = "train";
  if(vm.count("logbase")) {
    logbase = vm["logbase"].as<string>();
  }

  if(!vm.count("in")) {
    cerr << "Specify the input MSA!" << endl;
    exit(1);
  }
  
  cerr << "Input file was set to " << vm["in"].as<string>() << endl;
  string infile = vm["in"].as<string>();
  Msa msa(infile);
  Model model(msa);

  if(vm.count("ipair")) {
    cerr << "#reading ipair " << (vm["ipair"].as<string>()) << endl;
    model.read_ipairs(vm["ipair"].as<string>());
    if(vm.count("rmap")) {
      cerr << "#reading rmap " << (vm["rmap"].as<string>()) << endl;
      model.read_reference_3Dalignment(vm["rmap"].as<string>());
    } else {
      cerr << "Give --rmap for 3D contact mapping!" << endl;
      exit(2);
    }
  } else {
    cerr << "Setting default ipairs." << endl;
    model.set_default_ipairs();
  }
  
  if(vm.count("param")) {
    model.read_parameters(vm["param"].as<string>());
    cerr << "Parameters read from " << vm["param"].as<string>() << endl;
  } else {
    cerr << "Parameter file (with optimized J) must be specified!" << endl;
    exit(1);
  }
    
  OptParam optpar(model);
  int nloop = 20;
  if(vm.count("nloop")) {
    nloop = vm["nloop"].as<int>();
  }
  cerr << "# optimizing for " << nloop << " loops." << endl;
    
  if(vm.count("alphaH")) {
    optpar.alpha_H = vm["alphaH"].as<double>();
  }
  if(vm.count("alphaJ")) {
    optpar.alpha_J = vm["alphaJ"].as<double>();
  }
  if(vm.count("alphaK")) {
    optpar.alpha_K = vm["alphaK"].as<double>();
  }
  cerr << "# alpha_{H,J,K}= "
       << optpar.alpha_H << "\t"
       << optpar.alpha_J << "\t"
       << optpar.alpha_K << endl;

  if(vm.count("betaH")) {
    optpar.beta_H = vm["betaH"].as<double>();
  }
  if(vm.count("betaJ")) {
    optpar.beta_J = vm["betaJ"].as<double>();
  }
  if(vm.count("betaK")) {
    optpar.beta_K = vm["betaK"].as<double>();
  }
  cerr << "# beta_{H,J,K}= "
       << optpar.beta_H << "\t"
       << optpar.beta_J << "\t"
       << optpar.beta_K << endl;

  cerr << "energy(obs)= " << model.get_energy() << endl;
  
  const string param_out = logbase + ".param";
  const clock_t begin_time = clock();
  optpar.optimize_plm(model,nloop, param_out);
  const clock_t end_time = clock();
  cerr << "Time " << float(end_time - begin_time)/CLOCKS_PER_SEC << endl;
    

  //  optpar.write_momentum(model, logbase + ".mmt");

  cerr << "energy(obs)= " << model.get_energy() << endl;
  cerr << "energy(exp)= " << model.get_energy_exp() << endl;

  cerr << "max adel (site) " << model.check_site(false) << endl;
  cerr << "max adel (bonded) " << model.check_bonded(false) << endl;
  cerr << "max adel (nonbonded) " << model.check_nonbonded(false) << endl;

  return 0;
}
