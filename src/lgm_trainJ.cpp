#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>

#include "msa.hpp"
#include "model.hpp"
#include "OptParam.hpp"
#include "PartFunc.hpp"
using namespace std;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h","produce help message")
    ("in,i", po::value<string>(), "set input file (multiple FASTA)")
    ("default,d", po::value<bool>(), "Set default (NULL) sequence for training")
    ("ipair", po::value<string>(), "interacting pairs of reference sequence")
    ("rmap", po::value<string>(), "HMMer3 alignment of reference sequence")
    ("param,p", po::value<string>(), "input LGM parameter file")
    ("alpha", po::value<double>(), "alpha (for update)")
    ("niter,n", po::value<int>(), "number of iterations")
    ("logbase,l", po::value<string>(), "basename for log files")
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
  string logbase = "exJ";
  if(vm.count("logbase")) {
    logbase = vm["logbase"].as<string>();
  }
  cerr << "# LOGBASE= " << logbase << endl;
  if(!vm.count("in")) {
    cerr << "Specify the input MSA!" << endl;
    exit(1);
  }
  
  cerr << "Input file was set to " << vm["in"].as<string>() << endl;
  string infile = vm["in"].as<string>();
  Msa msa(infile);
  Model model(msa);
  if(vm.count("param")) {
    model.read_parameters(vm["param"].as<string>());
    cerr << "Setting parameters from " << vm["param"].as<string>() << endl;
  }
  
  if(vm.count("default")) {
    if(vm["default"].as<bool>()) {
      //      model.set_null_sequence();
      model.set_uniform_sequence();
      cerr << "Setting default NULL sequence for training" << endl;
    }
  }

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
    model.init_KOO_corr();
  } else {
    cerr << "Setting NO ipairs." << endl;
  }

  PartFunc pfunc(model);
  double alpha = 0.1;
  if(vm.count("alpha")) {
    alpha = vm["alpha"].as<double>();
  }
  cerr << "alpha= " << alpha << endl;
  int niter = 300;
  if(vm.count("niter")) {
    niter = vm["niter"].as<int>();
  }
  cerr << "Optimizing for " << niter << " iterations." << endl;
  pfunc.learn_J(model, niter, alpha);
  model.write_parameters(logbase + ".param");

  cerr << "energy(obs)= " << model.get_energy() << endl;
  cerr << "energy(exp)= " << model.get_energy_exp() << endl;

  cerr << "max adel (site) " << model.check_site(false) << endl;
  cerr << "max adel (bonded) " << model.check_bonded(false) << endl;
  cerr << "max adel (nonbonded) " << model.check_nonbonded(false) << endl;

  return 0;
}
