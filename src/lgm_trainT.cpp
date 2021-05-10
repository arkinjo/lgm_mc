/*
  Training from MC trajectories
 */
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
    ("param,p", po::value<string>(), "input LGM parameter file")
    ("logbase,l", po::value<string>(), "basename for log files")
    ("ftraj", po::value<string>(), "trajectory file to be parsed")
    ("fmmt", po::value<string>(), "momentum file to be parsed")
    ("alphaH", po::value<double>(), "alphaH")
    ("alphaJ", po::value<double>(), "alphaJ")
    ("alphaK", po::value<double>(), "alphaK")
    ("betaH", po::value<double>(), "betaH")
    ("betaJ", po::value<double>(), "betaJ")
    ("betaK", po::value<double>(), "betaK")
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
  
  AnaTraj atraj(model);
  string ftraj;
  if(vm.count("ftraj")) {
    ftraj = vm["ftraj"].as<string>();
    cerr << "# Trajectory file: " <<  ftraj << endl;
    atraj.set_trajectory_file(ftraj);
  } else {
    cerr << "Specify a trajectory file!" << endl;
    exit(1);
  }
    
  OptParam optpar(model);
  string fmmt;
  if(vm.count("fmmt")) {
    fmmt = vm["fmmt"].as<string>();
    cerr << "reading momentum from " << fmmt << endl;
    optpar.read_momentum(model,fmmt);
  }
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

  optpar.optimize_traj(model, atraj);
  atraj.unset_trajectory_file();
  model.write_parameters(logbase + ".param");
  optpar.write_momentum(model, logbase + ".mmt");

  return 0;
}
