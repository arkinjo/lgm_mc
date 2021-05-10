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
    ("in,i", po::value<string>(), "set input file (multiple FASTA)")
    ("dos", po::value<string>(), "WL DOS file")
    ("param,p", po::value<string>(), "input LGM parameter file")
    ("param_out", po::value<string>(), "output LGM parameter file")
    ("logbase,l", po::value<string>(), "basename for log files")
    ("temp", po::value<double>(), "Temperature")
    ("mu_len", po::value<double>(), "Mu_len")
    ("ftraj", po::value<string>(), "trajectory file to be analyzed")
    ("feigen", po::value<string>(), "eigen components file.")
    ("seed", po::value<int>(), "random seed")
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

  if(vm.count("in")) {
    cerr << "Input file was set to " << vm["in"].as<string>() << endl;
    string infile = vm["in"].as<string>();
    Msa msa(infile);
    Model model(msa);
    int nseq = model.get_num_seq();
    if(vm.count("param")) {
      model.read_parameters(vm["param"].as<string>());
    } else {
      cerr << "Specify the parameter file!" << endl;
      exit(1);
    }

    
    MCseq mcs(model);

    if(vm.count("seed")) {
      mcs.set_seed(vm["seed"].as<int>());
      cerr << "random seed= " << (vm["seed"].as<int>()) << endl;;
    }
    if(vm.count("dos")) {
      mcs.STATS_MULTICANONICAL = true;
      mcs.read_DOS(vm["dos"].as<string>(), true);
    } else {
      mcs.STATS_MULTICANONICAL = false;
      cerr << "Trajectory treated as CANONICAL." << endl;
    }
    double temp = 1.0;
    if(vm.count("temp")) {
      temp = vm["temp"].as<double>();
    }
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
    
    int Nseqs=0;
    mcs.init_stats(model);
    mcs.stats_full_nonbonded = false;
    const clock_t begin_time = clock();
    cerr << "# Sampling from Trajectory file " << endl;
    vector<double> weights;
    while(atraj.pick_seq(model, mcs)) {
      const double w = mcs.count_stats(model, true);
      weights.push_back(w);
      Nseqs++;
      if(Nseqs % 100000 == 0) cerr << Nseqs << " sequences processed." << endl;
    }
    atraj.unset_trajectory_file();
    cerr << "total weight= " << mcs.WeightTotal << endl;
    mcs.finish_stats(model);
    const clock_t end_time = clock();

    if(vm.count("feigen")) {
      const string feigen = vm["feigen"].as<string>();
      cerr << "# Eigen-components file: " <<  feigen << endl;
      atraj.read_eigen_components(feigen);
    }
    const int nseqs = weights.size();
    for(int i = 0; i != nseqs; ++i) {
      weights[i] /= mcs.WeightTotal;
    }
    vector<int> inds = PCA_project_seqs(atraj, cout, 10000, weights);
    atraj.set_trajectory_file(ftraj);
    int ind=0;
    while(atraj.pick_seq(model, mcs)) {
      if(std::binary_search (inds.begin(), inds.end(), ind)) {
	cerr << "weight " << ind << "\t" << weights[ind] << endl;
	darray1 p = atraj.pca_project(mcs);
	cout << p[0] << "\t"
	     << p[1] << "\t"
	     << p[2] << "\t"
	     << mcs.Energy << "\t"
	     << weights[ind] << endl;
      }
      ++ind;
    }
    atraj.unset_trajectory_file();
  }
  return 0;
}
