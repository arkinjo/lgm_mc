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

const int CORRELATION_CUTOFF = 0;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h","produce help message")
    ("in,i", po::value<string>(), "set input MSA file (multiple FASTA)")
    ("hmm", po::value<bool>(), "input MSA file is Pfam HMM alignment[true]")
    ("ipair", po::value<string>(), "interacting pairs of reference sequence")
    ("rmap", po::value<string>(), "HMMer3 alignment of reference sequence")
    ("pseed,s", po::value<string>(), "LGM parameter file (for restart)")
    ("param,p", po::value<string>(), "LGM parameter file (for restart)")
    ("wlopt", po::value<bool>(), "use WL sampling for optimization")
    ("epep", po::value<double>(), "Default energy parameter for peptide bond")
    ("logbase,l", po::value<string>(), "basename for log files")
    ("alphaH", po::value<double>(), "alphaH")
    ("alphaJ", po::value<double>(), "alphaJ")
    ("alphaK", po::value<double>(), "alphaK")
    ("betaH", po::value<double>(), "betaH")
    ("betaJ", po::value<double>(), "betaJ")
    ("betaK", po::value<double>(), "betaK")
    ("gibbs", po::value<bool>(), "use Gibbs sampler or not")
    ("nsweep", po::value<int>(), "number of sweeps")
    ("nsweep_eq", po::value<int>(), "number of sweeps for equilibration")
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

  if(vm.count("in")) {
    cerr << "Input file was set to " << vm["in"].as<string>() << endl;
    string infile = vm["in"].as<string>();
    Msa msa(infile,false);
    bool hmm=true;
    if(vm.count("hmm")) {
      hmm = vm["hmm"].as<bool>();
    }

    Model model(msa,hmm);
    double epep = 0.0;
    if(vm.count("pseed")) {
      model.read_parameters(vm["pseed"].as<string>());
      cerr << "Seed parameters read from " << vm["pseed"].as<string>() << endl;
      model.stat_all();
    }

    if(vm.count("epep")) {
      epep = vm["epep"].as<double>();
    }
    cerr << "Default energy parameter for peptide bond: " << epep << endl;
    cerr << "Gamma = " << model.get_gamma() << endl;
    
    int nseq = model.get_num_seq();
    for (int k = 0; k < nseq; ++k) {
      cerr << msa.header(k) << "\t" << model.get_weight(k) << "\t"
	   << msa.get_weight(k) << endl;
    }
    //    model.dump_premodel(msa,"hogehoge_premodel.dat");
    if(vm.count("ipair")) {
      cerr << "#reading ipair " << (vm["ipair"].as<string>()) << endl;
      model.read_ipairs(vm["ipair"].as<string>());
      if(!hmm) {
	model.align_ipairs_1st(msa);
      }
      else if(vm.count("rmap")) {
      	cerr << "#reading rmap " << (vm["rmap"].as<string>()) << endl;
      	model.read_reference_3Dalignment(vm["rmap"].as<string>());
      } else {
      	cerr << "Give non-hmm alignment or --rmap for 3D contact mapping!" << endl;
      	exit(2);
      }
    } else {
      cerr << "Setting default ipairs." << endl;
      model.set_default_ipairs();
    }

    if(vm.count("param")) {
      model.read_parameters(vm["param"].as<string>());
      cerr << "Parameters read from " << vm["param"].as<string>() << endl;
    }

    if(CORRELATION_CUTOFF==1) {
      model.NS_CUTOFF = 1.0e-4;
      model.NL_CUTOFF = 1.0e-3;
    } else if (CORRELATION_CUTOFF==2) {
      model.NS_CUTOFF = 1.0;
      model.NL_CUTOFF = 1.0;
    } else{
      model.NS_CUTOFF = -1;
      model.NL_CUTOFF = -1;
    }
    model.init_params(epep);
    //    model.set_nl_cutoff(model.model_length*nstate_core*10);
    
    MCseq mcs(model);
    mcs.select_type = CANONICAL;
    mcs.set_random_sequence(model,0);
    mcs.get_energy(model);
    mcs.gibbs_sampler_p = true;
    if(vm.count("gibbs")) {
      mcs.gibbs_sampler_p = vm["gibbs"].as<bool>();
    }
    cerr << "# Sampling method: "
	 << ((mcs.gibbs_sampler_p) ? "Gibbs" : "Metropolis") << endl;

    OptParam optpar(model);
    cerr << "AASEQ0\t";
    mcs.dump_aaseq();
    optpar.nsweeps_eq = 10000;
    if(vm.count("nsweep_eq")) {
      optpar.nsweeps_eq = vm["nsweep_eq"].as<int>();
    }
    optpar.nsweeps = 100000;
    if(vm.count("nsweep")) {
      optpar.nsweeps = vm["nsweep"].as<int>();
    }
    mcs.NMONITOR = optpar.nsweeps;
    
    if(vm.count("wlopt")) {
      optpar.use_WL = vm["wlopt"].as<bool>();
    }
    if(optpar.use_WL) {
      cerr << "# using WL sampling for optimization" << endl;
      mcs.NMONITOR = optpar.nsweeps/10;
    }


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

    cerr << "# Warming up..." << endl;
    mcs.run(model, 10000, false);
    mcs.init_stats(model);
    mcs.run(model, 10000, true);
    mcs.finish_stats(model);
    const clock_t begin_time = clock();
    for(int i = 1; i <= nloop; ++i) {
      if(optpar.use_WL) {
	optpar.WLoptimize(model, mcs, 10);// DEBUGGING
      } else {
	//optpar.optimize(model, mcs, 10);
	optpar.optimize_mul(model, mcs, 10);
      }
      model.write_parameters(logbase + ".param0");
    }
    const clock_t end_time = clock();
    cerr << "Time " << float(end_time - begin_time)/CLOCKS_PER_SEC << endl;

    model.write_parameters(logbase + ".param");
    //    model.write_parameters(logbase + ".param", true);
    optpar.write_momentum(model, logbase + ".mmt");
    
    cerr << "max adel (site) " << model.check_site(false) << endl;
    cerr << "max adel (bonded) " << model.check_bonded(false) << endl;
    cerr << "max adel (nonbonded) " << model.check_nonbonded(false) << endl;
  }

  return 0;
}
