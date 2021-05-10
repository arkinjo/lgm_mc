#include "MCseq.hpp"

boost::uniform_01<> dist01; // for Metropolis
boost::random::uniform_int_distribution<> rand_core(0,nstate_core-1),
  rand_insert(0,nstate_insert-1);

ISite::ISite()
{
  O = terminal_state;
  vector<int> v;
  I= v;
}


bool MCseq::check_EHist(void)
{
  const bool verbose = false;
  double cmin=1.e+100,cmax=0.0, ave=0.0;
  int imin=-1, imax=-1;

  for(map<int,double>::iterator it = EHist.begin(); it != EHist.end(); ++it) {
    const int bind = it->first;
    const double cnt = it->second;
    ave += cnt;
    if(cnt < cmin) { cmin = cnt; imin = bind;}
    if(cnt > cmax) { cmax = cnt; imax = bind;}

  }
  ave /= EHist.size();
  MAX_DOS = 0;
  MIN0_ENERGY = 1e+100;
  for(map<int,double>::iterator it = DOS.begin(); it != DOS.end(); ++it) {
    const double ene = it->first * DELTA_ENERGY;
    if(ene < MIN0_ENERGY) MIN0_ENERGY = ene;
    if(it->second > MAX_DOS) MAX_DOS = it->second;
  }
  
  if(verbose) {
    cerr << "#Wang-Landau updated (imin,cmin),(imax,cmax)= ("
	 << imin << "," << cmin << ")\t(" << imax << "," << cmax << ")" << endl;
  }


  //  if(ave > 0.0 && cmin/ave > 0.9 && cmax/ave < 1.1)
  if(cmin > 0.0 && cmin/cmax > 0.8)
    {
      WLfactor *= 0.5; // this may be improved.
      cerr << "#MCseq::check_EHist WLfactor modified to " << WLfactor
	   << "\tc(min/max)/ave,cmin/cmax= " << cmin/ave << "\t" << cmax/ave
	   << "\t" << cmin/cmax
	   << endl;
      return true;
    }
  return false;
}


void MCseq::clear_DOS(void)
{
  DOS.clear();
  WLfactor=1.0;
  return;
}

void MCseq::clear_EHist(void)
{
  EHist.clear();
  EHistTraj.clear();
  for(map<int,double>::iterator it = DOS.begin(); it != DOS.end(); ++it) {
    EHist[it->first] = 0;
    EHistTraj[it->first] = 0;
  }
  return;
}

inline bool MCseq::Metropolis(const double dene, const double coef)
{
  return (dene <= 0.0 || dist01(rng) < coef*exp(-dene*Beta));
}

inline bool MCseq::WL(const double dene, const double coef)
{
  const double ene1 = Energy + dene;
  const int i0 = get_DOS_bindex(Energy);
  const int i1 = get_DOS_bindex(ene1);

  if(Energy < MIN_ENERGY || Energy >= MAX_ENERGY) {
    cerr << "MCseq::WL: Energy out of range " << Energy << " not in "
	 << MIN_ENERGY << " .. " << MAX_ENERGY << endl;
    exit(1);
  }
  if(DOS.count(i0) == 0) DOS[i0] = EHist[i0] = 0;
  if(ene1 >= MAX_ENERGY || ene1 < MIN_ENERGY) {
    EHist[i0] += 1;
    DOS[i0] += WLfactor;
    return false;
  }
  
  if(DOS.count(i1) == 0) DOS[i1] = EHist[i1] = 0;
  
  double dh = DOS[i1] - DOS[i0];

  if (dh <= 0.0 || dist01(rng) < coef*exp(-dh)) {
    EHist[i1] += 1;
    DOS[i1] += WLfactor;
    //    cerr << "WL1 " << i1 << "\t" << EHist[i1] << endl;
    return true;
  } else {
    EHist[i0] += 1;
    DOS[i0] += WLfactor;
    //    cerr << "WL2 " << i0 << "\t" << EHist[i0] << endl;
    return false;
  }
}

inline bool MCseq::select_criteria(const double dene, const double coef)
{
  if(select_type == CANONICAL) {
    return Metropolis(dene, coef);
  } else {
    return WL(dene, coef);
  }
}


MCseq::MCseq(Model& model)
{
  rng.seed(135193);
  set_model_length(model.length());
  WeightTotal = 0.0;
  trajectory_p = false;
  select_type = CANONICAL;
  NMONITOR=100;
  NTRAJ=1;
  MONITOR_p = true;
  Temperature = model.temperature; Beta = 1.0/Temperature;
  gibbs_sampler_p = false;
  MAX_ENERGY = -1.0e+100;
  DELTA_ENERGY = 1.0;
  MIN_ENERGY = 1e+100;
  WLfactor = 1.0;
  MAX_DOS = 0;
  WLprec = 1e-5;

  STATS_MULTICANONICAL = false;
  stats_full_nonbonded = false;

  set_random_sequence(model, 0);
  
  IntPairs &intpairs = model.get_intpairs();
  for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
    const int i = it->i, j = it-> j;
    if(IntList.count(i) == 0) {
      vector<int> v;
      IntList[i] = v;
    }
    if(IntList.count(j) == 0) {
      vector<int> v;
      IntList[j] = v;
    } 
    IntList[i].push_back(j);
    IntList[j].push_back(i);
  }
  
}

void MCseq::set_trajectory_file(string filename)
{
  ftraj.open(filename.c_str());
  ftraj.setf(std::ios::scientific);
  trajectory_p = true;
  return;
}

int MCseq::traj_min_count(void)
{
  int min_cnt = 1000000000, max_cnt=0;
  int imin=-1,imax=-1;
  for(map<int,double>::iterator it = EHistTraj.begin();
      it != EHistTraj.end(); ++it) {
    if(it->second < min_cnt) {
      min_cnt = it->second;
      imin = it->first;
    } else if (it->second > max_cnt) {
      max_cnt = it->second;
      imax = it->first;
    }
  }
  cerr << "traj_min_count(min/max): "
       << "\t" << imin
       << "\t" << imin*DELTA_ENERGY
       << "\t" << min_cnt
       << "\t" << imax
       << "\t" << imax*DELTA_ENERGY
       << "\t" << max_cnt
       << endl;
  
  return min_cnt;
}

void MCseq::save_trajectory(void)
{
  if(select_type == MULTI_CANONICAL &&
     Energy >= MIN_ENERGY && Energy < MAX_ENERGY) {
    const int bin = get_DOS_bindex(Energy);
    EHistTraj[bin] += 1;
  }
  
  ftraj << Energy << "\t";
  for(int i = 1; i <= Length; ++i) {
    ftraj << amino1_core[aaseq[i].O];
    if(i == Length) continue;
    for(vector<int>::iterator it = aaseq[i].I.begin(); it != aaseq[i].I.end(); ++it) {
      ftraj << amino1_insert[*it];
    }

  }
  ftraj << endl;
}

double MCseq::mutate_core(const Model& model, const int i, const int a1)
{
  int a0 = aaseq[i].O;

  if(a0 == a1) return 0.0;

  double ene0 = model.H_O[i][a0];
  double ene1 = model.H_O[i][a1];
  
  if(i > 1) {
    const int ilen = aaseq[i-1].I.size();
    if(ilen > 0) {
      int b = aaseq[i-1].I[ilen-1];
      ene0 += model.J_IO[i-1][b][a0];
      ene1 += model.J_IO[i-1][b][a1];
    } else {
      int b = aaseq[i-1].O;
      ene0 += model.J_OO[i-1][b][a0];
      ene1 += model.J_OO[i-1][b][a1];
    }
  }
  
  if(i < Length) {
    if(aaseq[i].I.size() > 0) {
      int b = aaseq[i].I[0];
      ene0 += model.J_OI[i][a0][b];
      ene1 += model.J_OI[i][a1][b];
    } else {
      int b = aaseq[i+1].O;
      ene0 += model.J_OO[i][a0][b];
      ene1 += model.J_OO[i][a1][b];
    }
  }
  
  if(IntList.count(i) > 0) {
    int n = IntList[i].size();
    for(int k = 0; k < n; ++k) {
      int j = IntList[i][k];
      int b = aaseq[j].O;
      ene0 += model.K_OO[i][j][a0][b];
      ene1 += model.K_OO[i][j][a1][b];
    }
  }

  return (ene0 - ene1);
}

double MCseq::mutate_core_gs(const Model& model, const int i, int& a1)
{
  const int ilen0 = aaseq[i-1].I.size();
  const int ilen1 = aaseq[i].I.size();
  const int a0 = aaseq[i].O;
  
  for(int a = 0; a < nstate_core; ++a) {
    ep_core[a] = model.H_O[i][a];
    if(i > 1) {
      if(ilen0 > 0) {
	int b = aaseq[i-1].I[ilen0-1];
	ep_core[a] += model.J_IO[i-1][b][a];
      } else {
	int b = aaseq[i-1].O;
	ep_core[a] += model.J_OO[i-1][b][a];
      }
    }
    
    if(i < Length) {
      if(ilen1 > 0) {
	int b = aaseq[i].I[0];
	ep_core[a] += model.J_OI[i][a][b];
      } else {
	int b = aaseq[i+1].O;
	ep_core[a] += model.J_OO[i][a][b];
      }
    }
    
    if(IntList.count(i) > 0) {
      int n = IntList[i].size();
      for(int k = 0; k < n; ++k) {
	int j = IntList[i][k];
	int b = aaseq[j].O;
	ep_core[a] += model.K_OO[i][j][a][b];
      }
    }
  }


  double t;
  double e0 = ep_core[a0];
  double w = 1.0;
  for(int a = 0; a < nstate_core; ++a) {
    t = exp(Beta*(ep_core[a] - e0));
    w += t;
    dist_core[a] = t;
  }
  
  for(int a = 0; a < nstate_core; ++a) {
    dist_core[a] /= w;
  }

  boost::random::discrete_distribution<> dist(dist_core);

  a1 = dist(rng);
  return (ep_core[a0] - ep_core[a1]);
}

double MCseq::mutate_insert(const Model& model, const int i, const int j, const int a1)
{
  const int iend = aaseq[i].I.size() - 1;
  const int a0 = aaseq[i].I[j];

  if(a0 == a1) return 0.0;
  
  double ene0, ene1;

  ene0 = model.H_I[i][a0];
  ene1 = model.H_I[i][a1];
  
  if(iend == 0) {
    ene0 += model.J_OI[i][aaseq[i].O][a0] + model.J_IO[i][a0][aaseq[i+1].O];
    ene1 += model.J_OI[i][aaseq[i].O][a1] + model.J_IO[i][a1][aaseq[i+1].O];
  } else if(j == 0) {
    ene0 += model.J_OI[i][aaseq[i].O][a0] + model.J_II[i][a0][aaseq[i].I[j+1]];
    ene1 += model.J_OI[i][aaseq[i].O][a1] + model.J_II[i][a1][aaseq[i].I[j+1]];
  } else if (j == iend) {
    ene0 += model.J_II[i][aaseq[i].I[j-1]][a0] + model.J_IO[i][a0][aaseq[i+1].O];
    ene1 += model.J_II[i][aaseq[i].I[j-1]][a1] + model.J_IO[i][a1][aaseq[i+1].O];
  } else {
    ene0 += model.J_II[i][aaseq[i].I[j-1]][a0] + model.J_II[i][a0][aaseq[i].I[j+1]];
    ene1 += model.J_II[i][aaseq[i].I[j-1]][a1] + model.J_II[i][a1][aaseq[i].I[j+1]];
  }

  return (ene0 - ene1);
}

double MCseq::mutate_insert_gs(const Model& model, const int i, const int j,
			       int& a1)
{
  const int iend = aaseq[i].I.size() - 1;
  const int a0 = aaseq[i].I[j];
  
  for(int a = 0; a < nstate_insert; ++a) {
    ep_insert[a] = model.H_I[i][a];
    if(iend == 0) {
      ep_insert[a] += model.J_OI[i][aaseq[i].O][a] +  model.J_IO[i][a][aaseq[i+1].O];
    } else if (j == 0) {
      ep_insert[a] += model.J_OI[i][aaseq[i].O][a] + model.J_II[i][a][aaseq[i].I[j+1]];
    } else if (j == iend) {
      ep_insert[a] += model.J_II[i][aaseq[i].I[j-1]][a] + model.J_IO[i][a][aaseq[i+1].O];
    } else {
      ep_insert[a] += model.J_II[i][aaseq[i].I[j-1]][a] + model.J_II[i][a][aaseq[i].I[j+1]];
    }
  }

  double w=0.0,t;
  for(int a = 0; a < nstate_insert; ++a) {
    t = exp(Beta*(ep_insert[a] - ep_insert[a0]));
    dist_insert[a] = t;
    w += t;
  }
  
  for(int a = 0; a <= nstate_insert; ++a) {
    dist_insert[a] /= w;
  }
  boost::random::discrete_distribution<> dist(dist_insert);

  a1 = dist(rng);
  return (ep_insert[a0] - ep_insert[a1]);
}

double MCseq::extend_insert(const Model& model, const int i, const int j, const int a1)
{ 
  const int ilen = aaseq[i].I.size();
  double ene0, ene1;

  ene0 = 0.0;
  ene1 = model.H_I[i][a1];
  if(ilen == 0) {
    ene0 = model.J_OO[i][aaseq[i].O][aaseq[i+1].O];
    ene1 += model.J_OI[i][aaseq[i].O][a1] + model.J_IO[i][a1][aaseq[i+1].O];
  } else if(j == 0) {
    int a0 = aaseq[i].I[0];
    int b = aaseq[i].O;
    ene0 += model.J_OI[i][b][a0];
    ene1 += model.J_OI[i][b][a1] + model.J_II[i][a1][a0];
  } else if(j == ilen) {
    int a0 = aaseq[i].I[ilen-1];
    int b = aaseq[i+1].O;
    ene0 += model.J_IO[i][a0][b];
    ene1 += model.J_II[i][a0][a1] + model.J_IO[i][a1][b];
  } else {
    int a0 = aaseq[i].I[j];
    int b = aaseq[i].I[j-1];
    ene0 += model.J_II[i][b][a0];
    ene1 += model.J_II[i][b][a1] + model.J_II[i][a1][a0];
  }

  return (ene0 - ene1);
}

double MCseq::extend_insert_gs(const Model& model, const int i, const int j, 
			       int& a1)
{// This routine includes the case where ilen == 0. 
  const int ilen = aaseq[i].I.size();
  double e0=0.0;

  
  if(ilen == 0) {
    for(int a = 0; a < nstate_insert; ++a) {
      ep_insert[a] = model.H_I[i][a]
	+ model.J_OI[i][aaseq[i].O][a] + model.J_IO[i][a][aaseq[i+1].O];
    }
    e0 = model.J_OO[i][aaseq[i].O][aaseq[i+1].O];
  }
  else if(j == 0) {
    const int b = aaseq[i].O;
    const int a0 = aaseq[i].I[0];
    for(int a = 0; a < nstate_insert; ++a) {
      ep_insert[a] = model.H_I[i][a]
	+ model.J_OI[i][b][a] + model.J_II[i][a][a0];
    }
    e0 = model.J_OI[i][b][a0];
  } else if(j == ilen) {
    const int b = aaseq[i+1].O;
    const int a0 = aaseq[i].I[ilen-1];
    for(int a = 0; a < nstate_insert; ++a) {
      ep_insert[a] = model.H_I[i][a]
	+ model.J_II[i][a0][a] + model.J_IO[i][a][b];
    }
    e0 = model.J_IO[i][a0][b];
  } else {
    int a0 = aaseq[i].I[j];
    int b = aaseq[i].I[j-1];
    for(int a = 0; a < nstate_insert; ++a) {
      ep_insert[a] = model.H_I[i][a]
	+ model.J_II[i][b][a] + model.J_II[i][a][a0];
    }
    e0 = model.J_II[i][b][a0];
  }
  
  double w=0.0, t;
  for(int a = 0; a < nstate_insert; ++a) {
    t = exp(Beta*(ep_insert[a] - e0));
    w += t;
    dist_insert[a] = t;
  }

  for(int a = 0; a < nstate_insert; ++a) {
    dist_insert[a] /= w;
  }
  boost::random::discrete_distribution<> dist(dist_insert);
  a1 = dist(rng);
  return (e0 - ep_insert[a1]);
}

double MCseq::shorten_insert(const Model& model, const int i, const int j)
{
  const int iend = aaseq[i].I.size() - 1;
  int a0 = aaseq[i].I[j];
  
  double ene0 = model.H_I[i][a0];
  double ene1 = 0.0;
  if(iend == 0) { // ilen = 1; 1-residue insert
    ene0 += model.J_OI[i][aaseq[i].O][a0] + model.J_IO[i][a0][aaseq[i+1].O];
    ene1 += model.J_OO[i][aaseq[i].O][aaseq[i+1].O];
  } else if(j == 0) { // ilen > 1
    ene0 += model.J_OI[i][aaseq[i].O][a0] + model.J_II[i][a0][aaseq[i].I[j+1]];
    ene1 += model.J_OI[i][aaseq[i].O][aaseq[i].I[j+1]];
  } else if(j == iend) {
    ene0 += model.J_II[i][aaseq[i].I[j-1]][a0] + model.J_IO[i][a0][aaseq[i+1].O];
    ene1 += model.J_IO[i][aaseq[i].I[j-1]][aaseq[i+1].O];
  } else {
    ene0 += model.J_II[i][aaseq[i].I[j-1]][a0] + model.J_II[i][a0][aaseq[i].I[j+1]];
    ene1 += model.J_II[i][aaseq[i].I[j-1]][aaseq[i].I[j+1]];
  }

  return (ene0 - ene1);
}


double MCseq::sweep(const Model& model)
{
  int nacc = 0;
  double dene;
  boost::random::uniform_int_distribution<> randO(1,Length),randI(1,Length-1), imove(0,1);
  int i;
  for(int ii = 1; ii <= Length; ++ii) {
    i = randO(rng);
    const int a0 = aaseq[i].O;
    const int a1 = rand_core(rng);
    dene = mutate_core(model, i, a1);
    if(select_criteria(dene, 1.0)) {
      if(a0 != a1) nacc++;
      aaseq[i].O = a1;
      Energy += dene;
      if(a0 == terminal_state && a1 != a0) SeqLength++;
      else if(a1 == terminal_state && a1 != a0) SeqLength--;
    }

    if(ii == Length) continue;
    i = randI(rng);
    const int move = imove(rng);
    const int ilen = aaseq[i].I.size();
    const int iend = ilen - 1;
    const double qe = 1.0/20, qd = 1.0/21.0;
    const double p = (double)(ilen+1)/(2*ilen+1);
    if (dist01(rng) < p) { // open/extend
      boost::random::uniform_int_distribution<> rand_pos(0,ilen);
      const int j = rand_pos(rng);
      const int a1 = rand_insert(rng);
      dene = extend_insert(model, i, j, a1);
      const double pp = (double)(ilen+1)/(2*ilen+3);
      const double pins = pp*qd/(p*qe);
      if(select_criteria(dene, pins)) {
	aaseq[i].I.insert(aaseq[i].I.begin()+j, a1);
	Energy += dene;
	SeqLength++;
	nacc++;
      }
    } else if(dist01(rng) < qd) {
      boost::random::uniform_int_distribution<> rand_pos(0,iend);
      const int j = rand_pos(rng);
      dene = shorten_insert(model, i, j);
      const double pm = ilen/(2*ilen - 1);
      const double pdel = pm*qe/((1.0-p)*qd);
      if(select_criteria(dene, pdel)) {
	aaseq[i].I.erase(aaseq[i].I.begin() + j);
	Energy += dene;
	SeqLength--;
	nacc++;
      }
    } else {
      boost::random::uniform_int_distribution<> rand_pos(0,iend);
      const int j = rand_pos(rng);
      const int a1 = rand_insert(rng);
      dene = mutate_insert(model, i, j, a1);
      if(select_criteria(dene, 1.0)) {
	if(aaseq[i].I[j] != a1) nacc++;
	aaseq[i].I[j] = a1;
	Energy += dene;
      }
    }
  }
  return nacc;
}

double MCseq::sweep_gs(const Model& model)
{
  int nacc = 0;
  double dene;
  boost::random::uniform_int_distribution<> randO(1,Length),randI(1,Length-1);
  int i;
  int a1;
  for(int ii = 1; ii <= Length; ++ii) {
    i = randO(rng);
    const int a0 = aaseq[i].O;
    dene = mutate_core_gs(model, i, a1);
    if(a0 != a1) nacc++;
    if(a0 == terminal_state && a0 != a1) SeqLength++;
    else if(a1 == terminal_state && a0 != a1) SeqLength--;
    aaseq[i].O = a1;
    Energy += dene;

    if(ii == Length) continue;
    i = randI(rng);
    const int ilen = aaseq[i].I.size();
    const int iend = ilen - 1;
    const double qe = 1.0/20, qd = 1.0/21.0;
    const double p = (double)(ilen+1)/(2*ilen+1);
    if (dist01(rng) < p) { // open/extend
      boost::random::uniform_int_distribution<> rand_pos(0,ilen);
      const int j = rand_pos(rng);
      const int a1 = rand_insert(rng);
      dene = extend_insert(model, i, j, a1);
      const double pp = (double)(ilen+1)/(2*ilen+3);
      const double pins = pp*qd/(p*qe);
      if(select_criteria(dene,pins)) {
	aaseq[i].I.insert(aaseq[i].I.begin()+j, a1);
	Energy += dene;
	SeqLength++;
	nacc++;
      }
    } else if(dist01(rng) < qd) { 
      boost::random::uniform_int_distribution<> rand_pos(0,iend);
      const int j = rand_pos(rng);
      dene = shorten_insert(model, i, j);
      const double pm = ilen/(2*ilen - 1);
      const double pdel = pm*qe/((1.0-p)*qd);
      if(select_criteria(dene,pdel)) {
	aaseq[i].I.erase(aaseq[i].I.begin() + j);
	Energy += dene;
	SeqLength--;
	nacc++;
      }
    } else {
      boost::random::uniform_int_distribution<> rand_pos(0,iend);
      const int j = rand_pos(rng);
      dene = mutate_insert_gs(model, i, j, a1);
      if(aaseq[i].I[j] != a1) nacc++;
      aaseq[i].I[j] = a1;
      Energy += dene;
    }
  }

  return nacc;
}

void MCseq::scan_natural_sequences(Model& model)
{
  const int nseq = model.get_num_seq();
  emin_natural = 1e+100;
  emax_natural = -1e+100;
  imin_natural = imax_natural = -1;
  eave_natural = evar_natural = 0;
  for(int i = 0; i < nseq; ++i) {
    double w = model.get_weight(i);
    set_natural_sequence(model,i);
    //    cerr << "SCAN_NATURAL " << i << "\t" << Energy << endl;
    if(Energy < emin_natural) {
      emin_natural = Energy;
      imin_natural = i;
    }
    if(Energy > emax_natural) {
      emax_natural = Energy;
      imax_natural = i;
    }
    eave_natural += w*Energy;
    evar_natural += w*Energy*Energy;
  }

  evar_natural -= eave_natural*eave_natural;
}

int MCseq::run(Model& model, const int nsweep, const bool samplep)
{
  double nacc_total=0;
  int isweep;
  int ncount=0;
  if(samplep) {
    WLfactor = 0;
  }
  // en_{O|I}, ens_{OO|OI|IO|II}, enl_OO are supposed
  // to be initialized appropriately.
  get_energy(model);
  for(isweep = 1; isweep <= nsweep; ++isweep) {
    const double nacc = (select_type == CANONICAL && gibbs_sampler_p) ?
      sweep_gs(model) : sweep(model);
    nacc_total += nacc;
    if(select_type == MULTI_CANONICAL && WLfactor > 0) {
      if(check_EHist() && WLfactor >= WLprec)
	clear_EHist();
    }
    if(select_type == MULTI_CANONICAL && WLfactor < WLprec && !samplep)
      return isweep;

    if(MONITOR_p && (isweep % NMONITOR == 0))
      cerr << "MCseq::run i,acc,E,slen,div\t"
	   << isweep << "\t" 
	   << 100.0*nacc/(2*Length+1) << "\t"
	   << Energy << "\t"
	   << SeqLength << "\t"
	   << get_divergence(model) << endl;

    if(samplep) {
      count_stats(model);
      ncount++;
    }
    if(trajectory_p && (isweep % NTRAJ == 0)) {
      save_trajectory();
    }
  }
  isweep--;
  const double mutation_rate = 100*nacc_total/((2*Length+1)*isweep);
  get_energy(model);
  if(MONITOR_p) cerr << "# T,Energy,mut_rate(Final)\t"
		     << Temperature << "\t"
		     << Energy << "\t"
		     << mutation_rate << endl;

  return isweep;
}

double MCseq::count_stats(Model& model, const bool just_weight)
{
  double w = 1.0;
  if(STATS_MULTICANONICAL) {
    const int iene = get_DOS_bindex(Energy);
    if(DOS.count(iene) == 0) DOS[iene]=1.0;
    w =  exp(-(Energy - MIN0_ENERGY)*Beta + (DOS[iene] - MAX_DOS));
  }
  WeightTotal += w;
  const double we=w*Energy, wl=w*SeqLength;
  Eave += we;
  Evar += we*Energy;
  Lave += wl;
  Lvar += wl*SeqLength;
  
  if(just_weight) return w;
  /////////////////////////////////////////////
  for(int i = 1; i <= Length; ++i) {
    const int a = aaseq[i].O;
    model.en_O[i][a] += w;
    model.nE_O[i][a] += w*Energy;
    const int ilen = aaseq[i].I.size();
    if(ilen > 0) {
      for(int j = 0; j < ilen; ++j) {
    	model.en_I[i][aaseq[i].I[j]] += w;
	model.nE_I[i][aaseq[i].I[j]] += w*Energy;
      }
      model.ens_OI[i][a][aaseq[i].I[0]] += w;
      model.ens_IO[i][aaseq[i].I[ilen-1]][aaseq[i+1].O] += w;
      for(int j = 0; j <= ilen - 2; ++j) {
    	model.ens_II[i][aaseq[i].I[j]][aaseq[i].I[j+1]] += w;
      }
    } else {
      model.ens_OO[i][a][aaseq[i+1].O] += w;
    }
    
    if(stats_full_nonbonded) {
      vector<int> insa = aaseq[i].I;
      for(int j = i; j < Length; ++j) {
	vector<int> insb = aaseq[j].I;
	for(vector<int>::iterator a = insa.begin(); a != insa.end(); ++a) {
	  for(vector<int>::iterator b = insb.begin(); b != insb.end(); ++b) {
	    model.enl_II[i][j][*a][*b] += w;
	  }
	}
      }
      model.enl_OO[i][i][a][a] += w;
      model.nE_OO[i][i][a][a] += w*Energy;
      for(int j = i; j <= Length; ++j) {
	model.enl_OO[i][j][a][aaseq[j].O] += w;
	model.nE_OO[i][j][a][aaseq[j].O] += w*Energy;
      }
      for(int j = 1; j < Length; ++j) {
	vector<int> insb = aaseq[j].I;
	for(vector<int>::iterator b = insb.begin(); b != insb.end(); ++b) {
	  model.enl_OI[i][j][a][*b] += w;
	}
      }
    }
  }

  if(!stats_full_nonbonded) {
    IntPairs &intpairs = model.get_intpairs();
    for(IntPairs::iterator it = intpairs.begin(); it != intpairs.end(); ++it) {
      const int i = it->i, j = it->j;
      model.enl_OO[i][j][aaseq[i].O][aaseq[j].O] += w;
      model.nE_OO[i][j][aaseq[i].O][aaseq[j].O] += w*Energy;
    }
  }

  return w;
}

void MCseq::save_DOS(const string filename)
{
  ofstream fout(filename.c_str());

  double dmin = 1e+100;
  for(map<int,double>::iterator it = DOS.begin(); it != DOS.end(); ++it) {
    if(it->second < dmin) dmin = it->second;
  }

  dmin -= 1.0;
  fout << "#MIN_ENERGY,MAX_ENERGY,DELTA_ENERGY,WLfactor= "
       << MIN_ENERGY << "\t"
       << MAX_ENERGY << "\t"
       << DELTA_ENERGY << "\t"
       << WLfactor << endl;
  fout << "#NBINS= " << DOS.size() << endl;
  for(map<int,double>::iterator it = DOS.begin(); it != DOS.end(); ++it) {
    fout << it->first << "\t"
	 << it->first * DELTA_ENERGY << "\t"
	 << it->second - dmin << "\t" 
	 << EHist[it->first] << endl;
  }

  fout.close();
  return;
}

void MCseq::setup_WL(const double min_energy, const double max_energy,
		     const double delta_energy) {
  MIN_ENERGY = min_energy;
  MAX_ENERGY = max_energy;
  DELTA_ENERGY = delta_energy;
  MAX_DOS = 0;
  MIN0_ENERGY = 1e+10;
  clear_DOS();
  clear_EHist();
  WLfactor = 1.0;
  cerr << "MCseq::setup_WL: MIN_ENERGY,MAX_ENERGY,DELTA_ENERGY= " << "\t"
       << MIN_ENERGY << "\t"
       << MAX_ENERGY << "\t"
       << DELTA_ENERGY << endl;
  return;
}

int MCseq::read_DOS(const string filename, const bool reset_max)
{
  ifstream fin(filename.c_str());
  string header;
  int nbins,bind;
  double ene,count,min_energy,max_energy,ehist;
  MAX_DOS = -1;
  MIN0_ENERGY = 1e+100;
  clear_DOS();

  fin >> header >> min_energy >> max_energy >> DELTA_ENERGY >> WLfactor;
  if(reset_max) {
    MIN_ENERGY = min_energy;
    MAX_ENERGY = max_energy;
    cerr << "MCseq::read_DOS: resetting min/max/delta to " << min_energy
	 << "\t" << max_energy
	 << "\t" << DELTA_ENERGY << endl;
  }
  fin >> header >> nbins;
  for(int i=0; i < nbins; ++i) {
    fin >> bind >> ene >> count >> ehist;
    if(!reset_max && ene >= MAX_ENERGY) continue;
    if(!reset_max && ene < MIN_ENERGY) continue;
    if(count > MAX_DOS) MAX_DOS = count;
    if(ene < MIN0_ENERGY) MIN0_ENERGY = ene;
    DOS[bind] = count;
  }
  fin.close();
  clear_EHist();
  return nbins;
}

int MCseq::anneal(Model& model, const int nsweeps, const int nsteps,
		  const double temp_high, const double temp_low)
{
  double dt = (temp_high - temp_low)/nsteps;
  for (int n = 0; n <= nsteps; ++n) {
    const double tt = temp_high - dt*n;
    cerr << "# Annealing: T = " << tt << endl;
    set_temperature(tt);
    run(model, nsweeps, false);
  }

  return nsweeps*nsteps;
}

