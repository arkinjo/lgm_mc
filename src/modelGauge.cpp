#include "model.hpp"

void Model::set_Ising_gauge_K(void)
{
  darray1 Kooa(boost::extents[nstate_core]),Koob(boost::extents[nstate_core]);

  for(int i = 1; i <= model_length - 2; ++i) {
    for(int j = i + 2; j <= model_length; ++j) {
      std::fill(Kooa.data(), Kooa.data() + Kooa.num_elements(),0.0);
      std::fill(Koob.data(), Koob.data() + Koob.num_elements(),0.0);
      double koo=0.0;
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  Kooa[a] += K_OO[i][j][a][b];
	  Koob[b] += K_OO[i][j][a][b];
	  koo += K_OO[i][j][a][b];
	}
      }
      koo /= nstate_core*nstate_core;
      for(int a = 0; a < nstate_core; ++a) {
	Kooa[a] /= nstate_core;
	Koob[a] /= nstate_core;
      }
      for(int a = 0; a < nstate_core; ++a) {
	H_O[i][a] += Kooa[a];
	H_O[j][a] += Koob[a];
	for(int b = 0; b < nstate_core; ++b) {
	  K_OO[i][j][a][b] += koo - Kooa[a] - Koob[b];
	  K_OO[j][i][b][a] = K_OO[i][j][a][b];
	}
      }
    }
  }
  return;
}

double Model::zero_Jp(const int i)
{
  darray1 JOp(boost::extents[nstate_core]),JIp(boost::extents[nstate_insert]);
  double del = 0.0;
  
  set_JOp(i, JOp);
  set_JIp(i, JIp);

  if(i == 0) {
    const int a = terminal_state;
    const double sca = 1.0/(nstate_core + nstate_insert);

    del += fabs(JOp[a]);
    for(int b = 0; b < nstate_core; ++b) {
      J_OO[i][a][b] -= sca*JOp[a];
    }
    for(int b = 0; b < nstate_insert; ++b) {
      J_OI[i][a][b] -= sca*JOp[a];
    }
    
    for(int a = 0; a < nstate_insert; ++a) {
      del += fabs(JIp[a]);
      H_I[i][a] += JIp[a];
      for(int b = 0; b < nstate_core; ++b) {
	J_IO[i][a][b] -= sca*JIp[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i][a][b] -= sca*JOp[a];
      }
    }	
  } else if (i == model_length) {
    const int b = terminal_state;
    const double sca = 1.0/(1.0 + nstate_insert);
    for(int a = 0; a < nstate_core; ++a) {
      del += fabs(JOp[a]);
      H_O[i][a] += JOp[a];
      J_OO[i][a][b] -= sca*JOp[a];
      for(int b = 0; b < nstate_insert; ++b) {
	J_OI[i][a][b] -= sca*JOp[a];
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      del += fabs(JIp[a]);
      H_I[i][a] += JIp[a];
      J_IO[i][a][b] -= sca*JIp[a];
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i][a][b] -= sca*JIp[a];
      }
    }
  } else {
    const double sca = 1.0/(nstate_core + nstate_insert);
    for(int a = 0; a < nstate_core; ++a) {
      del += fabs(JOp[a]);
      H_O[i][a] += JOp[a];
      for(int b = 0; b < nstate_core; ++b) {
	J_OO[i][a][b] -= sca*JOp[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_OI[i][a][b] -= sca*JOp[a];
      }
    }
    for(int a = 0; a < nstate_insert; ++a) {
      del += fabs(JIp[a]);
      H_I[i][a] += JIp[a];
      for(int b = 0; b < nstate_core; ++b) {
	J_IO[i][a][b] -= sca*JIp[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i][a][b] -= sca*JIp[a];
      }
    }
  }

  return del;
}

double Model::zero_Jm(const int i)
{
  darray1 JO(boost::extents[nstate_core]),JI(boost::extents[nstate_insert]);
  double del = 0.0;
  
  set_JOm(i, JO);
  set_JIm(i-1, JI);

  if(i == model_length + 1) {
    const int a = terminal_state;
    const double sca = 1.0/(nstate_core + nstate_insert);
    del += fabs(JO[a]);
    for(int b = 0; b < nstate_core; ++b) {
      J_OO[i-1][b][a] -= sca*JO[a];
    }
    for(int b = 0; b < nstate_insert; ++b) {
      J_IO[i-1][b][a] -= sca*JO[a];
    }

    for(int a = 0; a < nstate_insert; ++a) {
      del += fabs(JI[a]);
      H_I[i-1][a] += JI[a];
      for(int b = 0; b < nstate_core; ++b) {
	J_OI[i-1][b][a] -= sca*JI[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i-1][b][a] -= sca*JI[a];
      }
    }
  } if (i == 1) {
    const int b = terminal_state;
    const double sca = 1.0/(1 + nstate_insert);
    for(int a = 0; a < nstate_core; ++a) {
      del += fabs(JO[a]);
      H_O[i][a] += JO[a];
      J_OO[i-1][b][a] -= sca*JO[a];
    }
    for(int a = 0; a < nstate_insert; ++a) {
      H_I[i-1][a] += JI[a];
      J_OI[i-1][b][a] -= sca*JI[a];
    }

    for(int a = 0; a < nstate_insert; ++a) {
      del += fabs(JI[a]);
      H_I[i-1][a] += JI[a];
      J_OI[i-1][b][a] -= sca*JI[a];
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i-1][b][a] -= sca*JI[a];
      }
    }
  } else {
    const double sca = 1.0/(nstate_core + nstate_insert);
    for(int a = 0; a < nstate_core; ++a) {
      del += fabs(JO[a]);
      H_O[i][a] += JO[a];
      for(int b = 0; b < nstate_core; ++b) {
	J_OO[i-1][b][a] -= sca*JO[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_IO[i-1][b][a] -= sca*JO[a];
      }
    }
    
    for(int a = 0; a < nstate_insert; ++a) {
      del += fabs(JI[a]);
      H_I[i-1][a] += JI[a];
      for(int b = 0; b < nstate_core; ++b) {
	J_OI[i-1][b][a] -= sca*JI[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i-1][b][a] -= sca*JI[a];
      }
    }
  }

  return del;
}

void Model::set_Ising_gauge_J(void)
{
  darray1 JOm(boost::extents[nstate_core]),JOp(boost::extents[nstate_core]);
  darray1 JIm(boost::extents[nstate_core]),JIp(boost::extents[nstate_core]);


  for(int i = 0; i <= model_length; ++i) {
    set_JOp(i,JOp);
    set_JIp(i, JIp);
    set_JIm(i, JIm);
    set_JOm(i+1, JOm);
    double com=0.0;
    if(i == 0) {
      const int ao = terminal_state;
      double Qp = nstate_core + nstate_insert;
      double Qm = 1 + nstate_insert;
      com = 0.0;
      /*
      for(int bo = 0; bo < nstate_core; ++bo) {
	com += JOm[bo];
      }
      for(int bi = 0; bi < nstate_insert; ++bi) {
	com += JIm[bi];
      }
      */
      com = JOp[ao];
      for(int ai = 0; ai < nstate_insert; ++ai) {
	com += JIp[ai];
      }
      com /= (Qp*Qm);
      for(int bo = 0; bo < nstate_core; ++bo) {
	H_O[i+1][bo] += JOm[bo]/Qm;
	J_OO[i][ao][bo] += -JOp[ao]/Qp -JOm[bo]/Qm + com;
      }
      for(int bi = 0; bi < nstate_insert; ++bi) {
	J_OI[i][ao][bi] += -JOp[ao]/Qp -JIm[bi]/Qm + com;
      }
      for(int ai = 0; ai < nstate_insert; ++ai) {
	H_I[i][ai] += JIp[ai]/Qp + JIm[ai]/Qm - com;
	for(int bo = 0; bo < nstate_core; ++bo) {
	  J_IO[i][ai][bo] += -JIp[ai]/Qp -JOm[bo]/Qm + com;
	}
	for(int bi = 0; bi < nstate_insert; ++bi) {
	  J_II[i][ai][bi] += -JIp[ai]/Qp -JIm[bi]/Qm + com;
	}
      }
    }
    else if(i == model_length) {
      const int bo = terminal_state;
      double Qp = 1 + nstate_insert;
      double Qm = nstate_core + nstate_insert;
      com += JOm[bo];
      for(int bi = 0; bi < nstate_insert; ++bi) {
	com += JIm[bi];
      }
      com /= (Qp*Qm);

      for(int ao = 0; ao < nstate_core; ++ao) {
	H_O[i][ao] += JOp[ao]/Qp ;
	J_OO[i][ao][bo] += -JOp[ao]/Qp - JOm[bo]/Qm + com;
	for(int bi = 0; bi < nstate_insert; ++bi) {
	  J_OI[i][ao][bi] += -JOp[ao]/Qp - JIm[bi]/Qm + com;
	}
      }
      for(int ai = 0; ai < nstate_insert; ++ai) {
	H_I[i][ai] += JIp[ai]/Qp + JIm[ai]/Qm - com;
	J_IO[i][ai][bo] += -JIp[ai]/Qp - JOm[bo]/Qm + com;
	for(int bi = 0; bi < nstate_insert; ++bi) {
	  J_II[i][ai][bi] += -JIp[ai]/Qp - JIm[bi]/Qm + com;
	}
      }
    }
    else {
      double Qp = nstate_core + nstate_insert;
      double Qm = nstate_core + nstate_insert;
      for(int bo = 0; bo < nstate_core; ++bo) {
	com += JOm[bo];
      }
      for(int bi = 0; bi < nstate_insert; ++bi) {
	com += JIm[bi];
      }
      com /= (Qp*Qm);

      for(int ao = 0; ao < nstate_core; ++ao) {
	H_O[i][ao] += JOp[ao]/Qp ;
	H_O[i+1][ao] += JOm[ao]/Qm ;
	for(int bo = 0; bo < nstate_core; ++bo) {
	  J_OO[i][ao][bo] += -JOp[ao]/Qp - JOm[bo]/Qm + com;
	}
	for(int bi = 0; bi < nstate_insert; ++bi) {
	  J_OI[i][ao][bi] += -JOp[ao]/Qp - JIm[bi]/Qm + com;
	}
      }
      for(int ai = 0; ai < nstate_insert; ++ai) {
	H_I[i][ai] += JIp[ai]/Qp + JIm[ai]/Qm - com;
	for(int bo = 0; bo < nstate_core; ++bo) {
	  J_IO[i][ai][bo] += -JIp[ai]/Qp - JOm[bo]/Qm + com;
	}
	for(int bi = 0; bi < nstate_insert; ++bi) {
	  J_II[i][ai][bi] += -JIp[ai]/Qp - JIm[bi]/Qm + com;
	}
      }
    }
  }

  return;
}

void Model::set_Ising_gauge_H(void)
{
  for(int i = 1; i <= model_length; ++i) {
    double aO = 0.0;
    for(int a = 0; a < nstate_core; ++a) {
      aO += H_O[i][a];
    }
    aO /= nstate_core;
    for(int a = 0; a < nstate_core; ++a) {
      H_O[i][a] -= aO;
    }
  }
  
  return;
}

void Model::eliminate_H(void)
{
  for(int i = 0; i <= model_length; ++i) {
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i][a][b] += H_I[i][b];
      }
      for(int b = 0; b < nstate_core; ++b) {
	if(i==model_length && b != terminal_state) continue;
	J_IO[i][a][b] += H_O[i+1][b];
      }
    }
    for(int a = 0; a < nstate_core; ++a) {
      if(i == 0 && a != terminal_state) continue;
      for(int b = 0; b < nstate_core; ++b) {
	if(i==model_length && b != terminal_state) continue;
	J_OO[i][a][b] += H_O[i+1][b];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_OI[i][a][b] += H_I[i][b];
      }
    }
  }

  std::fill(H_O.data(), H_O.data() + H_O.num_elements(), 0.0);
  std::fill(H_I.data(), H_I.data() + H_I.num_elements(), 0.0);
}

void Model::set_JOp(const int i, darray1& jv)
{
  std::fill(jv.data(), jv.data() + jv.num_elements(),0.0);

  if(i > model_length) return;

  if(i == 0) {
    const int a = terminal_state;
    for(int b = 0; b < nstate_core; ++b) {
      jv[a] += J_OO[i][a][b];
    }
    for(int b = 0; b < nstate_insert; ++b) {
      jv[a] += J_OI[i][a][b];
    }
  } else if (i == model_length)  {
    const int b = terminal_state;
    for(int a = 0; a < nstate_core; ++a) {
      jv[a] += J_OO[i][a][b];
    }
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	jv[a] += J_OI[i][a][b];
      }
    }
  } else {
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	jv[a] += J_OO[i][a][b];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	jv[a] += J_OI[i][a][b];
      }
    }
  }

  return;
}

void Model::set_JOm(const int i, darray1& jv)
{
  std::fill(jv.data(), jv.data() + jv.num_elements(),0.0);

  if(i == 0) return;

  if(i == model_length+1) {
    const int a = terminal_state;
    for(int b = 0; b < nstate_core; ++b) {
      jv[a] += J_OO[i-1][b][a];
    }
    for(int b = 0; b < nstate_insert; ++b) {
      jv[a] += J_IO[i-1][b][a];
    }
  } else if (i == 1) {
    const int b = terminal_state;
    for(int a = 0; a < nstate_core; ++a) {
      jv[a] += J_OO[i-1][b][a];
    }
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	jv[a] += J_IO[i-1][b][a];
      }
    }
  } else {
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	jv[a] += J_OO[i-1][b][a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	jv[a] += J_IO[i-1][b][a];
      }
    }
  }

  return;
}

void Model::set_JIp(const int i, darray1& jv)
{
  std::fill(jv.data(), jv.data() + jv.num_elements(),0.0);

  if(i > model_length) return;

  for(int a = 0; a < nstate_insert; ++a) {
    for(int b = 0; b < nstate_insert; ++b) {
      jv[a] += J_II[i][a][b];
    }
  }

  if(i == model_length) {
    const int b = terminal_state;
    for(int a = 0; a < nstate_insert; ++a) {
      jv[a] += J_IO[i][a][b];
    }
  } else {
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	jv[a] += J_IO[i][a][b];
      }
    }
  }
  

  return;
}

void Model::set_JIm(const int i, darray1& jv)
{
  std::fill(jv.data(), jv.data() + jv.num_elements(),0.0);

  if (i > model_length) return;
  
  if(i == 0) {
    const int b = terminal_state;
    for(int a = 0; a < nstate_insert; ++a) {
      jv[a] += J_OI[i][b][a];
    }
  } else {
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	jv[a] += J_OI[i][b][a];
      }
    }
  }
  for(int a = 0; a < nstate_insert; ++a) {
    for(int b = 0; b < nstate_insert; ++b) {
      jv[a] += J_II[i][b][a];
    }
  }

  return;
}

double Model::symmetrize_JI(const int i)
{
  darray1 JIm(boost::extents[nstate_insert]),JIp(boost::extents[nstate_insert]),
    csI(boost::extents[nstate_insert]);

  const int qo0 = (i == 0 || i == model_length + 1) ? 1 : nstate_core;
  const int qop = (i == model_length) ? 1 : nstate_core;
  const int qom = (i == 1) ? 1 : nstate_core;
  const int qi = nstate_insert;

  double delta = 0.0;
  
  set_JIp(i, JIp);
  set_JIm(i, JIm);

  std::fill(csI.data(), csI.data() + csI.num_elements(),0.0);
  double c = 0.0;
  for(int a = 0; a < nstate_insert; ++a) {
    double d = JIp[a] - JIm[a];
    c += d;
    delta += fabs(d);
  }
  c *= 2.0/(qo0 + qop);
  for(int a = 0; a < nstate_insert; ++a) {
    csI[a] = (JIp[a] - JIm[a] + c)/(qo0 + qop + 2*qi);
  }

  if(i == 0) {
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	J_IO[i][a][b] -= csI[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i][a][b] += csI[b] - csI[a];
      }
      const int b = terminal_state;
      J_OI[i][b][a] += csI[a];      
    }
  } else if (i == model_length) {
    for(int a = 0; a < nstate_insert; ++a) {
      const int b = terminal_state;
      J_IO[i][a][b] -= csI[a];
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i][a][b] += csI[b] - csI[a];
      }
      for(int b = 0; b < nstate_core; ++b) {
	J_OI[i][b][a] += csI[a];
      }
    }
  } else {
    for(int a = 0; a < nstate_insert; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	J_IO[i][a][b] -= csI[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_II[i][a][b] += csI[b] - csI[a];
      }
      for(int b = 0; b < nstate_core; ++b) {
	J_OI[i][b][a] += csI[a];
      }
    }
  }

  return delta;
}

double Model::symmetrize_JO(const int i)
{
  darray1 JOm(boost::extents[nstate_core]),JOp(boost::extents[nstate_core]),
    csO(boost::extents[nstate_core]);

  double delta = 0.0;
  
  const int qo0 = (i == 0 || i == model_length + 1) ? 1 : nstate_core;
  const int qop = (i == model_length) ? 1 : nstate_core;
  const int qom = (i == 1) ? 1 : nstate_core;
  const int qi = nstate_insert;
 
  set_JOp(i, JOp);
  set_JOm(i, JOm);
  std::fill(csO.data(), csO.data() + csO.num_elements(),0.0);
  for(int a = 0; a < nstate_core; a++) {
    double d = (JOp[a] - JOm[a]);
    csO[a] = d/(qom + qop + 2*qi);
    delta += fabs(d);
  }
  if(i == 1) {
    for(int a = 0; a < nstate_core; a++) {
      for(int b = 0; b < nstate_core; b++) {
	J_OO[i][a][b] -= csO[a];
      }
      for(int b = 0; b < nstate_insert; b++) {
	J_OI[i][a][b] -= csO[a];
      }
      const int b = terminal_state;
      J_OO[i-1][b][a] += csO[a];
      for(int b = 0; b < nstate_insert; ++b) {
	J_IO[i-1][b][a] += csO[a];
      }
    }
  } else if (i == model_length) {
    const int b = terminal_state;
    for(int a = 0; a < nstate_core; a++) {
      J_OO[i][a][b] -= csO[a];
      for(int b = 0; b < nstate_insert; b++) {
	J_OI[i][a][b] -= csO[a];
      }
      for(int b = 0; b < nstate_core; ++b) {
	J_OO[i-1][b][a] += csO[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_IO[i-1][b][a] += csO[a];
      }
    }
  } else {
    for(int a = 0; a < nstate_core; a++) {
      for(int b = 0; b < nstate_core; ++b) {
	J_OO[i][a][b] -= csO[a];
      }
      for(int b = 0; b < nstate_insert; b++) {
	J_OI[i][a][b] -= csO[a];
      }
      for(int b = 0; b < nstate_core; ++b) {
	J_OO[i-1][b][a] += csO[a];
      }
      for(int b = 0; b < nstate_insert; ++b) {
	J_IO[i-1][b][a] += csO[a];
      }
    }
  }

  return delta;
}

double Model::normalize_JOp(const int i)
{
  darray1 JOp(boost::extents[nstate_core]);
  const int qo0 = (i == 0 || i == model_length + 1) ? 1 : nstate_core;
  const int qop = (i == model_length) ? 1 : nstate_core;
  const int qi = nstate_insert;

  set_JOp(i, JOp);

  double c=0.0;
  if(i == 0) {
    const int a = terminal_state;
    c = JOp[a];
    c /= qo0*(qop + qi);
    for(int b = 0; b < nstate_core; ++b) {
      J_OO[i][a][b] -= c;
    }
    for(int b = 0; b < nstate_insert; ++b) {
      J_OI[i][a][b] -= c;
    }
  } else if (i == model_length) {
    for(int a = 0; a < nstate_core; ++a) {
      c += JOp[a];
    }
    c /= qo0*(qop + qi);
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	J_OI[i][a][b] -= c;
	J_IO[i-1][b][a] += c;
      }
    }
    for(int a = 0; a < nstate_core; ++a) {
      J_OO[i][a][terminal_state] -= c;
      for(int b = 0; b < nstate_core; ++b) {
	J_OO[i-1][b][a] += c;
      }
    }
  } else {
    for(int a = 0; a < nstate_core; ++a) {
      c += JOp[a];
    }
    c /= qo0*(qop + qi);
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	J_OI[i][a][b] -= c;
	J_IO[i-1][b][a] += c;
      }
    }
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_core; ++b) {
	J_OO[i][a][b] -= c;
      }
    }
    if(i == 1) {
      for(int a = 0; a < nstate_core; ++a) {
	J_OO[i-1][terminal_state][a] += c;
      }
    } else {
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  J_OO[i-1][b][a] += c;
	}
      }
    }
  }

  return (c);
}

double Model::normalize_JOm(const int i)
{
  darray1 JOm(boost::extents[nstate_core]);
  const int qo0 = (i == model_length + 1) ? 1 : nstate_core;
  const int qom = (i == 1) ? 1 : nstate_core;
  const int qi = nstate_insert;


  set_JOm(i, JOm);
  
  double c=0.0;

  if(i == model_length+1) {
    const int a = terminal_state;
    c = JOm[a]/(nstate_core + nstate_insert);
    for(int b = 0; b < nstate_core; ++b) {
      J_OO[i-1][b][a] -= c;
    }
    for(int b = 0; b < nstate_insert; ++b) {
      J_IO[i-1][b][a] -= c;
    }
  } else if (i == 1) {
    for(int a = 0; a < nstate_core; ++a) {
      c += JOm[a];
    }
    c /= qo0*(qom + qi);
    const int b = terminal_state;
    for(int a = 0; a < nstate_core; ++a) {
      J_OO[i-1][b][a] -= c;
      for(int b = 0; b < nstate_core; ++b) {
	J_OO[i][a][b] += c;
      }
    }
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	J_IO[i-1][b][a] -= c;
	J_OI[i][a][b] += c;
      }
    }
  } else {
    for(int a = 0; a < nstate_core; ++a) {
      c += JOm[a];
    }
    c /= qo0*(qom + qi);
    if (i == model_length) {
      for(int a = 0; a < nstate_core; ++a) {
	J_OO[i][a][terminal_state] += c;
      }
    } else {
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  J_OO[i][a][b] += c;
	}
      }
    }

    for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  J_OO[i-1][b][a] -= c;
	}
    }
    for(int a = 0; a < nstate_core; ++a) {
      for(int b = 0; b < nstate_insert; ++b) {
	J_IO[i-1][b][a] -= c;
	J_OI[i][a][b] += c;
      }
    }
  }

  return fabs(c);
}

double Model::update_J_gauge(void)
{
  double Delta=0.0;
  // Forward
  for(int i = 0; i <= model_length; ++i) {
    Delta += normalize_JOp(i);
    Delta += symmetrize_JI(i);
    if(i > 0) Delta += symmetrize_JO(i);
  }

  double DeltaB = 0.0;
  DeltaB += normalize_JOm(model_length+1);
  // Backward
  for(int i = model_length+1; i >= 1; --i) {
    DeltaB += normalize_JOm(i);
    if(i <= model_length) {
      DeltaB += symmetrize_JI(i);
      DeltaB += symmetrize_JO(i);
    }
  }
  return (Delta + DeltaB);
}

double Model::fix_B_gauge_J(void)
{
  const int max_iterJ = 10000;
  int iter;
  double deltaJ,dold;
  const double eps = 1.0e-9;
  
  cout.setf(std::ios::scientific);
  deltaJ = 0;
  for(iter = 1; iter < max_iterJ; ++iter) {
    dold = deltaJ;
    deltaJ = update_J_gauge();
    if(deltaJ < eps) break;
    //    if(iter % 100 == 0)
    //      cout << "delJ " << iter << "\t" << deltaJ << "\t" << dold << endl;
  }

  //  cout << "delJ " << iter << "\t" << deltaJ << endl;
  return deltaJ;
}


double Model::fix_B_gauge_K(void)
{
  darray1 Kooa(boost::extents[nstate_core]),Koob(boost::extents[nstate_core]);
  darray1 JOi(boost::extents[nstate_core]),JOj(boost::extents[nstate_core]);

  const double Q = 2*nstate_core + nstate_insert;
  
  double DeltaOO=0.0;
  
  for(int i = 1; i <= model_length - 2; ++i) {
    for(int j = i + 2; j <= model_length; ++j) {
      normalize_JOp(i); normalize_JOm(j);
      set_JOp(i, JOi);   set_JOm(j, JOj);
      std::fill(Kooa.data(), Kooa.data() + Kooa.num_elements(),0.0);
      std::fill(Koob.data(), Koob.data() + Koob.num_elements(),0.0);
      double koo=0.0;
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  Kooa[a] += K_OO[i][j][a][b];
	  Koob[b] += K_OO[i][j][a][b];
	  koo += K_OO[i][j][a][b];
	}
      }
      koo *= 1.0/(nstate_core*Q);
      for(int a = 0; a < nstate_core; ++a) {
	Kooa[a] = (Kooa[a] - JOi[a])/Q;
	Koob[a] = (Koob[a] - JOj[a])/Q;
      }
      for(int a = 0; a < nstate_core; ++a) {
	for(int b = 0; b < nstate_core; ++b) {
	  K_OO[i][j][a][b] += koo - Kooa[a] - Koob[b];
	  K_OO[j][i][b][a] = K_OO[i][j][a][b];

	  J_OO[i][a][b] += Kooa[a];
	  J_OO[j-1][a][b] += Koob[b];
	}
	for(int b = 0; b < nstate_insert; ++b) {
	  J_OI[i][a][b] += Kooa[a];
	  J_IO[j-1][b][a] += Koob[a];
	}
	DeltaOO += fabs(Kooa[a]) + fabs(Koob[a]);;
      }
      //      fix_B_gauge_J();
    }
  }

  return DeltaOO;
}


void Model::fix_B_gauge(void)
{
  double deltaJ = 0.0;
  const int max_iter = 5000;
  const double eps = 1.0e-9;

  cerr << "Fixing to B-gauge..." << endl;
  eliminate_H();
  for(int iter = 1; iter <= max_iter; ++iter) {
    deltaJ = fix_B_gauge_J();
    double deltaK = fix_B_gauge_K();
    if(iter % 100 == 0)
      cout << "dev(J,K) = " << iter << "\t" << deltaJ << "\t" << deltaK << endl;
    if(deltaK < eps) break;
  }
  deltaJ = fix_B_gauge_J();
  cout << "dev(J/final) = " << deltaJ << endl;

  return;
}
