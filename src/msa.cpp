#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "msa.hpp"

using namespace std;

Msa::Msa() {}
Msa::Msa(const string file)
{
  read_from_file(file);
}
Msa::Msa(const string file, const bool blast_outp) 
{
  if(blast_outp) {
    read_from_blast_out(file);
  } else {
    read_from_file(file);
  }
}

Msa::~Msa() {}


void Msa::read_from_file(const string file)
{
  read_from_file(file.c_str());
}

void Msa::read_from_file(const char *file)
{
  ifstream fin;
  fin.open(file);
  if (!fin.is_open()) {
    cerr << "Msa: could not open: " << file << endl;
    exit(1);
  }

  string line;
  while(getline(fin,line)) {
    if(line.length() == 0) continue;
    if(line[0] == '>') {
      headers.push_back(line);
      sequences.push_back("");
    }
    else {
      sequences.back() += line;
    }
  }

  fin.close();

  set_weights();

  return;
}

void Msa::read_from_blast_out(const string file)
{
  read_from_blast_out(file.c_str());
}

void Msa::read_from_blast_out(const char *file)
{
  ifstream fin;
  fin.open(file);
  if (!fin.is_open()) {
    cerr << "Msa: could not open: " << file << endl;
    exit(1);
  }

  string line,hd,qseq,sseq;
  int qs,qe;
  while(getline(fin,line)) {
    if(line.length() == 0) continue;
    istringstream iss(line);
    iss >> hd >> qs >> qe >> qseq >> sseq;
    headers.push_back(hd);
    qstarts.push_back(qs);
    qends.push_back(qe);
    qseqs.push_back(qseq);
    sequences.push_back(sseq);
  }

  fin.close();

  set_weights();

  return;
}


void Msa::set_weights(void)
{
  const int nseq = headers.size();

  weights.resize(nseq,1.0);
  return;
}
