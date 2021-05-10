#ifndef MSA_H_
#define MSA_H_

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

class Msa
{
private:
  vector<string> headers;
  vector<string> sequences;
  vector<string> qseqs;
  vector<string> sseqs;
  vector<int> qstarts;
  vector<int> qends;
  int base_length;
public:
  vector<double> weights;
  Msa();
  Msa(const string file);
  Msa(const string file, const bool blast_outp);
  ~Msa();
  void read_from_file(const string);
  void read_from_file(const char *);
  void read_from_blast_out(const string);
  void read_from_blast_out(const char *);
  
  void set_weights(void);
  inline  void set_base_length(const int len) { base_length = len; };
  inline  int get_base_length(void) const { return base_length; };
  inline int num_seq(void) const {return headers.size();};
  inline string header(const int i) const { return headers[i]; };
  inline string sequence(const int i) const { return sequences[i]; };
  inline size_t length(void) const { return sequences[0].length(); };
  inline double get_weight(const int i) const  { return weights[i]; };
  inline int qstart(const int i) const { return qstarts[i]; };
  inline int qend(const int i) const { return qends[i]; };
  inline string qseq(const int i) const { return qseqs[i];};
  inline string sseq(const int i) const { return sequences[i];};
};

  
#endif //MSA_H_
