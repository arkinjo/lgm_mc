#ifndef BASE_H_
#define BASE_H_

#include <cmath>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

typedef std::vector<std::string> string_vector;
typedef std::vector<bool> bool_vector;

typedef boost::multi_array<double,1> darray1;
typedef boost::multi_array<double,2> darray2;
typedef boost::multi_array<double,3> darray3;
typedef boost::multi_array<double,4> darray4;

typedef boost::multi_array<bool,2> barray2;
typedef boost::multi_array<bool,3> barray3;
typedef boost::multi_array<bool,4> barray4;

const int nstate = 20;
const int nstate_match = 20;
const int nstate_delete = 1;
const int nstate_core = nstate_match;
const int nstate_insert = 20;
const int terminal_state = nstate_match;
const char terminal_residue = '-';

const double match_fraction = 0.75;

const double neg_infinity = -std::numeric_limits<double>::infinity();
const double pos_infinity = std::numeric_limits<double>::infinity();

//                            01234567890123456789
const std::string amino1_core =  "ACDEFGHIKLMNPQRSTVWY-";
const std::string amino1_insert = "acdefghiklmnpqrstvwy";
int find_amino1_core(const char a);
int find_amino1_insert(const char a);

bool amino1_core_p(const char a);
bool amino1_insert_p(const char a);

double log_add(const double a, const double b);
double log_dot(darray1& a, darray1& b, const int n);
void log_mat_mult(darray2& a, darray2& b, darray2& c,
		  const int l, const int m, const int n);

#endif //BASE_H_
