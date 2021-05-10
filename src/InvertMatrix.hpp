#ifndef INVERT_MATRIX_HPP
#define INVERT_MATRIX_HPP

 // REMEMBER to update "lu.hpp" header includes from boost-CVS
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;
/* Matrix inversion routine.
   Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */

bool InvertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse);
bool InvertMatrixGen(const ublas::matrix<double>& input, ublas::matrix<double>& inverse);
#endif //INVERT_MATRIX_HPP
