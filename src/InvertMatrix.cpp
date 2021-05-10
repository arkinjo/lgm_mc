
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

namespace ublas = boost::numeric::ublas;
/* Matrix inversion routine.
   Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */

bool InvertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse)
{
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  matrix<double> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());
  // perform LU-factorization
  int res = lu_factorize(A,pm);
  if( res != 0 ) return false;
  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<double>(A.size1()));
  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);
  return true;
}


bool InvertMatrixGen(const ublas::matrix<double>& input, ublas::matrix<double>& inverse)
{
  const int n = input.size1();
  gsl_matrix *A = gsl_matrix_alloc(n,n);
  gsl_matrix *V = gsl_matrix_alloc(n,n);
  gsl_vector *S = gsl_vector_alloc(n);
  gsl_vector *work = gsl_vector_alloc(n);

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      gsl_matrix_set(A,i,j,input(i,j));
      inverse(i,j) = 0.0;
    }
  }

  gsl_linalg_SV_decomp(A,V,S,work);

  const double smax = gsl_vector_get(S,0);
  for(int k = 0; k < n; ++k) {
    const double s = gsl_vector_get(S,k);
    if(s/smax > 1.0e-8) {
      for(int i = 0; i < n; ++i) {
	for(int j = 0; j < n; ++j) {
	  inverse(i,j) += gsl_matrix_get(V,i,k) * gsl_matrix_get(A,j,k) / s;
	}
      }
    }
  }

  gsl_matrix_free(A);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);
  return true;
}
