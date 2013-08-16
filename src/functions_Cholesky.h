#ifndef JAMA_CHOLESKY_H
#define JAMA_CHOLESKY_H

#include "math.h"
/* needed for sqrt() below. */

namespace JAMA
{

	using namespace std; // use of "cout<< << endl" only
	using namespace TNT;

/**
   <P>
   For a symmetric, positive definite matrix A, this function
   computes the Cholesky factorization, i.e. it computes a lower 
   triangular matrix L such that A = L*L'.
   If the matrix is not symmetric or positive definite, the function
   computes only a partial decomposition.  This can be tested with
   the is_spd() flag.

   Typical usage looks like:
 
 	 	Matrix<double> A(n,n);
	 	Matrix<double> L;
	 	... 
	 	Cholesky chol(A);

	 	if (chol.is_spd())
			L = chol.getL();		
  	else
			cout << "factorization was not complete.\n";

	(Adapted from JAMA, a Java Matrix Library, developed by jointly 
	by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).

  05_July_04, downloaded this folder from TNT website http://math.nist.gov/tnt/
       			  changed 'Array1D, Arra2D and Real' to 'Vector, Matrix and double'
  */
	class Cholesky 
	{ // {{{
		Matrix<double> L_;	// lower triangular factor
		int isspd;		// 1 if matrix to be factored was SPD

		public:
		Cholesky();
		Cholesky(const Matrix<double> &A);
		Matrix<double> getL() const;
		Vector<double> solve(const Vector<double> &B);
		Matrix<double> solve(const Matrix<double> &B);
		Matrix<double> inv();
		int is_spd() const;
	}; // }}}


	Cholesky::Cholesky() : L_(0,0), isspd(0) {}


	// @return 1 if original matrix to be factored is symmetric positive-definite
	int Cholesky::is_spd() const
	{ // {{{
		return isspd;
	} // }}}


	// @return the lower triangular factor, L, such that L*L'=A.
	Matrix<double> Cholesky::getL() const
	{ // {{{
		return L_;
	} // }}}

  
	/** Constructs a lower triangular matrix L, such that L*L'= A.
			If A is not symmetric positive-definite (SPD), only a
			partial factorization is performed.  If is_spd() is true (1), 
			then the factorization was successful.
	*/
	Cholesky::Cholesky(const Matrix<double> &A)
	{ // {{{
		int m = A.num_rows(), n = A.num_cols();
		isspd = (m == n);
		if (m != n) {
			L_ = Matrix<double>(0,0);
			return;
		}
		L_ = Matrix<double>(n,n);

		// Main loop.
		for (int j = 0; j < n; j++) {
			double d = 0.0;
			for (int k = 0; k < j; k++) {
				double s = 0.0;
				for (int i = 0; i < k; i++) { 
					s += L_[k][i]*L_[j][i]; 
				}
				L_[j][k] = s = (A[j][k] - s)/L_[k][k];
				d = d + s*s;
				isspd = isspd && (A[k][j] == A[j][k]); 
			}
			d = A[j][j] - d;
			isspd = isspd && (d > 0.0);
			L_[j][j] = sqrt(d > 0.0 ? d : 0.0);
			for (int k = j+1; k < n; k++)  {
				L_[j][k] = 0.0; 
			}
		}
	} // }}}


	/** Solve a linear system A*x = b, using the previously computed
		cholesky factorization of A: L*L'.

		 @param  B   A Matrix with as many rows as A and any number of columns.
		 @return     x so that L*L'*x = b.  If b is nonconformat, or if A
									 was not symmetric posidtive definite, a null (0x0)
									 array is returned.
	*/
	Vector<double> Cholesky::solve(const Vector<double> &b)
	{ // {{{
		int n = L_.num_rows();
		if (b.size() != n)	return Vector<double>();

		Vector<double>  x=b;

		// Solve L*y = b;
		for (int k = 0; k < n; k++) {
			for (int i = 0; i < k; i++) {
				x[k] -= x[i]*L_[k][i];
			}
			x[k] /= L_[k][k];
		}

		// Solve L'*X = Y;
		for (int k = n-1; k >= 0; k--) {
			for (int i = k+1; i < n; i++) {
				x[k] -= x[i]*L_[i][k];
			}
			x[k] /= L_[k][k];
		}

		return x;
	} // }}}


	/** Solve a linear system A*X = B, using the previously computed
		cholesky factorization of A: L*L'.

		 @param  B   A Matrix with as many rows as A and any number of columns.
		 @return     X so that L*L'*X = B.  If B is nonconformat, or if A
				was not symmetric posidtive definite, a null (0x0)
				array is returned.
	*/
	Matrix<double> Cholesky::solve(const Matrix<double> &B)
	{ // {{{
		int n = L_.num_rows();
		if (B.num_rows() != n)   return Matrix<double>();

		Matrix<double> X=B;
		int nx = B.num_cols();

		// Solve L*y = b;
		for (int j=0; j< nx; j++) {
			for (int k = 0; k < n; k++) {
				for (int i = 0; i < k; i++) {
					X[k][j] -= X[i][j]*L_[k][i];
				}
				X[k][j] /= L_[k][k];
			}
		}

		// Solve L'*X = Y;
		for (int j=0; j<nx; j++) {
			for (int k = n-1; k >= 0; k--) {
				for (int i = k+1; i < n; i++) X[k][j] -= X[i][j]*L_[i][k];
				X[k][j] /= L_[k][k];
			}
		}

		return X;

	} // }}}

} // end of namespace JAMA

#endif 
// JAMA_CHOLESKY_H
