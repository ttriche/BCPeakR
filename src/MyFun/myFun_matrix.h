
#ifndef  _myFun_matrix_h
#define  _myFun_matrix_h

#include	<iostream>
#include	<fstream>
#include        <cassert>
#include	"../TNT/tnt.h"
#include	"../JAMA_C/my_jama_cholesky.h"
#include	"../JAMA_C/my_jama_lu.h"

using namespace std;
using namespace TNT;
using namespace JAMA;
/*  E = 2.7183, Pi = 3.1416 */

/**************************************************************************
 ***      Vectr/Matrix Calculations
 **************************************************************************/
/* In TNT,
    Vector: +, -, *(mul_element), dot_prod, 
    Matrix: +, -, mul_element, transpose, *(mat*mat, mat*vec)

   Here,
    Vector: mul_scalar, tensor_prod, 
    Matrix: mul_scalar, mul_vecmatvec, 


PutRowVectorInMatrix, GetRowVectorInMatrix, PutMatrixInArray3D,
GetMatrixInArray3D
*/



/* This function returns 'mat * scalar'. */
Matrix<double>  mul_scalar(const Matrix<double> mat, double scalar)
{
  Matrix<double> res(mat.num_rows(), mat.num_cols());
  for (int i=0; i<mat.num_rows(); i++)
    for (int j=0; j<mat.num_cols(); j++)
      res[i][j] = mat[i][j] * scalar;
  return res;
}

/* This function returns 'vec * scalar'. */
Vector<double>  mul_scalar(const Vector<double> vec, double scalar)
{
  Vector<double> res(vec.size());
  for (int i=0; i<vec.size(); i++)  res[i] = vec[i] * scalar;
  return res;
}


/* This function returns the tensor product of vecA and vecB. */
Matrix<double>  tensor_prod( const Vector<double> &vecA,
			     const Vector<double> &vecB)
{
  Matrix<double>  res(vecA.size(), vecB.size());
  for (int i=0; i<vecA.size(); i++)
    for (int j=0; j < vecB.size(); j++)
      res[i][j] = vecA[i] * vecB[j];
  return res;
}


double		mul_vecmatvec(const Vector<double> &vec, 
			const Matrix<double> &mat)
{
  assert( (vec.size()==mat.num_cols())&(mat.num_rows()==vec.size()) );
  double  res=0.0;
  for (int i=0; i<vec.size(); i++)    for (int j=0; j<vec.size(); j++)
    res += vec[i]*mat[i][j]*vec[j];
  return res;
}


/*  Matrix multiplication:  mat1 * mat2 * mat3 */
Matrix<double>   MatMulti(Matrix<double>  mat1, Matrix<double>  mat2,
			  Matrix<double>  mat3)
{
  assert( (mat1.num_cols()==mat2.num_rows()) &&
	  (mat2.num_cols()==mat3.num_rows()) );
  Matrix<double>   res(mat1.num_rows(),  mat3.num_cols());
  for (int i=0; i<mat1.num_rows(); i++)
    for (int j=0; j<mat3.num_cols(); j++) {
      res[i][j] = 0.0;
      for (int k=0; k<mat2.num_rows(); k++)
	for (int s=0; s<mat2.num_cols(); s++)
	  res[i][j] += mat1[i][k] * mat2[k][s] * mat3[s][j];
    }
  return  res;
}


/* This function returns the inverse of the matrix. If it is not square matrix, abort.
 * (Write on 05_Jul_04) (corrected it on 06_May_28)
 */
Matrix<double> inv(const Matrix<double> mat)
{
  assert( mat.num_rows() == mat.num_cols() );
  Cholesky   chol(mat);
  Matrix<double>  id(mat.num_rows(), mat.num_rows());
  for (int i=0; i<mat.num_rows(); i++) for(int j=0; j<mat.num_rows(); j++)
    id[i][j] = 0.0;
  for (int i=0; i<mat.num_rows(); i++)  id[i][i] = 1.0;
  LU lu(mat);
  return lu.solve(id);
/*
  return  chol.solve(id);
	this doesn't work for certain cases, the problem is solved by
  http://wiki.cs.princeton.edu/index.php/TNT
*/
}



/* This function returns the diterminant of the matrix */
double	det(const Matrix<double> mat)
{
  LU	  lu(mat);
  return  lu.det();
}

/*
double	trace(const Matrix<double> &mat)
{
  double value = 0.0;
  for (int i=0; i<mat.num_rows(); i++)  trace += mat[i][i];
  return value;
}
*/



/* These four functions put(get) rowVector (matrix) in matrix (array3d). 
   Notice the position variables are in the C++ convention. */
template <class T>
static void   PutRowVectorInMatrix(Matrix<T> &mat, int rowNum, 
				   Vector<T> vec)
{
  assert(mat.num_cols() == vec.size());
  for (int j=0; j<mat.num_cols(); j++)
    mat[rowNum][j] = vec[j];
}

template <class T>  
Vector<T>   GetRowVectorInMatrix(Matrix<T> &mat, int rowNum)
{
  Vector<T>  vec(mat.num_cols());
  for (int i=0; i<mat.num_cols(); i++)
    vec[i] = mat[rowNum][i];
  return vec;
}

template <class T>
Vector<T>   GetColVectorInMatrix(Matrix<T> &mat, int colNum)
{
  Vector<T>  vec(mat.num_rows());
  for (int i=0; i<mat.num_rows(); i++)
    vec[i] = mat[i][colNum];
  return vec;
}

template <class T>
static void   PutColVectorInMatrix(Matrix<T> &mat, int colNum,
                                   Vector<T> vec)
{
  assert(mat.num_rows() == vec.size());
  for (int j=0; j<mat.num_rows(); j++)
    mat[j][colNum] = vec[j];
}


template <class T>
Matrix<T>     ColBind(Vector<T> vec1, Vector<T> vec2)
{
  assert( vec1.size() == vec2.size() );
  int	      N = vec1.size();
  Matrix<T>   newMat(N, 2);
  for (int i=0; i<N; i++) {
    newMat[i][0] = vec1[i];     newMat[i][1] = vec2[i];
  }
  return        newMat;
}


template <class T>
Matrix<T>     ColBind(Matrix<T> mat,  Vector<T> vec)
{
  assert( mat.num_rows() == vec.size() );
  int	      nRows=mat.num_rows(), nCols=mat.num_cols()+1;
  Matrix<T>   newMat(nRows, nCols);
  for (int i=0; i<nRows; i++) for(int j=0; j<nCols-1; j++)
    newMat[i][j] = mat[i][j];
  for (int i=0; i<nRows; i++)    newMat[i][nCols-1] = vec[i];
  return	newMat;
}

template <class T>
Matrix<T>     ColBind(Matrix<T> mat1, Matrix<T> mat2)
{
  assert( mat1.num_rows() == mat2.num_rows() );
  int	      nRows = mat1.num_rows(), nCols1 = mat1.num_cols(), 
	      nCols2 = mat2.num_cols();
  Matrix<T>   newMat(nRows, nCols1+nCols2);

  for (int i=0; i<nRows; i++) {
    for (int j=0; j<nCols1; j++)	newMat[i][j] = mat1[i][j];
    for (int j=0; j<nCols2; j++)	newMat[i][nCols1+j] = mat2[i][j];
  }
  return      newMat;
}





/* This function assign matrix in Array3D at index [position][*][*] */
template <class T>  
static void   PutMatrixInArray3D(Array3D<T> &arr3d, int position,
				 Matrix<T> mat)
{
  assert( (arr3d.dim2()==mat.num_rows())&(arr3d.dim3()==mat.num_cols()) );
  for (int i=0; i<mat.num_rows(); i++)
    for (int j=0; j<mat.num_cols(); j++)
      arr3d[position][i][j] = mat[i][j];
}

template <class T>
Matrix<T>  GetMatrixInArray3D(Array3D<T> &arr3d, int position)
{
  Matrix<double> mat(arr3d.dim2(), arr3d.dim3());
  for (int i=0; i<arr3d.dim2(); i++)
    for (int j=0; j<arr3d.dim3(); j++)
      mat[i][j] = arr3d[position][i][j];
  return mat;
}

template <class T>
static void   PutRowVectorInArray3D(Array3D<T> &arr3d, int posInDim1,
		int posInDim2, Vector<T> vec)
{
  assert( arr3d.dim3() == vec.size() );
  for (int i=0; i<vec.size(); i++)
    arr3d[posInDim1][posInDim2][i] = vec[i];
}

template <class T>
Vector<T>  GetRowVectorInArray3D(Array3D<T> &arr3d, int posInDim1,
                int posInDim2)
{
  Vector<T>  vec(arr3d.dim3());
  for (int i=0; i<vec.size(); i++)
    vec[i] = arr3d[posInDim1][posInDim2][i];
  return  vec;
}


/* Put vector in  arr[ posInDim1, ..., posInDim3] */
template <class T>
static void   PutColVectorInArray3D(Array3D<T> &arr3d, int posInDim1,
                int posInDim3, Vector<T> vec)
{
  assert( arr3d.dim2() == vec.size() );
  for (int i=0; i<vec.size(); i++)
    arr3d[posInDim1][i][posInDim3] = vec[i];
}

/* Get vector in  arr[ posInDim1, ..., posInDim3] */
template <class T>
Vector<T>  GetColVectorInArray3D(Array3D<T> &arr3d, int posInDim1,
                int posInDim3)
{
  Vector<T>  vec(arr3d.dim2());
  for (int i=0; i<vec.size(); i++)
    vec[i] = arr3d[posInDim1][i][posInDim3];
  return  vec;
}

#endif

