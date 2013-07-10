
#ifndef _myFun_IO_h
#define _myFun_IO_h

#include        <iostream>
#include        <fstream>
#include        <cassert>
#include        "../TNT/tnt.h"
#include        "../JAMA_C/my_jama_cholesky.h"
#include        "../JAMA_C/my_jama_lu.h"

using namespace std;
using namespace TNT;
using namespace JAMA;


template <class T>
static void      ReadValue(string  filename, T &numSub)
{
  ifstream        infile;  
  infile.open( filename.c_str() );
  if (!infile)    {
    cout << "Couldn't open file!" << endl;
    return;
  }
  while (infile.good()) infile >> numSub;
  infile.close();
}

/* This function reads the vector into a file.
 */
template <class T>
static void   ReadVector(string filename, Vector<T> &vec)
{
  ifstream        infile;
  int             i=0;
  infile.open( filename.c_str() );
  if (!infile)    {
    cout << "Couldn't open file!" << endl;
    return;
  }
  while (infile.good())   infile >> vec[i++];
  infile.close();
}

Matrix<double>   ReadMatrix(string filename)
{
  ifstream        infile;
  infile.open( filename.c_str() );
  if (!infile) {
    cout << "Couldn't open file!" << endl;
    return  Matrix<double>();
  }

  int num_rows, num_cols;
  infile >> num_rows >> num_cols ;
  Matrix<double>  mat(num_rows, num_cols);
  if (infile.good())
    for (int i=0; i<num_rows; i++) for (int j=0; j<num_cols; j++)
      infile >> mat[i][j];

  infile.close();
  return mat;
}


Array3D<double>  ReadArray3D(string filename)
{
  ifstream  infile;
  infile.open( filename.c_str() );
  if (!infile) {
    cout << "Couldn't open file!" << endl;
    return  Array3D<double>();
  }

  int  dim1, dim2, dim3;
  infile >> dim1 >> dim2 >> dim3;
  Array3D<double> arr3d(dim1, dim2, dim3);
  if (infile.good())
    for (int i=0; i<dim1; i++) for(int j=0; j<dim2; j++) 
      for(int k=0; k<dim3; k++)  infile >> arr3d[i][j][k];
  infile.close();
  return arr3d;
}

/* This function writes the vector into a file.
 */

template <class T>
static void   WriteVector(string filename, const Vector<T> &vec)
{
  ofstream        outfile;

  outfile.open( filename.c_str() );
  if (!outfile)   {
    cout << "Couldn't open file!" << endl;
    return;
  }
  for (int i=0; i<vec.size(); i++)    outfile << vec[i] << "\n";
  outfile.close();
}

/* This function writes a matrix into a file.
 */

template <class T>
static void   WriteMatrix(string filename, const Matrix<T> mat)
{
  ofstream        outfile;
  
  outfile.open( filename.c_str() );
  if (!outfile)   {
    cout << "Couldn't open file!" << endl;
    return;
  }
  outfile << mat.num_rows() << "\t" << mat.num_cols() << "\n" << endl;
  for (int i=0; i<mat.num_rows(); i++) {
    for (int j=0; j < mat.num_cols(); j++) outfile << mat[i][j] << " ";
    outfile << "\n";
  }
  outfile.close();
}

template <class T>
static void   WriteArray3D(string filename, const Array3D<T> &arr3d)
{
  ofstream        outfile;
  
  outfile.open( filename.c_str() );
  if (!outfile)   {
    cout << "Couldn't open file!" << endl;
    return;
  }
  outfile << arr3d.dim1() << "\t" << arr3d.dim2() << "\t" <<
    arr3d.dim3() << "\n" << endl;
  for (int i=0; i<arr3d.dim1(); i++) {
    for (int j=0; j<arr3d.dim2(); j++) {
      for (int k=0; k<arr3d.dim3(); k++) outfile<<arr3d[i][j][k]<<" ";
      outfile << "\n";
    }
    outfile << "\n" << endl;
  }
}

/* Initialze vectors, matrices, and arrays */
template <class T>
static void     Initialize(Vector<T>  &vec)
{
  T   zero=vec[0]-vec[0];
  for (int i=0; i<vec.dim(); i++)   vec[i] = zero;
}

template <class T>
static void     Initialize(Matrix<T>  &mat)
{
  T   zero=mat[0][0]-mat[0][0];
  for (int i=0; i<mat.num_rows(); i++) for(int j=0;j<mat.num_cols();j++)
    mat[i][j] = zero;
}


template <class T>
static void     Initialize(Array3D<T>  &arr)
{
  T  zero=arr[0][0][0]-arr[0][0][0];
  for (int i=0; i<arr.dim1(); i++)  for(int j=0; j<arr.dim2(); j++)
    for (int k=0; k<arr.dim3(); k++)     arr[i][j][k] = zero;
}


#endif


