/*
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/



#ifndef TNT_ARRAY4D_H
#define TNT_ARRAY4D_H

#include <cstdlib>
#include <iostream>
#ifdef TNT_BOUNDS_CHECK
#include <assert.h>
#endif

#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "tnt_array3d.h"

namespace TNT
{

template <class T>
class Array4D 
{
  private:
  	Array1D<T> data_;
//	Array2D<T*> v_;
	Array3D<T*> u_;
	int m_, n_, g_, h_;

  public:
    typedef         T   value_type;
	Array4D();
	Array4D(int m, int n, int g, int h);

	inline Array4D(const Array4D &A);
	inline Array4D & operator=(const Array4D &A);
	inline Array4D & ref(const Array4D &A);
		   Array4D & inject(const Array4D & A);
	inline T*** operator[](int i);
	~Array4D();
};

template <class T>
Array4D<T>::Array4D() : data_(),u_(),m_(0),n_(0),g_(0),h_(0){}

template <class T>
Array4D<T>::Array4D(const Array4D<T> &A) : data_(A.data_), 
  u_(A.u_), m_(A.m_), n_(A.n_), g_(A.g_), h_(A.h_) {}

template <class T>
Array4D<T>::Array4D(int m, int n, int g, int h) : data_(m*n*g*h),
  u_(m,n,g), m_(m), n_(n), g_(g), h_(h)
{
  if (m>0 && n>0 && g>0 && h>0)  {	T* p = & (data_[0]);
    for (int i=0; i<m_; i++)	{	T* pingh = p+ i*n_*g_*h_;
      for (int j=0; j<n; j++) {		T* pjgh = pingh + j*g_*h_;
	for (int k=0; k<g; k++)		u_[i][j][k]=pjgh+k*h_;
      }
    }
  }
}

template <class T>
inline T*** Array4D<T>::operator[](int i)
{
#ifdef TNT_BOUNDS_CHECK
        assert(i >= 0);
        assert(i < m_);
#endif

return u_[i];

}

template <class T>
Array4D<T> & Array4D<T>::inject(const Array4D &A)
{
  if (A.m_ == m_ &&  A.n_ == n_ && A.g_ == g_ && A.h_ == h_)
  for (int i=0; i<m_; i++)  for (int j=0; j<n_; j++)
    for (int k=0; k<g_; k++)  for (int l=0; l<h_; l++)
      u_[i][j][k][l] = A.u_[i][j][k][l];
  return *this;
}

template <class T>
Array4D<T> & Array4D<T>::ref(const Array4D<T> &A)
{
  if (this != &A) {
    m_ = A.m_;	n_ = A.n_;    g_ = A.g_;  h_ = A.h_;
    u_ = A.u_;	data_ = A.data_;
  }
  return *this;
}

template <class T>
Array4D<T> & Array4D<T>::operator=(const Array4D<T> &A)
{
  return ref(A);
}

template <class T>
Array4D<T>::~Array4D() {}

} /* namespace TNT */

#endif
/* TNT_ARRAY4D_H */

