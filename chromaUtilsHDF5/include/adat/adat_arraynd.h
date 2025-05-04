// -*- C++ -*-
/*! @file
 * @brief Multi-dimensional arrays
 * 
 * Support for reference semantic multi-dimensional arrays
 */

#ifndef __adat_arraynd_h__
#define __adat_arraynd_h__

#include <iostream>
#include "xml_array.h"
#include "xml_array2d.h"
#include "xml_array3d.h"

namespace ADAT
{
  using namespace XMLArray;

  //----------------------------------------------------------------------------------
  //! Container for a generic N dimensional array
  template<class T> class ArrayNd
  {
  public:
    ArrayNd() {F=0;}
    explicit ArrayNd(const Array<int>& _nz) {F=0;resize(_nz);}
    ~ArrayNd() {delete[] F;}

    //! Copy constructor
    ArrayNd(const ArrayNd& s): nz(s.nz), sz(s.sz), F(0)
      {
	resize(nz);

	for(int i=0; i < sz; ++i)
	  F[i] = s.F[i];
      }

    //! Allocate mem for the array
    void resize(const Array<int>& _nz) 
      {
	delete[] F; 
	nz = _nz;
	sz = nz[0];
	for(int i=1; i < nz.size(); ++i)
	  sz *= nz[i];
//	F = new(nothrow) T[sz];
	F = new T[sz];
	if ( F==0x0 ) { 
	  std::cerr << "Unable to new memory in ArrayNd::resize():  sz= " << sz << "  size= ";
	  for(int i=0; i < _nz.size(); ++i) {
	    std::cerr << " " << _nz[i];
	  }
	  std::cerr << std::endl;
	  exit(1);
	}
      }

    //! Allocate mem for the array
    void resize(int n1) 
      {
	Array<int> _nz(1);
	_nz[0] = n1;

	resize(_nz);
      }

    //! Allocate mem for the array
    void resize(int n2, int n1) 
      {
	Array<int> _nz(2);
	_nz[0] = n1;
	_nz[1] = n2;

	resize(_nz);
      }

    //! Allocate mem for the array
    void resize(int n3, int n2, int n1) 
      {
	Array<int> _nz(3);
	_nz[0] = n1;
	_nz[1] = n2;
	_nz[2] = n3;

	resize(_nz);
      }

    //! Allocate mem for the array
    void resize(int n4, int n3, int n2, int n1) 
      {
	Array<int> _nz(4);
	_nz[0] = n1;
	_nz[1] = n2;
	_nz[2] = n3;
	_nz[3] = n4;

	resize(_nz);
      }

    //! Rank of the array - number of dimensions
    int rank() const {return nz.size();}

    //! Size of i-th array index. Indices run from left to right in operator() 
    /*! Note, the last/right index is the fastest varying index */
    int size(int i) const {return nz[i];}

    //! Size of an array containing sizes of each index.
    /*! Note, the last/right index is the fastest varying index */
    const Array<int>& size() const {return nz;}

    //! Number of elements in the array
    /*! The number of elements is the product of the sizes */
    int numElem() const {return sz;}

    //! Equal operator uses underlying = of T
    ArrayNd<T>& operator=(const ArrayNd<T>& s1)
      {
	resize(s1.size());

	for(int i=0; i < sz; ++i)
	  F[i] = s1.F[i];
	return *this;
      }

    //! Equal operator on an Array
    ArrayNd<T>& operator=(const Array<T>& s1)
      {
	resize(s1.size());

	for(int i=0; i < sz; ++i)
	  F[i] = s1[i];
	return *this;
      }

    //! Equal operator on an Array2d
    ArrayNd<T>& operator=(const Array2d<T>& s1)
      {
	resize(s1.nrows(), s1.ncols());

	for(int i=0; i < s1.nrows(); ++i)
	  for(int j=0; j < s1.ncols(); ++j)
	    (*this)(i,j) = s1(i,j);

	return *this;
      }

    //! Equal operator on an Array3d
    ArrayNd<T>& operator=(const Array3d<T>& s1)
      {
	resize(s1.leftSize(), s1.middleSize(), s1.rightSize());

	for(int i=0; i < s1.leftSize(); ++i)
	  for(int j=0; j < s1.middleSize(); ++j)
	    for(int k=0; k < s1.rightSize(); ++k)
	      (*this)(i,j,k) = s1(i,j,k);

	return *this;
      }

    //! Equal operator uses underlying = of T
    template<class T1>
    ArrayNd<T>& operator=(const T1& s1)
      {
	if (F == 0)
	{
	  std::cerr << "ArrayNd: left hand side not initialized in =" << std::endl;
	  exit(1);
	}

	for(int i=0; i < sz; ++i)
	  F[i] = s1;
	return *this;
      }

    //! Return ref to an element
    T& operator()(int i)
      {
	if (nz.size() != 1)
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	return F[i];
      }

    //! Return const ref to an element
    const T& operator()(int i) const
      {
	if (nz.size() != 1)
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	return F[i];
      }

    //! Return ref to an element
    T& operator()(int j, int i)
      {
	if (nz.size() != 2)
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	return F[i+nz[0]*j];
      }

    //! Return const ref to an element
    const T& operator()(int j, int i) const
      {
	if (nz.size() != 2)
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	return F[i+nz[0]*j];
      }

    //! Return ref to an element
    T& operator()(int k, int j, int i) 
      {
	if (nz.size() != 3)
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	return F[i+nz[0]*(j+nz[1]*(k))];
      }

    //! Return const ref to an element
    const T& operator()(int k, int j, int i) const
      {
	if (nz.size() != 3)
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	return F[i+nz[0]*(j+nz[1]*(k))];
      }

    //! Return ref to an element
    T& operator()(int l, int k, int j, int i) 
      {
	if (nz.size() != 4)
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	return F[i+nz[0]*(j+nz[1]*(k+nz[2]*l))];
      }

    //! Return const ref to an element
    const T& operator()(int l, int k, int j, int i) const
      {
	if (nz.size() != 4)
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	return F[i+nz[0]*(j+nz[1]*(k+nz[2]*l))];
      }

    //! Return ref to an element via indices packed in a Array array
    T& operator[](const Array<int>& ind)
      {
	if (ind.size() != nz.size())
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	int off = ind[0];
	for(int i=1; i < nz.size(); ++i)
	  off = off*nz[i] + ind[i];

	return F[off];
      }

    //! Return ref to an element via indices packed in a Array array
    const T& operator[](const Array<int>& ind) const
      {
	if (ind.size() != nz.size())
	{
	  std::cerr << "ArrayNd: improper rank of array indices" << std::endl;
	  exit(1);
	}

	int off = ind[0];
	for(int i=1; i < nz.size(); ++i)
	  off = off*nz[i] + ind[i];

	return F[off];
      }

    //! Return ref to an element with index flattened over indices
    /*! Right index is fastest varying */
    T& getElem(int off)
      {
	if (off < 0 || off >= sz)
	{
	  std::cerr << "ArrayNd: index out of bounds" << std::endl;
	  exit(1);
	}

	return F[off];
      }

    //! Return const-ref to an element with index flattened over indices
    /*! Right index is fastest varying */
    const T& getElem(int off) const
      {
	if (off < 0 || off >= sz)
	{
	  std::cerr << "ArrayNd: index out of bounds" << std::endl;
	  exit(1);
	}

	return F[off];
      }

  private:
    Array<int> nz;
    int sz;
    T *F;
  };


}

#endif
