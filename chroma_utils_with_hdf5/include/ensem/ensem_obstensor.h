// -*- C++ -*-
// $Id: ensem_obstensor.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! \file
 * \brief Observable Tensor
 */


namespace ENSEM {

//-------------------------------------------------------------------------------------
/*! \addtogroup obstensor Tensor observable
 * \ingroup fiber
 *
 * Observable type that transforms like a tensor
 *
 * @{
 */

//! Observable Tensor class
/*!
 * All tensor classes inherit this class
 */
template <class T> class OTensor
{
public:
  OTensor() {sz=0;F=0;}
  ~OTensor() {delete[] F;}

  //---------------------------------------------------------
  //! Conversion constructor
  template<class T1>
  OTensor(const OTensor<T1>& rhs) : sz(0), F(0)
    {
      checkResize("convert OTensor", rhs.size());
      for(int i=0; i < sz; ++i) 
	elem(i) = rhs.elem(i);
    }

  //---------------------------------------------------------
  //! OTensor = zero
  inline
  OTensor& operator=(const Zero& rhs)
    {
      checkSize("OTensor = zero");
      for(int i=0; i < sz; ++i)
	elem(i) = zero;
      return *this;
    }

  //! OTensor = OScalar
  template<class T1>
  inline
  OTensor& operator=(const OScalar<T1>& rhs) 
    {
      checkSize("OTensor = const");
      for(int i=0; i < sz; ++i)
	elem(i) = rhs.elem();
      return *this;
    }

  //! OTensor += OScalar
  template<class T1>
  inline
  OTensor& operator+=(const OScalar<T1>& rhs) 
    {
      checkSize("OTensor += OTensor");
      for(int i=0; i < sz; ++i)
	elem(i) += rhs.elem();

      return *this;
    }

  //! OTensor -= OScalar
  template<class T1>
  inline
  OTensor& operator-=(const OScalar<T1>& rhs) 
    {
      checkSize("OTensor -= OTensor");
      for(int i=0; i < sz; ++i)
	elem(i) -= rhs.elem();

      return *this;
    }

  //! OTensor *= OScalar
  template<class T1>
  inline
  OTensor& operator*=(const OScalar<T1>& rhs) 
    {
      checkSize("OTensor *= OTensor");
      for(int i=0; i < sz; ++i)
	elem(i) *= rhs.elem();

      return *this;
    }

  //! OTensor /= OScalar
  template<class T1>
  inline
  OTensor& operator/=(const OScalar<T1>& rhs) 
    {
      checkSize("OTensor /= OScalar");
      for(int i=0; i < sz; ++i)
	elem(i) /= rhs.elem();

      return *this;
    }

  //------------------------------------------------
  //! OTensor = OTensor
  /*! Set equal to another OTensor */
  inline
  OTensor& operator=(const OTensor& rhs) 
    {
      checkResize("OTensor = OTensor", rhs.size());
      for(int i=0; i < sz; ++i)
	elem(i) = rhs.elem(i);

      return *this;
    }

  //! OTensor = OTensor
  /*! Set equal to another OTensor */
  template<class T1>
  inline
  OTensor& operator=(const OTensor<T1>& rhs) 
    {
      checkResize("OTensor = OTensor<T1>", rhs.size());
      for(int i=0; i < sz; ++i)
	elem(i) = rhs.elem(i);

      return *this;
    }

  //! OTensor += OTensor
  template<class T1>
  inline
  OTensor& operator+=(const OTensor<T1>& rhs) 
    {
      checkSize("OTensor += OTensor", rhs.size());
      for(int i=0; i < sz; ++i)
	elem(i) += rhs.elem(i);

      return *this;
    }

  //! OTensor -= OTensor
  template<class T1>
  inline
  OTensor& operator-=(const OTensor<T1>& rhs) 
    {
      checkSize("OTensor -= OTensor", rhs.size());
      for(int i=0; i < sz; ++i)
	elem(i) -= rhs.elem(i);

      return *this;
    }

  //! OTensor *= OTensor
  template<class T1>
  inline
  OTensor& operator*=(const OTensor<T1>& rhs) 
    {
      checkSize("OTensor *= OTensor", rhs.size());
      for(int i=0; i < sz; ++i)
	elem(i) *= rhs.elem(i);

      return *this;
    }

  //! OTensor /= OTensor
  template<class T1>
  inline
  OTensor& operator/=(const OTensor<T1>& rhs) 
    {
      checkSize("OTensor /= OTensor", rhs.size());
      for(int i=0; i < sz; ++i)
	elem(i) /= rhs.elem(i);

      return *this;
    }


  //! Deep copy constructor
  OTensor(const OTensor& rhs) : nz(rhs.nz), sz(rhs.sz), F(0)
    {
      F = new(std::nothrow) T[sz];
      if ( F==0x0 ) { 
	std::cerr << "Unable to new memory in OTensor::copy()" << std::endl;
	exit(1);
      }
      
      for(int i=0; i < sz; ++i)
	elem(i) = rhs.elem(i);
    }

  //---------------------------------------------------------
public:
  //! The backdoor
  /*! 
   * Used by optimization routines (e.g., SSE) that need the memory address of data.
   * BTW: to make this a friend would be a real pain since functions are templatized.
   */
  inline const T* getF() const {return F;}
  inline T* getF() {return F;}

public:
  //---------------------------------------------------------
  inline void checkSize(const char *s) const
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkSize: %s\n", s);
#endif

      if (sz == 0)
      {
	std::cerr << s << ": Invalid OTensor size" << std::endl;
	exit(1);
      }
    }

  inline void checkSize(const char *s, const Array<int>& n1) const
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkSize: %s, OTensor[%d]\n", s, n1);
#endif

      if (sz == 0 || nz.size() == 0 || nz.size() != n1.size())
      {
	std::cerr << s << ": Invalid OTensor dest size" << std::endl;
	exit(1);
      }

      for(int i=0; i < n1.size(); ++i)
      {
	if (nz[i] != n1[i])
	{
	  std::cerr << s << ": Invalid OTensor dest and/or source size" << std::endl;
	  exit(1);
	}
      }
    }

  inline void checkResize(const char *s, const Array<int>& n1)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, OTensor[rank=%d]\n", s, n1.size());
#endif

      if (n1.size() == 0)
      {
	std::cerr << "checkResize: " << s << ": invalid OTensor source size" << std::endl;
	exit(1);
      }
      resize(n1);
    }

  inline void checkResize(const char *s, const Array<int>& n1, const Array<int>& n2)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, OTensor[%d]\n", s, n1);
#endif

      if (n1 != n2)
      {
	std::cerr << "checkResize: " << s << ": invalid OTensor source size" << std::endl;
	exit(1);
      }
      resize(n1);
    }

  //! Has this object been initialized (resized and such)
  bool initP() const 
    {
      return (sz != 0) ? true : false;
    }

  //! Size of i-th array index. Indices run from left to right in operator() 
  /*! Note, the last/right index is the fastest varying index */
  int size(int i) const {return nz[i];}

  //! Size of an array containing sizes of each index.
  /*! Note, the last/right index is the fastest varying index */
  const Array<int>& size() const {return nz;}

  //! Total length of array
  int numElem() const {return sz;}

  //! Resize the number of configs
  void resize(const Array<int>& _nz) 
    {
      delete[] F; 
      nz = _nz;
      sz = nz[0];
      for(int i=1; i < nz.size(); ++i)
	sz *= nz[i];
      F = new(std::nothrow) T[sz];
      if ( F==0x0 ) { 
	std::cerr << "Unable to new memory in OTensor::resize()" << std::endl;
	exit(1);
      }
    }

public:
  inline void resize(const OTensor& a) {resize(a.nz);}

public:
  //! Return ref to an element via indices packed in a Array array
  T& operator[](const Array<int>& ind)
    {
      if (ind.size() != nz.size())
      {
	std::cerr << "OTensor: improper rank of array indices\n";
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
	std::cerr << "OTensor: improper rank of array indices\n";
	exit(1);
      }

      int off = ind[0];
      for(int i=1; i < nz.size(); ++i)
	off = off*nz[i] + ind[i];

      return F[off];
    }

public:
  inline T& elem(int i) {return F[i];}
  inline const T& elem(int i) const {return F[i];}

private:
  Array<int> nz;
  int sz;
  T  *F;
};


// I/O
//! Binary input
template<class T>  
inline
void read(ADATIO::BinaryReader& bin, OTensor<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.numElem(); ++i)
    read(bin, d.elem(i));
}

//! Binary output
template<class T>
inline
void write(ADATIO::BinaryWriter& bin, const OTensor<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.numElem(); ++i)
    write(bin, d.elem(i));
}

//! Stream input
template<class T>  
inline
std::istream& operator>>(std::istream& s, OTensor<T>& d)
{
  d.checkSize(__func__);

  for(int i=0; i < d.numElem(); ++i)
    s >> d.elem(i);

  return s;
}

//! Stream output
template<class T>  
inline
std::ostream& operator<<(std::ostream& s, const OTensor<T>& d)
{
  d.checkSize(__func__);

  for(int k = 0; k < d.numElem(); ++k)
  {
    s << d.elem(k);
    if (k < d.numElem()-1)
      s << "\n";
  }

  return s;
}


//! Text input
template<class T>  
inline
ENSEM::TextReader& operator>>(ENSEM::TextReader& txt, OTensor<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.numElem(); ++i)
  {
    int k;
    txt >> k >> d.elem(i);
    if (k != i)
    {
      std::cerr << "error reading OTensor" << std::endl;
      exit(1);
    }
  }

  return txt;
}

//! Text output
template<class T>  
inline
ENSEM::TextWriter& operator<<(ENSEM::TextWriter& txt, const OTensor<T>& d)
{
  d.checkSize(__func__);

  for(int k = 0; k < d.numElem(); ++k)
  {
    txt << k << " " << d.elem(k);
    if (k < d.numElem()-1)
      txt << "\n";
  }

  return txt;
}


//! XML output
template<class T> 
inline
ADATXML::XMLWriter& operator<<(ADATXML::XMLWriter& xml, const OTensor<T>& d)
{
  xml.openTag("Tensor");

  XMLWriterAPI::AttributeList alist;

  // Copy into another array first
  for(int i=0; i < d.numElem(); ++i)
  {
    alist.clear();
    alist.push_back(XMLWriterAPI::Attribute("row", i));

    xml.openTag("elem", alist);
    xml << d.elem(i);
    xml.closeTag();
  }

  xml.closeTag();  // Tensor
  return xml;
}

/*! @} */  // end of group obstensor


//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1>
struct WordType<OTensor<T1> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

// Internally used scalars
template<class T>
struct InternalScalar<OTensor<T> > {
  typedef OScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive scalar leaving grid alone
template<class T>
struct PrimitiveScalar<OTensor<T> > {
  typedef OScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving observable indices alone
template<class T>
struct EnsemScalar<OTensor<T> > {
  typedef OTensor<typename EnsemScalar<T>::Type_t>  Type_t;
};

// Traits class to label IO types
template<class T> 
struct EnsbcIO<OTensor<T> > {
  enum {type = 2 + EnsbcIO<T>::type};
};

//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(OTensor) -> OTensor
template<class T1, class Op>
struct UnaryReturn<OTensor<T1>, Op> {
  typedef OTensor<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};
// Default binary(OScalar,OTensor) -> OTensor
template<class T1, class T2, class Op>
struct BinaryReturn<OScalar<T1>, OTensor<T2>, Op> {
  typedef OTensor<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(OTensor,OScalar) -> OTensor
template<class T1, class T2, class Op>
struct BinaryReturn<OTensor<T1>, OScalar<T2>, Op> {
  typedef OTensor<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(OTensor,OTensor) -> OTensor
template<class T1, class T2, class Op>
struct BinaryReturn<OTensor<T1>, OTensor<T2>, Op> {
  typedef OTensor<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};


#if 0
template<class T1, class T2>
struct UnaryReturn<OScalar<T2>, OpCast<T1> > {
  typedef OScalar<typename UnaryReturn<T, OpCast>::Type_t>  Type_t;
//  typedef T1 Type_t;
};
#endif


// Assignment is different
template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OTensor<T2>, OpAssign > {
  typedef OTensor<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OTensor<T2>, OpAddAssign > {
  typedef OTensor<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OTensor<T2>, OpSubtractAssign > {
  typedef OTensor<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OScalar<T2>, OpMultiplyAssign > {
  typedef OTensor<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OScalar<T2>, OpDivideAssign > {
  typedef OTensor<T1> &Type_t;
};
 



//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

/*! \addtogroup obstensor */
/*! @{ */

// Observable Tensors

template<class T1>
inline typename UnaryReturn<OTensor<T1>, OpUnaryPlus>::Type_t
operator+(const OTensor<T1>& l)
{
  typename UnaryReturn<OTensor<T1>, OpUnaryPlus>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = +l.elem(i);
  return d;
}


template<class T1>
inline typename UnaryReturn<OTensor<T1>, OpUnaryMinus>::Type_t
operator-(const OTensor<T1>& l)
{
  typename UnaryReturn<OTensor<T1>, OpUnaryMinus>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = -l.elem(i);
  return d;
}


// OTensor + OTensor
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OTensor<T2>, OpAdd>::Type_t
operator+(const OTensor<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OTensor<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) + r.elem(i);
  return d;
}

// OTensor + OScalar
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OScalar<T2>, OpAdd>::Type_t
operator+(const OTensor<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OScalar<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) + r.elem();
  return d;
}

// OScalar + OTensor
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OTensor<T2>, OpAdd>::Type_t
operator+(const OScalar<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OTensor<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem() + r.elem(i);
  return d;
}


// OTensor - OTensor
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OTensor<T2>, OpSubtract>::Type_t
operator-(const OTensor<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OTensor<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) - r.elem(i);
  return d;
}

// OTensor - OScalar
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OScalar<T2>, OpSubtract>::Type_t
operator-(const OTensor<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OScalar<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) - r.elem();
  return d;
}

// OScalar - OTensor
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OTensor<T2>, OpSubtract>::Type_t
operator-(const OScalar<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OTensor<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem() - r.elem(i);
  return d;
}


//! OTensor = OTensor * OTensor
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OTensor<T2>, OpMultiply>::Type_t
operator*(const OTensor<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OTensor<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) * r.elem(i);

  return d;
}

//! OTensor = OTensor * OScalar
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OScalar<T2>, OpMultiply>::Type_t
operator*(const OTensor<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OScalar<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) * r.elem();

  return d;
}

//! OTensor = OScalar * OTensor
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OTensor<T2>, OpMultiply>::Type_t
operator*(const OScalar<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OTensor<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem() * r.elem(i);

  return d;
}


//! OTensor / OTensor
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OTensor<T2>, OpDivide>::Type_t
operator/(const OTensor<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OTensor<T2>, OpDivide>::Type_t  d;
  Array<int> nz(concat(l.size(), r.size()));
  d.resize(nz);
  const int nr = r.numElem();

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) / r.elem(i);
  return d;
}

//! OTensor / OScalar
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OScalar<T2>, OpDivide>::Type_t
operator/(const OTensor<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OScalar<T2>, OpDivide>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) / r.elem();
  return d;
}

//! OScalar / OTensor
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OTensor<T2>, OpDivide>::Type_t
operator/(const OScalar<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OTensor<T2>, OpDivide>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem() / r.elem(i);
  return d;
}



//-----------------------------------------------------------------------------
// Functions


// OTensor = adj(OTensor)
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnAdjoint>::Type_t
adj(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnAdjoint>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = adj(s1.elem(i));
  return d;
}


// OTensor = conj(OTensor)
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnConjugate>::Type_t
conj(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnConjugate>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = conj(s1.elem(i));
  return d;
}


// OTensor = transpose(OTensor)
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnTranspose>::Type_t
transpose(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnTranspose>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = transpose(s1.elem(i));
  return d;
}


// OTensor = Trace(OTensor)
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnTrace>::Type_t
trace(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnTrace>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = trace(s1.elem(i));
  return d;
}


// OTensor = Re(Trace(OTensor))
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnRealTrace>::Type_t
trace_real(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnRealTrace>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = trace_real(s1.elem(i));
  return d;
}


// OTensor = Im(Trace(OTensor))
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnImagTrace>::Type_t
trace_imag(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnImagTrace>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = trace_imag(s1.elem(i));
  return d;
}


// OTensor = Re(OTensor)
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnReal>::Type_t
real(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnReal>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = real(s1.elem(i));
  return d;
}


// OTensor = Im(OTensor)
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnImag>::Type_t
imag(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnImag>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = imag(s1.elem(i));
  return d;
}


//! OTensor<T> = (OTensor<T> , OTensor<T>)
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OTensor<T2>, FnCmplx>::Type_t
cmplx(const OTensor<T1>& s1, const OTensor<T2>& s2)
{
  typename BinaryReturn<OTensor<T1>, OTensor<T2>, FnCmplx>::Type_t  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = cmplx(s1.elem(i), s2.elem(i));

  return d;
}


// ArcCos
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnArcCos>::Type_t
acos(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnArcCos>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = acos(s1.elem(i));
  return d;
}

// ArcSin
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnArcSin>::Type_t
asin(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnArcSin>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = asin(s1.elem(i));
  return d;
}

// ArcTan
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnArcTan>::Type_t
atan(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnArcTan>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = atan(s1.elem(i));
  return d;
}

// Cos
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnCos>::Type_t
cos(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnCos>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = cos(s1.elem(i));
  return d;
}

// Exp
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnExp>::Type_t
exp(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnExp>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = exp(s1.elem(i));
  return d;
}

// Fabs
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnFabs>::Type_t
fabs(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnFabs>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = fabs(s1.elem(i));
  return d;
}

// Log
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnLog>::Type_t
log(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnLog>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = log(s1.elem(i));
  return d;
}

// Sin
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnSin>::Type_t
sin(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnSin>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = sin(s1.elem(i));
  return d;
}

// Sqrt
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnSqrt>::Type_t
sqrt(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnSqrt>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = sqrt(s1.elem(i));
  return d;
}

// Tan
template<class T1>
inline typename UnaryReturn<OTensor<T1>, FnTan>::Type_t
tan(const OTensor<T1>& s1)
{
  typename UnaryReturn<OTensor<T1>, FnTan>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = tan(s1.elem(i));
  return d;
}


//! OTensor = pow(OTensor, OTensor)
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OTensor<T2>, FnPow>::Type_t
pow(const OTensor<T1>& s1, const OTensor<T2>& s2)
{
  typename BinaryReturn<OTensor<T1>, OTensor<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem(i));
  return d;
}

//! OTensor = pow(OTensor, OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OScalar<T2>, FnPow>::Type_t
pow(const OTensor<T1>& s1, const OScalar<T2>& s2)
{
  typename BinaryReturn<OTensor<T1>, OScalar<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem());
  return d;
}

//! OTensor = pow(OScalar, OTensor)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OTensor<T2>, FnPow>::Type_t
pow(const OScalar<T1>& s1, const OTensor<T2>& s2)
{
  typename BinaryReturn<OScalar<T1>, OTensor<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s2.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem(i));
  return d;
}


//! OTensor = atan2(OTensor, OTensor)
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OTensor<T2>, FnArcTan2>::Type_t
atan2(const OTensor<T1>& s1, const OTensor<T2>& s2)
{
  typename BinaryReturn<OTensor<T1>, OTensor<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = atan2(s1.elem(i), s2.elem(i));
  return d;
}

//! OTensor = atan2(OTensor, OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OScalar<T2>, FnArcTan2>::Type_t
atan2(const OTensor<T1>& s1, const OScalar<T2>& s2)
{
  typename BinaryReturn<OTensor<T1>, OScalar<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = atan2(s1.elem(i), s2.elem());
  return d;
}

//! OTensor = atan2(OScalar, OTensor)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OTensor<T2>, FnArcTan2>::Type_t
atan2(const OScalar<T1>& s1, const OTensor<T2>& s2)
{
  typename BinaryReturn<OScalar<T1>, OTensor<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s2.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = atan2(s1.elem(), s2.elem(i));
  return d;
}


//! OTensor = i * OTensor
template<class T>
inline typename UnaryReturn<OTensor<T>, FnTimesI>::Type_t
timesI(const OTensor<T>& s1)
{
  typename UnaryReturn<OTensor<T>, FnTimesI>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = timesI(s1.elem(i));

  return d;
}

//! OTensor = -i * OTensor
template<class T>
inline typename UnaryReturn<OTensor<T>, FnTimesMinusI>::Type_t
timesMinusI(const OTensor<T>& s1)
{
  typename UnaryReturn<OTensor<T>, FnTimesMinusI>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = timesMinusI(s1.elem(i));

  return d;
}


//! bool = isZero(OTensor)
template<class T>
bool
isZero(const OTensor<T>& s1)
{
  bool d = true;

  for(int i=0; i < s1.size(); ++i)
    d &= isZero(s1.elem(i));

  return d;
}


//! bool = isNaN(OTensor)
template<class T>
bool
isNaN(const OTensor<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.numElem(); ++i)
    d |= isNaN(s1.elem(i));

  return d;
}


//! bool = isInf(OTensor)
template<class T>
bool
isInf(const OTensor<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.numElem(); ++i)
    d |= isInf(s1.elem(i));

  return d;
}


//! bool = isFinite(OTensor)
template<class T>
bool
isFinite(const OTensor<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.numElem(); ++i)
    d |= isFinite(s1.elem(i));

  return d;
}


//! OTensor = outerProduct(OTensor, OTensor)
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OTensor<T2>, FnOuterProduct>::Type_t
outerProduct(const OTensor<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OTensor<T2>, FnOuterProduct>::Type_t  d;
  Array<int> nz(concat(l.size(), r.size()));
  d.resize(nz);
  const int nr = r.numElem();

  for(int i=0; i < l.numElem(); ++i)
    for(int j=0; j < r.numElem(); ++j)
      d.elem(i*nr+j) = l.elem(i) * r.elem(j);
  return d;
}

//! OTensor = outerProduct(OTensor, OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OTensor<T1>, OScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const OTensor<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OTensor<T1>, OScalar<T2>, FnOuterProduct>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem(i) * r.elem();
  return d;
}

//! OTensor = outerProduct(OScalar, OTensor)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OTensor<T2>, FnOuterProduct>::Type_t
outerProduct(const OScalar<T1>& l, const OTensor<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OTensor<T2>, FnOuterProduct>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = l.elem() * r.elem(i);
  return d;
}


//! Contract over a specific list of indices
template<class T>
struct UnaryReturn<OTensor<T>, FnContract> {
  typedef OTensor<typename UnaryReturn<T, FnContract>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OTensor<T>, FnContract>::Type_t
contract(const OTensor<T>& s1, const Array<int>& nn)
{
  typename UnaryReturn<OTensor<T>, FnContract>::Type_t  d;
  if (s1.size().size() == 0 || nn.size() == 0)
  {
    std::cerr << __func__ << ": Invalid OTensor size" << std::endl;
    exit(1);
  }
  int nd = s1.size().size() - nn.size();
  if (nd <= 0)
  {
    std::cerr << __func__ << ": Invalid OTensor contract" << std::endl;
    exit(1);
  }

  // Check indices
  for(int i=0; i < nn.size(); ++i)
  {
    if (nn[i] >= s1.size().size())
    {
      std::cerr << "Contract index out of bounds" << std::endl;
      exit(1);
    }
    if (s1.size()[nn[i]] != s1.size()[nn[0]])
    {
      std::cerr << "Contract indices not equal lengths" << std::endl;
      exit(1);
    }
  }

  Array<int> free_ind(nd);  // complement of mask

  // New array size
  Array<int> nz(nd);
  for(int i=0,j=0,k=0; i < s1.size().size(); ++i)
  {
    if (nn[j] == i)
    {
      j++;
    }
    else
    {
      free_ind[k] = i;
      nz[k++] = s1.size()[i];
    }
  }

  // Resize
  d.resize(nz);

  // Reset
  for(int i=0; i < d.numElem(); ++i)
    zero_rep(d.elem(i));

  // Do the contraction
  Array<int> ind(s1.size().size());

  for(int ipos0=0; ipos0 < s1.numElem(); ++ipos0)
  {
    bool diag = true;
    {
      // Decompose the indices of the source
      int ipos = ipos0;
      for(int i=0; i < s1.size().size(); ++i)
      {
	ind[i] = ipos % s1.size()[i];
	ipos /= s1.size()[i];
      }

      // Check if the summed indices are in fact the same
      for(int i=0; i < nn.size(); ++i)
      {
	if (ind[nn[i]] != ind[nn[0]])
	{
	  diag = false;
	  break;
	}
      }
    }
    
    // If this is a diagonal term, add it onto the result
    if (diag)
    {
      // Find the location inside the target
      int order = 0;
      for(int i=nz.size()-1; i >= 1; --i)
	order = nz[i-1]*(ind[free_ind[i]] + order);

      order += ind[free_ind[0]];

      d.elem(order) += s1.elem(ipos0);
    }
  }

  return d;
}


//! Extract observable tensor components 
/*! Generically, this is an identity operation. Defined differently under observable */
template<class T1>
struct UnaryReturn<OTensor<T1>, FnPeekObsTensor> {
  typedef OScalar<typename UnaryReturn<T1, FnPeekObsTensor>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OTensor<T>, FnPeekObsTensor>::Type_t
peekObs(const OTensor<T>& l, const Array<int>& nn)
{
  l.checkSize(__func__);

  return l[nn];
}

//! Extract observable tensor components 
/*! Generically, this is an identity operation. Defined differently under observable */
template<class T>
inline typename UnaryReturn<OTensor<T>, FnPeekObsTensor>::Type_t
peekObs(const OTensor<T>& l, int row)
{
  l.checkSize(__func__);

  return l.elem(row);
}

//! Extract observable vector components 
/*! Generically, this is an identity operation. Defined differently under observable */
template<class T1, class T2>
inline OTensor<T1>&
pokeObs(OTensor<T1>& l, const OScalar<T2>& r, const Array<int>& nn)
{
  l.checkSize(__func__);

  l[nn] = r.elem();
  return l;
}

//! Extract spin tensor components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<OTensor<T>, FnPeekSpinVector>::Type_t
peekSpin(const OTensor<T>& l, int row)
{
  typename UnaryReturn<OTensor<T>, FnPeekSpinVector>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = peekSpin(l.elem(i),row);
  return d;
}

//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<OTensor<T>, FnPeekSpinMatrix>::Type_t
peekSpin(const OTensor<T>& l, int row, int col)
{
  typename UnaryReturn<OTensor<T>, FnPeekSpinMatrix>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = peekSpin(l.elem(i),row,col);
  return d;
}

//! Insert spin tensor components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline typename UnaryReturn<OTensor<T1>, FnPokeSpinVector>::Type_t&
pokeSpin(OTensor<T1>& l, const OTensor<T2>& r, int row)
{
  typedef typename UnaryReturn<OTensor<T1>, FnPokeSpinVector>::Type_t  Return_t;
  l.checkSize(__func__, r.size());

  for(int i=0; i < l.size(); ++i)
    pokeSpin(l.elem(i),r.elem(i),row);
  return l;
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline typename UnaryReturn<OTensor<T1>, FnPokeSpinVector>::Type_t&
pokeSpin(OTensor<T1>& l, const OTensor<T2>& r, int row, int col)
{
  typedef typename UnaryReturn<OTensor<T1>, FnPokeSpinVector>::Type_t  Return_t;
  l.checkSize(__func__, r.size());

  for(int i=0; i < l.size(); ++i)
    pokeSpin(l.elem(i),r.elem(i),row,col);
  return l;
}


//! dest = 0
template<class T> 
inline void 
zero_rep(OTensor<T>& d) 
{
  d.checkSize(__func__);
  for(int i=0; i < d.numElem(); ++i)
    zero_rep(d.elem(i));
}

//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline void 
copymask(OTensor<T>& d, const OScalar<T1>& mask, const OTensor<T>& s1) 
{
  d.checkSize(__func__, s1.size());
  for(int i=0; i < d.numElem(); ++i)
    copymask(d.elem(i),mask.elem(),s1.elem(i));
}


//! dest  = random  
template<class T>
inline void
fill_random(OTensor<T>& d)
{
  d.checkSize(__func__);
  // Loop over rows the slowest
  for(int i=0; i < d.numElem(); ++i)
    fill_random(d.elem(i));
}


//! dest  = gaussian
template<class T>
inline void
fill_gaussian(OTensor<T>& d, OTensor<T>& r1, OTensor<T>& r2)
{
  d.checkResize(__func__, r1.size(), r2.size());
  for(int i=0; i < d.numElem(); ++i)
    fill_gaussian(d.elem(i), r1.elem(i), r2.elem(i));
}


#if 1
// Global sum over site indices only
template<class T>
struct UnaryReturn<OTensor<T>, FnSum > {
  typedef OTensor<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OTensor<T>, FnSum>::Type_t
sum(const OTensor<T>& s1)
{
  typename UnaryReturn<OTensor<T>, FnSum>::Type_t  d;
  d.checkSize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = sum(s1.elem(i));

  return d;
}
#endif


// OTensor<T> = localNorm2(OTensor<T>) = adj(OTensor<T>)*OTensor<T>)
template<class T>
struct UnaryReturn<OTensor<T>, FnNorm2 > {
  typedef OTensor<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<OTensor<T>, FnLocalNorm2 > {
  typedef OTensor<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OTensor<T>, FnLocalNorm2>::Type_t
localNorm2(const OTensor<T>& s1)
{
  typename UnaryReturn<OTensor<T>, FnLocalNorm2>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = localNorm2(s1.elem(i));

  return d;
}


//! OTensor<T> = InnerProduct(adj(OTensor<T1>)*OTensor<T1>)
template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OTensor<T2>, FnInnerProduct > {
  typedef OTensor<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OTensor<T2>, FnLocalInnerProduct > {
  typedef OTensor<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>
localInnerProduct(const OTensor<T1>& s1, const OTensor<T2>& s2)
{
  OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem() = localInnerProduct(s1.elem(i), s2.elem(i));

  return d;
}


//! OTensor<T> = InnerProductReal(adj(OTensor<T1>)*OTensor<T1>)
/*!
 * return  realpart of InnerProduct(adj(s1)*s2)
 */
template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OTensor<T2>, FnInnerProductReal > {
  typedef OTensor<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<OTensor<T1>, OTensor<T2>, FnLocalInnerProductReal > {
  typedef OTensor<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
inline OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>
localInnerProductReal(const OTensor<T1>& s1, const OTensor<T2>& s2)
{
  OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = localInnerProductReal(s1.elem(i), s2.elem(i));

  return d;
}


//! OTensor<T> = where(OScalar, OTensor, OTensor)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<OScalar<T1>, OTensor<T2>, OTensor<T3>, FnWhere> {
  typedef OTensor<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<OScalar<T1>, OTensor<T2>, OTensor<T3>, FnWhere>::Type_t
where(const OScalar<T1>& a, const OTensor<T2>& b, const OTensor<T3>& c)
{
  typename TrinaryReturn<OScalar<T1>, OTensor<T2>, OTensor<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, b.size(), c.size());

  // Not optimal - want to have where outside assignment
  for(int i=0; i < d.numElem(); ++i)
    d.elem(i) = where(a.elem(), b.elem(i), c.elem(i));

  return d;
}

/*! @} */  // end of group obstensor

} // namespace ENSEM

