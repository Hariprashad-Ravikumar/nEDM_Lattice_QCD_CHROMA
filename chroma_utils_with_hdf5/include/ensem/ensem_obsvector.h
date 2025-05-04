// -*- C++ -*-
// $Id: ensem_obsvector.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! \file
 * \brief Observable Vector
 */


namespace ENSEM {

//-------------------------------------------------------------------------------------
/*! \addtogroup obsvector Vector observable
 * \ingroup fiber
 *
 * Observable type that transforms like a vector
 *
 * @{
 */

//! Observable Vector class
/*!
 * All vector classes inherit this class
 * NOTE: For efficiency, there can be no virtual methods, so the data
 * portion is a part of the generic class, hence it is called a domain
 * and not a category
 */
template <class T> class OVector
{
public:
  OVector() {n1=0;F=0;}
  ~OVector() {delete[] F;}

  //---------------------------------------------------------
  //! Conversion constructor
  template<class T1>
  OVector(const OVector<T1>& rhs) : n1(0), F(0)
    {
      checkResize("convert OVector", rhs.size());
      for(int i=0; i < size(); ++i)
	elem(i) = rhs.elem(i);
    }

  //---------------------------------------------------------
  //! OVector = zero
  inline
  OVector& operator=(const Zero& rhs)
    {
      checkSize("OVector = zero");
      for(int i=0; i < size(); ++i)
	elem(i) = zero;
      return *this;
    }

  //! OVector = OScalar
  template<class T1>
  inline
  OVector& operator=(const OScalar<T1>& rhs) 
    {
      checkSize("OVector = const");
      for(int i=0; i < size(); ++i)
	elem(i) = rhs.elem();

      return *this;
    }

  //! OVector += OScalar
  template<class T1>
  inline
  OVector& operator+=(const OScalar<T1>& rhs) 
    {
      checkSize("OVector += OVector");
      for(int i=0; i < size(); ++i)
	elem(i) += rhs.elem();

      return *this;
    }

  //! OVector -= OScalar
  template<class T1>
  inline
  OVector& operator-=(const OScalar<T1>& rhs) 
    {
      checkSize("OVector -= OVector");
      for(int i=0; i < size(); ++i)
	elem(i) -= rhs.elem();

      return *this;
    }

  //! OVector *= OScalar
  template<class T1>
  inline
  OVector& operator*=(const OScalar<T1>& rhs) 
    {
      checkSize("OVector *= OVector");
      for(int i=0; i < size(); ++i)
	elem(i) *= rhs.elem();

      return *this;
    }

  //! OVector /= OScalar
  template<class T1>
  inline
  OVector& operator/=(const OScalar<T1>& rhs) 
    {
      checkSize("OVector /= OScalar");
      for(int i=0; i < size(); ++i)
	elem(i) /= rhs.elem();

      return *this;
    }

  //------------------------------------------------
  //! OVector = OVector
  /*! Set equal to another OVector */
  inline
  OVector& operator=(const OVector& rhs) 
    {
      checkResize("OVector = OVector", rhs.size());
      for(int i=0; i < size(); ++i)
	elem(i) = rhs.elem(i);

      return *this;
    }

  //! OVector = OVector
  /*! Set equal to another OVector */
  template<class T1>
  inline
  OVector& operator=(const OVector<T1>& rhs) 
    {
      checkResize("OVector = OVector<T1>", rhs.size());
      for(int i=0; i < size(); ++i)
	elem(i) = rhs.elem(i);

      return *this;
    }

  //! OVector += OVector
  template<class T1>
  inline
  OVector& operator+=(const OVector<T1>& rhs) 
    {
      checkSize("OVector += OVector", rhs.size());
      for(int i=0; i < size(); ++i)
	elem(i) += rhs.elem(i);

      return *this;
    }

  //! OVector -= OVector
  template<class T1>
  inline
  OVector& operator-=(const OVector<T1>& rhs) 
    {
      checkSize("OVector -= OVector", rhs.size());
      for(int i=0; i < size(); ++i)
	elem(i) -= rhs.elem(i);

      return *this;
    }

  //! OVector *= OVector
  template<class T1>
  inline
  OVector& operator*=(const OVector<T1>& rhs) 
    {
      checkSize("OVector *= OVector", rhs.size());
      for(int i=0; i < size(); ++i)
	elem(i) *= rhs.elem(i);

      return *this;
    }

  //! OVector /= OVector
  template<class T1>
  inline
  OVector& operator/=(const OVector<T1>& rhs) 
    {
      checkSize("OVector /= OVector", rhs.size());
      for(int i=0; i < size(); ++i)
	elem(i) /= rhs.elem(i);

      return *this;
    }


  //! Deep copy constructor
  OVector(const OVector& rhs) : n1(0), F(0)
    {
      checkResize("copy OVector", rhs.size());
      for(int i=0; i < size(); ++i)
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

      if (size() == 0)
      {
	std::cerr << s << ": Invalid OVector size" << std::endl;
	exit(1);
      }
    }

  inline void checkSize(const char *s, int n1) const
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkSize: %s, OVector[%d]\n", s, n1);
#endif

      if (size() == 0 && size() == n1)
      {
	std::cerr << s << ": Invalid OVector dest and/or source size" << std::endl;
	exit(1);
      }
    }

  inline void checkResize(const char *s, int n1)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, OVector[%d]\n", s, n1);
#endif

      if (n1 == 0)
      {
	std::cerr << "checkResize: " << s << ": invalid OVector source size" << std::endl;
	exit(1);
      }
      resize(n1);
    }

  inline void checkResize(const char *s, int n1, int n2)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, OVector[%d]\n", s, n1);
#endif

      if (n1 == 0 || n1 != n2)
      {
	std::cerr << "checkResize: " << s << ": invalid OVector source sizes" << std::endl;
	exit(1);
      }
      resize(n1);
    }

  inline void checkResize(const char *s, int n1, int n2, int n3)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, OVector[%d]\n", s, n1);
#endif

      if (n1 == 0 || n1 != n2 || n1 != n3)
      {
	std::cerr << "checkResize: " << s << ": invalid OVector source sizes" << std::endl;
	exit(1);
      }
      resize(n1);
    }

  //! Has this object been initialized (resized and such)
  bool initP() const 
    {
      return (n1 != 0) ? true : false;
    }

  //! Number of observables
  inline int size() const {return n1;}

  //! Number of observables
  inline int size1() const {return n1;}

  //! Number of observables
  inline int numElem() const {return n1;}

  //! Resize the number of observables
  inline void resize(int n)
    {
      if (n1 > 0)
      {
	delete[] F;
	n1 = 0;
      }
      n1 = n;
      F = new(std::nothrow) T[n1];
      if ( F == 0x0 ) { 
	std::cerr << "Unable to allocate memory in OVector::resize()" << std::endl;
	exit(1);
      }
    }

public:
  inline void resize(const OVector& a) {resize(a.n1);}

public:
  inline T& elem(int i) {return F[i];}
  inline const T& elem(int i) const {return F[i];}

private:
  int n1;
  T*  F;
};


// I/O
//! Binary input
template<class T>  
inline
void read(ADATIO::BinaryReader& bin, OVector<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.size(); ++i)
    read(bin, d.elem(i));
}

//! Binary output
template<class T>
inline
void write(ADATIO::BinaryWriter& bin, const OVector<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.size(); ++i)
    write(bin, d.elem(i));
}

//! Stream input
template<class T>  
inline
std::istream& operator>>(std::istream& s, OVector<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.size(); ++i)
    s >> d.elem(i);

  return s;
}

//! Stream output
template<class T>  
inline
std::ostream& operator<<(std::ostream& s, const OVector<T>& d)
{
  d.checkSize(__func__);

  for(int k = 0; k < d.size(); ++k)
  {
    s << d.elem(k);
    if (k < d.size()-1)
      s << "\n";
  }

  return s;
}


//! Text input
template<class T>  
inline
ENSEM::TextReader& operator>>(ENSEM::TextReader& txt, OVector<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.size(); ++i)
  {
    int k;
    txt >> k >> d.elem(i);
    if (k != i)
    {
      std::cerr << "error reading OVector" << std::endl;
      exit(1);
    }
  }

  return txt;
}

//! Text output
template<class T>  
inline
ENSEM::TextWriter& operator<<(ENSEM::TextWriter& txt, const OVector<T>& d)
{
  d.checkSize(__func__);
  typedef typename WordType<T>::Type_t Type_t;

  for(int k = 0; k < d.size(); ++k)
  {
    txt << k << " " << d.elem(k);
    if (k < d.size()-1)
      txt << "\n";
  }

  return txt;
}


//! XML output
template<class T> 
inline
ADATXML::XMLWriter& operator<<(ADATXML::XMLWriter& xml, const OVector<T>& d)
{
  xml.openTag("Vector");

  XMLWriterAPI::AttributeList alist;

  // Copy into another array first
  for(int i=0; i < d.size(); ++i)
  {
    alist.clear();
    alist.push_back(XMLWriterAPI::Attribute("row", i));

    xml.openTag("elem", alist);
    xml << d.elem(i);
    xml.closeTag();
  }

  xml.closeTag();  // Vector
  return xml;
}

/*! @} */  // end of group obsvector


//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1>
struct WordType<OVector<T1> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

// Internally used scalars
template<class T>
struct InternalScalar<OVector<T> > {
  typedef OScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive scalar leaving grid alone
template<class T>
struct PrimitiveScalar<OVector<T> > {
  typedef OScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving observable indices alone
template<class T>
struct EnsemScalar<OVector<T> > {
  typedef OVector<typename EnsemScalar<T>::Type_t>  Type_t;
};

// Traits class to label IO types
template<class T> 
struct EnsbcIO<OVector<T> > {
  enum {type = EnsbcIO<T>::type};
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(OVector) -> OVector
template<class T1, class Op>
struct UnaryReturn<OVector<T1>, Op> {
  typedef OVector<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};
// Default binary(OScalar,OVector) -> OVector
template<class T1, class T2, class Op>
struct BinaryReturn<OScalar<T1>, OVector<T2>, Op> {
  typedef OVector<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(OMatrix,OVector) -> OVector
template<class T1, class T2, class Op>
struct BinaryReturn<OMatrix<T1>, OVector<T2>, Op> {
  typedef OVector<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(OVector,OScalar) -> OVector
template<class T1, class T2, class Op>
struct BinaryReturn<OVector<T1>, OScalar<T2>, Op> {
  typedef OVector<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(OVector,OVector) -> OVector
template<class T1, class T2, class Op>
struct BinaryReturn<OVector<T1>, OVector<T2>, Op> {
  typedef OVector<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
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
struct BinaryReturn<OVector<T1>, OVector<T2>, OpAssign > {
  typedef OVector<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OVector<T1>, OVector<T2>, OpAddAssign > {
  typedef OVector<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OVector<T1>, OVector<T2>, OpSubtractAssign > {
  typedef OVector<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OVector<T1>, OScalar<T2>, OpMultiplyAssign > {
  typedef OVector<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OVector<T1>, OScalar<T2>, OpDivideAssign > {
  typedef OVector<T1> &Type_t;
};
 



//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

/*! \addtogroup obsvector */
/*! @{ */

// Observable Vectors

template<class T1>
inline typename UnaryReturn<OVector<T1>, OpUnaryPlus>::Type_t
operator+(const OVector<T1>& l)
{
  typename UnaryReturn<OVector<T1>, OpUnaryPlus>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = +l.elem(i);
  return d;
}


template<class T1>
inline typename UnaryReturn<OVector<T1>, OpUnaryMinus>::Type_t
operator-(const OVector<T1>& l)
{
  typename UnaryReturn<OVector<T1>, OpUnaryMinus>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = -l.elem(i);
  return d;
}


// OVector + OVector
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OVector<T2>, OpAdd>::Type_t
operator+(const OVector<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OVector<T1>, OVector<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) + r.elem(i);
  return d;
}

// OVector + OScalar
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OScalar<T2>, OpAdd>::Type_t
operator+(const OVector<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OVector<T1>, OScalar<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) + r.elem();
  return d;
}

// OScalar + OVector
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OVector<T2>, OpAdd>::Type_t
operator+(const OScalar<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OVector<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() + r.elem(i);
  return d;
}


// OVector - OVector
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OVector<T2>, OpSubtract>::Type_t
operator-(const OVector<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OVector<T1>, OVector<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) - r.elem(i);
  return d;
}

// OVector - OScalar
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OScalar<T2>, OpSubtract>::Type_t
operator-(const OVector<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OVector<T1>, OScalar<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) - r.elem();
  return d;
}

// OScalar - OVector
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OVector<T2>, OpSubtract>::Type_t
operator-(const OScalar<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OVector<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() - r.elem(i);
  return d;
}


// OVector * OVector
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OVector<T2>, OpMultiply>::Type_t
operator*(const OVector<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OVector<T1>, OVector<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) * r.elem(i);
  return d;
}

// OScalar * OVector
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OVector<T2>, OpMultiply>::Type_t
operator*(const OScalar<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OVector<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() * r.elem(i);
  return d;
}

// OVector * OScalar
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OScalar<T2>, OpMultiply>::Type_t
operator*(const OVector<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OVector<T1>, OScalar<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) * r.elem();
  return d;
}


// OMatrix * OVector
template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OVector<T2>, OpMultiply>::Type_t
operator*(const OMatrix<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OVector<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, l.size1());

  for(int i=0; i < l.size1(); ++i)
  {
    d.elem(i) = l.elem(i,0) * r.elem(0);
    for(int j=1; j < l.size2(); ++j)
      d.elem(i) += l.elem(i,j) * r.elem(j);
  }

  return d;
}


//! OVector / OVector
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OVector<T2>, OpDivide>::Type_t
operator/(const OVector<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OVector<T1>, OVector<T2>, OpDivide>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) / r.elem(i);
  return d;
}

//! OVector / OScalar
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OScalar<T2>, OpDivide>::Type_t
operator/(const OVector<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OVector<T1>, OScalar<T2>, OpDivide>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) / r.elem();
  return d;
}

//! OScalar / OVector
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OVector<T2>, OpDivide>::Type_t
operator/(const OScalar<T1>& l, const OVector<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OVector<T2>, OpDivide>::Type_t  d;
  d.checkResize(__func__, l.size(), r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() / r.elem(i);
  return d;
}



//-----------------------------------------------------------------------------
// Functions


// OVector = adj(OVector)
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnAdjoint>::Type_t
adj(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnAdjoint>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = adj(s1.elem(i));
  return d;
}


// OVector = conj(OVector)
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnConjugate>::Type_t
conj(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnConjugate>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = conj(s1.elem(i));
  return d;
}


// OVector = transpose(OVector)
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnTranspose>::Type_t
transpose(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnTranspose>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = transpose(s1.elem(i));
  return d;
}


// OVector = Trace(OVector)
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnTrace>::Type_t
trace(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnTrace>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = trace(s1.elem(i));
  return d;
}


// OVector = Re(Trace(OVector))
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnRealTrace>::Type_t
trace_real(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnRealTrace>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = trace_real(s1.elem(i));
  return d;
}


// OVector = Im(Trace(OVector))
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnImagTrace>::Type_t
trace_imag(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnImagTrace>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = trace_imag(s1.elem(i));
  return d;
}


// OVector = Re(OVector)
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnReal>::Type_t
real(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnReal>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = real(s1.elem(i));
  return d;
}


// OVector = Im(OVector)
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnImag>::Type_t
imag(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnImag>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = imag(s1.elem(i));
  return d;
}


//! OVector<T> = (OVector<T> , OVector<T>)
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OVector<T2>, FnCmplx>::Type_t
cmplx(const OVector<T1>& s1, const OVector<T2>& s2)
{
  typename BinaryReturn<OVector<T1>, OVector<T2>, FnCmplx>::Type_t  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = cmplx(s1.elem(i), s2.elem(i));

  return d;
}


// ArcCos
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnArcCos>::Type_t
acos(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnArcCos>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = acos(s1.elem(i));
  return d;
}

// ArcSin
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnArcSin>::Type_t
asin(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnArcSin>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = asin(s1.elem(i));
  return d;
}

// ArcTan
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnArcTan>::Type_t
atan(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnArcTan>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = atan(s1.elem(i));
  return d;
}

// Cos
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnCos>::Type_t
cos(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnCos>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = cos(s1.elem(i));
  return d;
}

// Exp
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnExp>::Type_t
exp(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnExp>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = exp(s1.elem(i));
  return d;
}

// Fabs
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnFabs>::Type_t
fabs(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnFabs>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = fabs(s1.elem(i));
  return d;
}

// Log
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnLog>::Type_t
log(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnLog>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = log(s1.elem(i));
  return d;
}

// Sin
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnSin>::Type_t
sin(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnSin>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = sin(s1.elem(i));
  return d;
}

// Sqrt
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnSqrt>::Type_t
sqrt(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnSqrt>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = sqrt(s1.elem(i));
  return d;
}

// Tan
template<class T1>
inline typename UnaryReturn<OVector<T1>, FnTan>::Type_t
tan(const OVector<T1>& s1)
{
  typename UnaryReturn<OVector<T1>, FnTan>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = tan(s1.elem(i));
  return d;
}


//! OVector = pow(OVector, OVector)
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OVector<T2>, FnPow>::Type_t
pow(const OVector<T1>& s1, const OVector<T2>& s2)
{
  typename BinaryReturn<OVector<T1>, OVector<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem(i));
  return d;
}

//! OVector = pow(OVector, OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OScalar<T2>, FnPow>::Type_t
pow(const OVector<T1>& s1, const OScalar<T2>& s2)
{
  typename BinaryReturn<OVector<T1>, OScalar<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem());
  return d;
}

//! OVector = pow(OScalar, OVector)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OVector<T2>, FnPow>::Type_t
pow(const OScalar<T1>& s1, const OVector<T2>& s2)
{
  typename BinaryReturn<OScalar<T1>, OVector<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s2.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem(i));
  return d;
}


//! OVector = atan2(OVector, OVector)
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OVector<T2>, FnArcTan2>::Type_t
atan2(const OVector<T1>& s1, const OVector<T2>& s2)
{
  typename BinaryReturn<OVector<T1>, OVector<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = atan2(s1.elem(i), s2.elem(i));
  return d;
}

//! OVector = atan2(OVector, OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OVector<T1>, OScalar<T2>, FnArcTan2>::Type_t
atan2(const OVector<T1>& s1, const OScalar<T2>& s2)
{
  typename BinaryReturn<OVector<T1>, OScalar<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = atan2(s1.elem(i), s2.elem());
  return d;
}

//! OVector = atan2(OScalar, OVector)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OVector<T2>, FnArcTan2>::Type_t
atan2(const OScalar<T1>& s1, const OVector<T2>& s2)
{
  typename BinaryReturn<OScalar<T1>, OVector<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s2.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = atan2(s1.elem(), s2.elem(i));
  return d;
}


//! OVector = i * OVector
template<class T>
inline typename UnaryReturn<OVector<T>, FnTimesI>::Type_t
timesI(const OVector<T>& s1)
{
  typename UnaryReturn<OVector<T>, FnTimesI>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = timesI(s1.elem(i));

  return d;
}

//! OVector = -i * OVector
template<class T>
inline typename UnaryReturn<OVector<T>, FnTimesMinusI>::Type_t
timesMinusI(const OVector<T>& s1)
{
  typename UnaryReturn<OVector<T>, FnTimesMinusI>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = timesMinusI(s1.elem(i));

  return d;
}


//! bool = isZero(OVector)
template<class T>
bool
isZero(const OVector<T>& s1)
{
  bool d = true;

  for(int i=0; i < s1.size(); ++i)
    d &= isZero(s1.elem(i));

  return d;
}


//! bool = isNaN(OVector)
template<class T>
bool
isNaN(const OVector<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.size(); ++i)
    d |= isNaN(s1.elem(i));

  return d;
}


//! bool = isInf(OVector)
template<class T>
bool
isInf(const OVector<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.size(); ++i)
    d |= isInf(s1.elem(i));

  return d;
}


//! bool = isFinite(OVector)
template<class T>
bool
isFinite(const OVector<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.size(); ++i)
    d |= isFinite(s1.elem(i));

  return d;
}


// OVector = shift(OVector, int)
template<class T1>
struct UnaryReturn<OVector<T1>, FnShift> {
  typedef OVector<typename UnaryReturn<T1, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OVector<T1>, FnShift>::Type_t
shift(const OVector<T1>& s1, int sh)
{
  typename UnaryReturn<OVector<T1>, FnShift>::Type_t  d;
  d.checkResize(__func__, s1.size());
  int len = d.size();

  for(int i=0; i < len; ++i)
  {
    int ii = (i + len + sh) % len;
    d.elem(i) = s1.elem(ii);
  }

  // Clean out the ends
  if (sh > 0)
  {
    for(int i = len-sh; i < len; ++i)
    {
      zero_rep(d.elem(i));
    }
  }
  else if (sh < 0)
  {
    for(int i = 0; i < -sh; ++i)
    {
      zero_rep(d.elem(i));
    }
  }

  return d;
}


// OVector = cshift(OVector, int)
template<class T1>
struct UnaryReturn<OVector<T1>, FnCshift> {
  typedef OVector<typename UnaryReturn<T1, FnCshift>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OVector<T1>, FnCshift>::Type_t
cshift(const OVector<T1>& s1, int sh)
{
  typename UnaryReturn<OVector<T1>, FnCshift>::Type_t  d;
  d.checkResize(__func__, s1.size());
  int len = d.size();

  for(int i=0; i < len; ++i)
  {
    int ii = (i + len + sh) % len;
    d.elem(i) = s1.elem(ii);
  }

  return d;
}


//! Extract observable vector components 
/*! Generically, this is an identity operation. Defined differently under observable */
template<class T1>
struct UnaryReturn<OVector<T1>, FnPeekObsVector> {
  typedef OScalar<typename UnaryReturn<T1, FnPeekObsVector>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OVector<T>, FnPeekObsVector>::Type_t
peekObs(const OVector<T>& l, int row)
{
  l.checkSize(__func__);

  return l.elem(row);
}

//! Extract spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<OVector<T>, FnPeekSpinVector>::Type_t
peekSpin(const OVector<T>& l, int row)
{
  typename UnaryReturn<OVector<T>, FnPeekSpinVector>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = peekSpin(l.elem(i),row);
  return d;
}

//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<OVector<T>, FnPeekSpinMatrix>::Type_t
peekSpin(const OVector<T>& l, int row, int col)
{
  typename UnaryReturn<OVector<T>, FnPeekSpinMatrix>::Type_t  d;
  d.checkResize(__func__, l.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = peekSpin(l.elem(i),row,col);
  return d;
}

//! Extract observable vector components 
/*! Generically, this is an identity operation. Defined differently under observable */
template<class T1, class T2>
inline OVector<T1>&
pokeObs(OVector<T1>& l, const OScalar<T2>& r, int row)
{
  l.checkSize(__func__);

  l.elem(row) = r.elem();
  return l;
}

//! Insert spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline typename UnaryReturn<OVector<T1>, FnPokeSpinVector>::Type_t&
pokeSpin(OVector<T1>& l, const OVector<T2>& r, int row)
{
  l.checkSize(__func__, r.size());

  for(int i=0; i < l.size(); ++i)
    pokeSpin(l.elem(i),r.elem(i),row);
  return l;
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline typename UnaryReturn<OVector<T1>, FnPokeSpinVector>::Type_t&
pokeSpin(OVector<T1>& l, const OVector<T2>& r, int row, int col)
{
  l.checkSize(__func__, r.size());

  for(int i=0; i < l.size(); ++i)
    pokeSpin(l.elem(i),r.elem(i),row,col);
  return l;
}


//-----------------------------------------------------------------------------
// Gamma algebra
template<int N, int m, class T, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, OVector<T>, OpGammaConstMultiply> {
  typedef OVector<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! OVector = Gamma<N,m> * OVector
template<class T, int N, int m>
inline typename BinaryReturn<GammaConst<N,m>, OVector<T>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<N,m>& l, const OVector<T>& r)
{
  typename BinaryReturn<GammaConst<N,m>, OVector<T>, OpGammaConstMultiply>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l * r.elem(i);
  return d;
}

template<class T, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<OVector<T>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef OVector<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! OVector = OVector * Gamma<N,m>
template<class T, int N, int m>
inline typename BinaryReturn<OVector<T>, GammaConst<N,m>, OpMultiplyGammaConst>::Type_t
operator*(const OVector<T>& l, const GammaConst<N,m>& r)
{
  typename BinaryReturn<OVector<T>, GammaConst<N,m>, OpMultiplyGammaConst>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) * r;
  return d;
}


//-----------------------------------------------------------------------------
// Gamma algebra
template<int N, int m, class T, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, OVector<T>, OpGammaConstDPMultiply> {
  typedef OVector<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! OVector = Gamma<N,m> * OVector
template<class T, int N, int m>
inline typename BinaryReturn<GammaConstDP<N,m>, OVector<T>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<N,m>& l, const OVector<T>& r)
{
  typename BinaryReturn<GammaConstDP<N,m>, OVector<T>, OpGammaConstDPMultiply>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l * r.elem(i);
  return d;
}

template<class T, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<OVector<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef OVector<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! OVector = OVector * Gamma<N,m>
template<class T, int N, int m>
inline typename BinaryReturn<OVector<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP>::Type_t
operator*(const OVector<T>& l, const GammaConstDP<N,m>& r)
{
  typename BinaryReturn<OVector<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP>::Type_t  d;
  d.checkResize(__func__, r.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) * r;
  return d;
}


//-----------------------------------------------------------------------------
//! dest = 0
template<class T> 
inline void 
zero_rep(OVector<T>& d) 
{
  d.checkSize(__func__);
  for(int i=0; i < d.size(); ++i)
    zero_rep(d.elem(i));
}

//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline void 
copymask(OVector<T>& d, const OScalar<T1>& mask, const OVector<T>& s1) 
{
  d.checkSize(__func__, s1.size());
  for(int i=0; i < d.size(); ++i)
    copymask(d.elem(i),mask.elem(),s1.elem(i));
}


//! dest  = random  
template<class T>
inline void
fill_random(OVector<T>& d)
{
  d.checkSize(__func__);
  // Loop over rows the slowest
  for(int i=0; i < d.size(); ++i)
    fill_random(d.elem(i));
}


//! dest  = gaussian
template<class T>
inline void
fill_gaussian(OVector<T>& d, OVector<T>& r1, OVector<T>& r2)
{
  d.checkResize(__func__, r1.size(), r2.size());
  for(int i=0; i < d.size(); ++i)
    fill_gaussian(d.elem(i), r1.elem(i), r2.elem(i));
}


#if 1
// Global sum over site indices only
template<class T>
struct UnaryReturn<OVector<T>, FnSum > {
  typedef OVector<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OVector<T>, FnSum>::Type_t
sum(const OVector<T>& s1)
{
  typename UnaryReturn<OVector<T>, FnSum>::Type_t  d;
  d.checkSize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = sum(s1.elem(i));

  return d;
}
#endif


// OVector<T> = localNorm2(OVector<T>) = adj(OVector<T>)*OVector<T>)
template<class T>
struct UnaryReturn<OVector<T>, FnNorm2 > {
  typedef OVector<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<OVector<T>, FnLocalNorm2 > {
  typedef OVector<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OVector<T>, FnLocalNorm2>::Type_t
localNorm2(const OVector<T>& s1)
{
  typename UnaryReturn<OVector<T>, FnLocalNorm2>::Type_t  d;
  d.checkResize(__func__, s1.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = localNorm2(s1.elem(i));

  return d;
}


//! OVector<T> = InnerProduct(adj(OVector<T1>)*OVector<T1>)
template<class T1, class T2>
struct BinaryReturn<OVector<T1>, OVector<T2>, FnInnerProduct > {
  typedef OVector<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<OVector<T1>, OVector<T2>, FnLocalInnerProduct > {
  typedef OVector<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>
localInnerProduct(const OVector<T1>& s1, const OVector<T2>& s2)
{
  OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.size(); ++i)
    d.elem() = localInnerProduct(s1.elem(i), s2.elem(i));

  return d;
}


//! OVector<T> = InnerProductReal(adj(OVector<T1>)*OVector<T1>)
/*!
 * return  realpart of InnerProduct(adj(s1)*s2)
 */
template<class T1, class T2>
struct BinaryReturn<OVector<T1>, OVector<T2>, FnInnerProductReal > {
  typedef OVector<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<OVector<T1>, OVector<T2>, FnLocalInnerProductReal > {
  typedef OVector<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
inline OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>
localInnerProductReal(const OVector<T1>& s1, const OVector<T2>& s2)
{
  OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  d;
  d.checkResize(__func__, s1.size(), s2.size());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = localInnerProductReal(s1.elem(i), s2.elem(i));

  return d;
}


//! OVector<T> = where(OScalar, OVector, OVector)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<OScalar<T1>, OVector<T2>, OVector<T3>, FnWhere> {
  typedef OVector<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<OScalar<T1>, OVector<T2>, OVector<T3>, FnWhere>::Type_t
where(const OScalar<T1>& a, const OVector<T2>& b, const OVector<T3>& c)
{
  typename TrinaryReturn<OScalar<T1>, OVector<T2>, OVector<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, b.size(), c.size());

  // Not optimal - want to have where outside assignment
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = where(a.elem(), b.elem(i), c.elem(i));

  return d;
}

/*! @} */  // end of group obsvector

} // namespace ENSEM

