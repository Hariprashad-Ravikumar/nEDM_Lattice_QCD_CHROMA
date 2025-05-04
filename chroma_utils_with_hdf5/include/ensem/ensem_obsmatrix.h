// -*- C++ -*-
// $Id: ensem_obsmatrix.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! \file
 * \brief Observable Matrix
 */

namespace ENSEM {

//-------------------------------------------------------------------------------------
/*! \addtogroup obsmatrix Matrix primitive
 * \ingroup fiber
 *
 * Observable type that transforms like a matrix
 *
 * @{
 */


//! Observable Matrix class
/*!
 * All Matrix classes inherit this class
 * NOTE: For efficiency, there can be no virtual methods, so the data
 * portion is a part of the generic class, hence it is called a domain
 * and not a category
 */
template <class T> class OMatrix
{
public:
  OMatrix() {}
  ~OMatrix() {}

  //! OMatrix = zero
  inline
  OMatrix& assign(const Zero& rhs)
    {
      checkSize("OMatrix = zero");
      for(int i=0; i < N; ++i)
	for(int j=0; j < N; ++j)
	  zero_rep(elem(i,j));

      return *this;
    }

  //! OMatrix = OScalar
  /*! Fill with primitive scalar */
  template<class T1>
  inline
  OMatrix& assign(const OScalar<T1>& rhs)
    {
      for(int i=0; i < N; ++i)
	for(int j=0; j < N; ++j)
	  if (i == j)
	    elem(i,j) = rhs.elem();
	  else
	    zero_rep(elem(i,j));

      return *this;
    }

  //! OMatrix = OMatrix
  /*! Set equal to another OMatrix */
  template<class T1>
  inline
  OMatrix& assign(const OMatrix<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	for(int j=0; j < N; ++j)
	  elem(i,j) = rhs.elem(i,j);

      return *this;
    }

  //! OMatrix += OMatrix
  template<class T1>
  inline
  OMatrix& operator+=(const OMatrix<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	for(int j=0; j < N; ++j)
	  elem(i,j) += rhs.elem(i,j);

      return *this;
    }

  //! OMatrix -= OMatrix
  template<class T1>
  inline
  OMatrix& operator-=(const OMatrix<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	for(int j=0; j < N; ++j)
	  elem(i,j) -= rhs.elem(i,j);

      return *this;
    }

  //! OMatrix += OScalar
  template<class T1>
  inline
  OMatrix& operator+=(const OScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i,i) += rhs.elem();

      return *this;
    }

  //! OMatrix -= OScalar
  template<class T1>
  inline
  OMatrix& operator-=(const OScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i,i) -= rhs.elem();

      return *this;
    }

  //! OMatrix *= OScalar
  template<class T1>
  inline
  OMatrix& operator*=(const OScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	for(int j=0; j < N; ++j)
	  elem(i,j) *= rhs.elem();

      return *this;
    }

  //! OMatrix /= OScalar
  template<class T1>
  inline
  OMatrix& operator/=(const OScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	for(int j=0; j < N; ++j)
	  elem(i,j) /= rhs.elem();

      return *this;
    }

  //! Deep copy constructor
  OMatrix(const OMatrix& a) : n1(0), F(0)
    { 
      checkResize("copy OMatrix", a.size());
      for(int i=0; i < N*N; ++i)
	F[i] = a.F[i];
    }

public:
  //---------------------------------------------------------
  inline void checkSize(const char *s) const
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkSize: %s\n", s);
#endif

      if (size() == 0)
      {
	cerr << s << ": Invalid Ensem size" << endl;
	exit(1);
      }
    }

  inline void checkSize(const char *s, int n1) const
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkSize: %s, Ensem[%d]\n", s, n1);
#endif

      if (size() == 0 && size() == n1)
      {
	cerr << s << ": Invalid Ensem dest and/or source size" << endl;
	exit(1);
      }
    }

  inline void checkResize(const char *s, int n1)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, Ensem[%d]\n", s, n1);
#endif

      if (n1 == 0)
      {
	cerr << "checkResize: " << s << ": invalid Ensem source size" << endl;
	exit(1);
      }
      resize(n1);
    }

  inline void checkResize(const char *s, int n1, int n2)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, Ensem[%d]\n", s, n1);
#endif

      if (n1 == 0 || n1 != n2)
      {
	cerr << "checkResize: " << s << ": invalid Ensem source sizes" << endl;
	exit(1);
      }
      resize(n1);
    }

  inline void checkResize(const char *s, int n1, int n2, int n3)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, Ensem[%d]\n", s, n1);
#endif

      if (n1 == 0 || n1 != n2 || n1 != n3)
      {
	cerr << "checkResize: " << s << ": invalid Ensem source sizes" << endl;
	exit(1);
      }
      resize(n1);
    }

  //! Number of configurations
  inline int size() const {return n1;}

  //! Size of observable bits
  inline int size1() const {return n1;}

  //! Resize the number of configs
  inline void resize(int n)
    {
      if (n1 > 0)
      {
	delete[] F;
	n1 = 0;
      }
      n1 = n;
      F = new(nothrow) T[n1];
      if ( F == 0x0 ) { 
	cerr << "Unable to allocate memory in OVector::resize()" << endl;
	exit(1);
      }
    }

public:
  T& elem(int i, int j) {return F[j+N*i];}
  const T& elem(int i, int j) const {return F[j+N*i];}

private:
  int n1;
  T* F;
};


//! Text input
template<class T>  
inline
TextReader& operator>>(TextReader& txt, OMatrix<T>& d)
{
  for(int j=0; j < N; ++j)
    for(int i=0; i < N; ++i)
      txt >> d.elem(i,j);

  return txt;
}

//! Text output
template<class T>  
inline
TextWriter& operator<<(TextWriter& txt, const OMatrix<T>& d)
{
  for(int j=0; j < N; ++j)
    for(int i=0; i < N; ++i)
      txt << d.elem(i,j);

  return txt;
}


//! XML output
template<class T>  
inline
ADAT::XMLWriter& operator<<(ADAT::XMLWriter& xml, const OMatrix<T>& d)
{
  xml.openTag("Matrix");

  XMLWriterAPI::AttributeList alist;

  for(int i=0; i < N; ++i)
  {
    for(int j=0; j < N; ++j)
    {
      alist.clear();
      alist.push_back(XMLWriterAPI::Attribute("row", i));
      alist.push_back(XMLWriterAPI::Attribute("col", j));

      xml.openTag("elem", alist);
      xml << d.elem(i,j);
      xml.closeTag();
    }
  }

  xml.closeTag(); // Matrix
  return xml;
}

/*! @} */  // end of group obsmatrix

//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1>
struct WordType<OMatrix<T1> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

// Internally used scalars
template<class T>
struct InternalScalar<OMatrix<T> > {
  typedef OScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive scalar leaving grid alone
template<class T>
struct PrimitiveScalar<OMatrix<T> > {
  typedef OScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T>
struct EnsemScalar<OMatrix<T> > {
  typedef OMatrix<typename EnsemScalar<T>::Type_t>  Type_t;
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

/*
 * NOTE***: no Op defaults - they cause conflicts with specialized versions.
 * Avoid them.
 */


#if 0
template<class T1, class T2>
struct UnaryReturn<OScalar<T2>, OpCast<T1> > {
  typedef OScalar<typename UnaryReturn<T, OpCast>::Type_t>  Type_t;
//  typedef T1 Type_t;
};
#endif

template<class T>
struct UnaryReturn<OMatrix<T>, OpIdentity> {
  typedef OMatrix<typename UnaryReturn<T, OpIdentity>::Type_t>  Type_t;
};


// Assignment is different
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpAssign > {
  typedef OMatrix<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpAddAssign > {
  typedef OMatrix<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpSubtractAssign > {
  typedef OMatrix<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpAssign > {
  typedef OMatrix<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpAddAssign > {
  typedef OMatrix<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpSubtractAssign > {
  typedef OMatrix<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpMultiplyAssign > {
  typedef OMatrix<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpDivideAssign > {
  typedef OMatrix<T1> &Type_t;
};
 


//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------
/*! \addtogroup obsmatrix */
/*! @{ */

// Observable Matrices

// OMatrix = + OMatrix
template<class T>
struct UnaryReturn<OMatrix<T>, OpUnaryPlus> {
  typedef OMatrix<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OMatrix<T1>, OpUnaryPlus>::Type_t
operator+(const OMatrix<T1>& l)
{
  typename UnaryReturn<OMatrix<T1>, OpUnaryPlus>::Type_t  d;
  d.resize(l.size1(), l.size2());

  for(int i=0; i < d.size1(); ++i)
    for(int j=0; j < d.size2(); ++j)
      d.elem(i,j) = +l.elem(i,j);

  return d;
}


// OMatrix = - OMatrix
template<class T>
struct UnaryReturn<OMatrix<T>, OpUnaryMinus> {
  typedef OMatrix<typename UnaryReturn<T, OpUnaryMinus>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OMatrix<T1>, OpUnaryMinus>::Type_t
operator-(const OMatrix<T1>& l)
{
  typename UnaryReturn<OMatrix<T1>, OpUnaryMinus>::Type_t  d;
  d.resize(l.size1(), l.size2());

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = -l.elem(i,j);

  return d;
}


// OMatrix = OMatrix + OMatrix
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpAdd> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpAdd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpAdd>::Type_t
operator+(const OMatrix<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpAdd>::Type_t  d;
  d.resize(l.size1(), l.size2());

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = l.elem(i,j) + r.elem(i,j);

  return d;
}

// OMatrix = OMatrix + OScalar
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpAdd> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpAdd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OScalar<T2>, OpAdd>::Type_t
operator+(const OMatrix<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OScalar<T2>, OpAdd>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = (i == j) ? l.elem(i,i) + r.elem() : l.elem(i,j);

  return d;
}

// OMatrix = OScalar + OMatrix
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OMatrix<T2>, OpAdd> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpAdd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OMatrix<T2>, OpAdd>::Type_t
operator+(const OScalar<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OMatrix<T2>, OpAdd>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = (i == j) ? l.elem() + r.elem(i,i) : r.elem(i,j);

  return d;
}


// OMatrix = OMatrix - OMatrix
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpSubtract> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpSubtract>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpSubtract>::Type_t
operator-(const OMatrix<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpSubtract>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = l.elem(i,j) - r.elem(i,j);

  return d;
}

// OMatrix = OMatrix - OScalar
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpSubtract> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpSubtract>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OScalar<T2>, OpSubtract>::Type_t
operator-(const OMatrix<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OScalar<T2>, OpSubtract>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = (i == j) ? l.elem(i,i) - r.elem() : l.elem(i,j);

  return d;
}

// OMatrix = OScalar - OMatrix
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OMatrix<T2>, OpSubtract> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpSubtract>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OMatrix<T2>, OpSubtract>::Type_t
operator-(const OScalar<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OMatrix<T2>, OpSubtract>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = (i == j) ? l.elem() - r.elem(i,i) : r.elem(i,j);

  return d;
}


// OMatrix = OMatrix * OMatrix
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpMultiply> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpMultiply>::Type_t
operator*(const OMatrix<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
    {
      d.elem(i,j) = l.elem(i,0) * r.elem(0,j);
      for(int k=1; k < N; ++k)
	d.elem(i,j) += l.elem(i,k) * r.elem(k,j);
    }

  return d;
}

// OMatrix = OMatrix * OScalar
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpMultiply> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OScalar<T2>, OpMultiply>::Type_t
operator*(const OMatrix<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OScalar<T2>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = l.elem(i,j) * r.elem();
  return d;
}

// OMatrix = OScalar * OMatrix
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OMatrix<T2>, OpMultiply> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OMatrix<T2>, OpMultiply>::Type_t
operator*(const OScalar<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OMatrix<T2>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = l.elem() * r.elem(i,j);
  return d;
}


// OMatrix = OMatrix / OScalar
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, OpDivide> {
  typedef OMatrix<typename BinaryReturn<T1, T2, OpDivide>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OScalar<T2>, OpDivide>::Type_t
operator/(const OMatrix<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OScalar<T2>, OpDivide>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = l.elem(i,j) / r.elem();
  return d;
}



//-----------------------------------------------------------------------------
// Functions

// Adjoint
template<class T>
struct UnaryReturn<OMatrix<T>, FnAdjoint> {
  typedef OMatrix<typename UnaryReturn<T, FnAdjoint>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OMatrix<T1>, FnAdjoint>::Type_t
adj(const OMatrix<T1>& l)
{
  typename UnaryReturn<OMatrix<T1>, FnAdjoint>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = adj(l.elem(j,i));

  return d;
}


// Conjugate
template<class T>
struct UnaryReturn<OMatrix<T>, FnConjugate> {
  typedef OMatrix<typename UnaryReturn<T, FnConjugate>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OMatrix<T1>, FnConjugate>::Type_t
conj(const OMatrix<T1>& l)
{
  typename UnaryReturn<OMatrix<T1>, FnConjugate>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = conj(l.elem(i,j));

  return d;
}


// Transpose
template<class T>
struct UnaryReturn<OMatrix<T>, FnTranspose> {
  typedef OMatrix<typename UnaryReturn<T, FnTranspose>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OMatrix<T1>, FnTranspose>::Type_t
transpose(const OMatrix<T1>& l)
{
  typename UnaryReturn<OMatrix<T1>, FnTranspose>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = transpose(l.elem(j,i));

  return d;
}


// TRACE
// OScalar = Trace(OMatrix)
template<class T>
struct UnaryReturn<OMatrix<T>, FnTrace> {
  typedef OScalar<typename UnaryReturn<T, FnTrace>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnTrace>::Type_t
trace(const OMatrix<T>& s1)
{
  typename UnaryReturn<OMatrix<T>, FnTrace>::Type_t  d;

  d.elem() = trace(s1.elem(0,0));
  for(int i=1; i < N; ++i)
    d.elem() += trace(s1.elem(i,i));

  return d;
}


// OScalar = Re(Trace(OMatrix))
template<class T>
struct UnaryReturn<OMatrix<T>, FnRealTrace> {
  typedef OScalar<typename UnaryReturn<T, FnRealTrace>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OMatrix<T1>, FnRealTrace>::Type_t
realTrace(const OMatrix<T1>& s1)
{
  typename UnaryReturn<OMatrix<T1>, FnRealTrace>::Type_t  d;

  d.elem() = realTrace(s1.elem(0,0));
  for(int i=1; i < N; ++i)
    d.elem() += realTrace(s1.elem(i,i));

  return d;
}


//! OScalar = Im(Trace(OMatrix))
template<class T>
struct UnaryReturn<OMatrix<T>, FnImagTrace> {
  typedef OScalar<typename UnaryReturn<T, FnImagTrace>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OMatrix<T1>, FnImagTrace>::Type_t
imagTrace(const OMatrix<T1>& s1)
{
  typename UnaryReturn<OMatrix<T1>, FnImagTrace>::Type_t  d;

  d.elem() = imagTrace(s1.elem(0,0));
  for(int i=1; i < N; ++i)
    d.elem() += imagTrace(s1.elem(i,i));

  return d;
}


//! OMatrix = traceSpin(OMatrix)   [this is an identity in general]
template<class T>
struct UnaryReturn<OMatrix<T>, FnTraceSpin> {
  typedef OMatrix<typename UnaryReturn<T, FnTraceSpin>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnTraceSpin>::Type_t
traceSpin(const OMatrix<T>& s1)
{
  typename UnaryReturn<OMatrix<T>, FnTraceSpin>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = traceSpin(s1.elem(i,j));

  return d;
}


// OScalar = traceMultiply(OMatrix,OMatrix)
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnTraceMultiply> {
  typedef OScalar<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnTraceMultiply>::Type_t
traceMultiply(const OMatrix<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnTraceMultiply>::Type_t  d;

  d.elem() = traceMultiply(l.elem(0,0), r.elem(0,0));
  for(int k=1; k < N; ++k)
    d.elem() += traceMultiply(l.elem(0,k), r.elem(k,0));

  for(int j=1; j < N; ++j)
    for(int k=0; k < N; ++k)
      d.elem() += traceMultiply(l.elem(j,k), r.elem(k,j));

  return d;
}

// OScalar = traceMultiply(OMatrix,OScalar)
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, FnTraceMultiply> {
  typedef OScalar<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OScalar<T2>, FnTraceMultiply>::Type_t
traceMultiply(const OMatrix<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OScalar<T2>, FnTraceMultiply>::Type_t  d;

  d.elem() = traceMultiply(l.elem(0,0), r.elem());
  for(int k=1; k < N; ++k)
    d.elem() += traceMultiply(l.elem(k,k), r.elem());

  return d;
}

// OScalar = traceMultiply(OScalar,OMatrix)
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OMatrix<T2>, FnTraceMultiply> {
  typedef OScalar<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OMatrix<T2>, FnTraceMultiply>::Type_t
traceMultiply(const OScalar<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OMatrix<T2>, FnTraceMultiply>::Type_t  d;

  d.elem() = traceMultiply(l.elem(), r.elem(0,0));
  for(int k=1; k < N; ++k)
    d.elem() += traceMultiply(l.elem(), r.elem(k,k));

  return d;
}



//! OMatrix = traceSpinMultiply(OMatrix,OMatrix)   [the trace is an identity in general]
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnTraceSpinMultiply> {
  typedef OMatrix<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const OMatrix<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnTraceSpinMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
    {
      d.elem(i,j) = traceSpinMultiply(l.elem(i,0), r.elem(0,j));
      for(int k=1; k < N; ++k)
	d.elem(i,j) += traceSpinMultiply(l.elem(i,k), r.elem(k,j));
    }

  return d;
}

// OScalar = traceSpinMultiply(OMatrix,OScalar)   [the trace is an identity in general]
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, FnTraceSpinMultiply> {
  typedef OMatrix<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OScalar<T2>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const OMatrix<T1>& l, const OScalar<T2>& r)
{
  typename BinaryReturn<OMatrix<T1>, OScalar<T2>, FnTraceSpinMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = traceSpinMultiply(l.elem(i,j), r.elem());

  return d;
}

// OScalar = traceSpinMultiply(OScalar,OMatrix)   [the trace is an identity in general]
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OMatrix<T2>, FnTraceSpinMultiply> {
  typedef OMatrix<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OMatrix<T2>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const OScalar<T1>& l, const OMatrix<T2>& r)
{
  typename BinaryReturn<OScalar<T1>, OMatrix<T2>, FnTraceSpinMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = traceSpinMultiply(l.elem(), r.elem(i,j));

  return d;
}


//! OMatrix = Re(OMatrix)
template<class T>
struct UnaryReturn<OMatrix<T>, FnReal> {
  typedef OMatrix<typename UnaryReturn<T, FnReal>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnReal>::Type_t
real(const OMatrix<T>& s1)
{
  typename UnaryReturn<OMatrix<T>, FnReal>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = real(s1.elem(i,j));

  return d;
}


//! OMatrix = Im(OMatrix)
template<class T>
struct UnaryReturn<OMatrix<T>, FnImag> {
  typedef OMatrix<typename UnaryReturn<T, FnImag>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnImag>::Type_t
imag(const OMatrix<T>& s1)
{
  typename UnaryReturn<OMatrix<T>, FnImag>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = imag(s1.elem(i,j));

  return d;
}


//! OMatrix<T> = (OMatrix<T> , OMatrix<T>)
template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnCmplx>::Type_t
cmplx(const OMatrix<T1>& s1, const OMatrix<T2>& s2)
{
  typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnCmplx>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = cmplx(s1.elem(i,j), s2.elem(i,j));

  return d;
}




// Functions
//! OMatrix = i * OMatrix
template<class T>
struct UnaryReturn<OMatrix<T>, FnTimesI> {
  typedef OMatrix<typename UnaryReturn<T, FnTimesI>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnTimesI>::Type_t
timesI(const OMatrix<T>& s1)
{
  typename UnaryReturn<OMatrix<T>, FnTimesI>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = timesI(s1.elem(i,j));

  return d;
}

//! OMatrix = -i * OMatrix
template<class T>
struct UnaryReturn<OMatrix<T>, FnTimesMinusI> {
  typedef OMatrix<typename UnaryReturn<T, FnTimesMinusI>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnTimesMinusI>::Type_t
timesMinusI(const OMatrix<T>& s1)
{
  typename UnaryReturn<OMatrix<T>, FnTimesMinusI>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = timesMinusI(s1.elem(i,j));

  return d;
}


//! Extract spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
struct UnaryReturn<OMatrix<T>, FnPeekSpinVector> {
  typedef OMatrix<typename UnaryReturn<T, FnPeekSpinVector>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnPeekSpinVector>::Type_t
peekSpin(const OMatrix<T>& l, int row)
{
  typename UnaryReturn<OMatrix<T>, FnPeekSpinVector>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = peekSpin(l.elem(i,j),row);
  return d;
}

//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
struct UnaryReturn<OMatrix<T>, FnPeekSpinMatrix> {
  typedef OMatrix<typename UnaryReturn<T, FnPeekSpinMatrix>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnPeekSpinMatrix>::Type_t
peekSpin(const OMatrix<T>& l, int row, int col)
{
  typename UnaryReturn<OMatrix<T>, FnPeekSpinMatrix>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = peekSpin(l.elem(i,j),row,col);
  return d;
}

//! Insert spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
struct UnaryReturn<OMatrix<T>, FnPokeSpinMatrix> {
  typedef OMatrix<typename UnaryReturn<T, FnPokeSpinMatrix>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename UnaryReturn<OMatrix<T1>, FnPokeSpinMatrix>::Type_t&
pokeSpin(OMatrix<T1>& l, const OMatrix<T2>& r, int row)
{
  typedef typename UnaryReturn<OMatrix<T1>, FnPokeSpinMatrix>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      pokeSpin(l.elem(i,j),r.elem(i,j),row);
  return static_cast<Return_t&>(l);
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline typename UnaryReturn<OMatrix<T1>, FnPokeSpinMatrix>::Type_t&
pokeSpin(OMatrix<T1>& l, const OMatrix<T2>& r, int row, int col)
{
  typedef typename UnaryReturn<OMatrix<T1>, FnPokeSpinMatrix>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      pokeSpin(l.elem(i,j),r.elem(i,j),row,col);
  return static_cast<Return_t&>(l);
}



//! dest = 0
template<class T> 
inline void 
zero_rep(OMatrix<T>& dest) 
{
  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      zero_rep(dest.elem(i,j));
}


//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline void 
copymask(OMatrix<T>& d, const OScalar<T1>& mask, const OMatrix<T>& s1) 
{
  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      copymask(d.elem(i,j),mask.elem(),s1.elem(i,j));
}


//! dest  = random  
template<class T>
inline void
fill_random(OMatrix<T>& d)
{
  // The skewed_seed is the starting seed to use
  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      fill_random(d.elem(i,j));
}

//! dest  = gaussian
template<class T>
inline void
fill_gaussian(OMatrix<T>& d, OMatrix<T>& r1, OMatrix<T>& r2)
{
  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      fill_gaussian(d.elem(i,j), r1.elem(i,j), r2.elem(i,j));
}



#if 0
// Global sum over site indices only
template<class T>
struct UnaryReturn<OMatrix<T>, FnSum> {
  typedef OMatrix<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnSum>::Type_t
sum(const OMatrix<T>& s1)
{
  typename UnaryReturn<OMatrix<T>, FnSum>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = sum(s1.elem(i,j));

  return d;
}
#endif


// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T>
struct UnaryReturn<OMatrix<T>, FnNorm2> {
  typedef OScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<OMatrix<T>, FnLocalNorm2> {
  typedef OScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OMatrix<T>, FnLocalNorm2>::Type_t
localNorm2(const OMatrix<T>& s1)
{
  typename UnaryReturn<OMatrix<T>, FnLocalNorm2>::Type_t  d;

  d.elem() = localNorm2(s1.elem(0,0));
  for(int j=1; j < N; ++j)
    d.elem() += localNorm2(s1.elem(0,j));

  for(int i=1; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem() += localNorm2(s1.elem(i,j));

  return d;
}


//! OScalar = innerProduct(OMatrix,OMatrix)
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnInnerProduct> {
  typedef OScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

//! OScalar = localInnerProduct(OMatrix,OMatrix)
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnLocalInnerProduct> {
  typedef OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const OMatrix<T1>& s1, const OMatrix<T2>& s2)
{
  typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnLocalInnerProduct>::Type_t  d;

  d.elem() = localInnerProduct(s1.elem(0,0), s2.elem(0,0));
  for(int k=1; k < N; ++k)
    d.elem() += localInnerProduct(s1.elem(k,0), s2.elem(k,0));

  for(int j=1; j < N; ++j)
    for(int k=0; k < N; ++k)
      d.elem() += localInnerProduct(s1.elem(k,j), s2.elem(k,j));

  return d;
}

//! OScalar = localInnerProduct(OMatrix,OScalar)
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, FnLocalInnerProduct> {
  typedef OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OScalar<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const OMatrix<T1>& s1, const OScalar<T2>& s2)
{
  typename BinaryReturn<OMatrix<T1>, OScalar<T2>, FnLocalInnerProduct>::Type_t  d;

  d.elem() = localInnerProduct(s1.elem(0,0), s2.elem());
  for(int k=1; k < N; ++k)
    d.elem() += localInnerProduct(s1.elem(k,k), s2.elem());

  return d;
}

//! OScalar = localInnerProduct(OScalar,OMatrix)
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OMatrix<T2>, FnLocalInnerProduct> {
  typedef OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OMatrix<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const OScalar<T1>& s1, const OMatrix<T2>& s2)
{
  typename BinaryReturn<OScalar<T1>, OMatrix<T2>, FnLocalInnerProduct>::Type_t  d;

  d.elem() = localInnerProduct(s1.elem(), s2.elem(0,0));
  for(int k=1; k < N; ++k)
    d.elem() += localInnerProduct(s1.elem(), s2.elem(k,k));

  return d;
}


//! OScalar = innerProductReal(OMatrix,OMatrix)
/*!
 * return  realpart of InnerProduct(adj(s1)*s2)
 */
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnInnerProductReal > {
  typedef OScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

//! OScalar = innerProductReal(OMatrix,OMatrix)
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnLocalInnerProductReal > {
  typedef OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnLocalInnerProductReal>::Type_t
localInnerProductReal(const OMatrix<T1>& s1, const OMatrix<T2>& s2)
{
  typename BinaryReturn<OMatrix<T1>, OMatrix<T2>, FnLocalInnerProductReal>::Type_t  d;

  d.elem() = localInnerProductReal(s1.elem(0,0), s2.elem(0,0));
  for(int k=1; k < N; ++k)
    d.elem() += localInnerProductReal(s1.elem(k,0), s2.elem(k,0));

  for(int j=1; j < N; ++j)
    for(int k=0; k < N; ++k)
      d.elem() += localInnerProductReal(s1.elem(k,j), s2.elem(k,j));

  return d;
}

//! OScalar = localInnerProductReal(OMatrix,OScalar)
template<class T1, class T2>
struct BinaryReturn<OMatrix<T1>, OScalar<T2>, FnLocalInnerProductReal > {
  typedef OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OMatrix<T1>, OScalar<T2>, FnLocalInnerProductReal>::Type_t
localInnerProductReal(const OMatrix<T1>& s1, const OScalar<T2>& s2)
{
  typename BinaryReturn<OMatrix<T1>, OScalar<T2>, FnLocalInnerProductReal>::Type_t  d;

  d.elem() = localInnerProductReal(s1.elem(0,0), s2.elem());
  for(int k=1; k < N; ++k)
    d.elem() += localInnerProductReal(s1.elem(k,0), s2.elem(k,k));

  return d;
}

//! OScalar = localInnerProductReal(OScalar,OMatrix)
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OMatrix<T2>, FnLocalInnerProductReal > {
  typedef OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OMatrix<T2>, FnLocalInnerProductReal>::Type_t
localInnerProductReal(const OScalar<T1>& s1, const OMatrix<T2>& s2)
{
  typename BinaryReturn<OScalar<T1>, OMatrix<T2>, FnLocalInnerProductReal>::Type_t  d;

  d.elem() = localInnerProductReal(s1.elem(), s2.elem(0,0));
  for(int k=1; k < N; ++k)
    d.elem() += localInnerProductReal(s1.elem(), s2.elem(k,k));

  return d;
}


//! OMatrix<T> = where(OScalar, OMatrix, OMatrix)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<OScalar<T1>, OMatrix<T2>, OMatrix<T3>, FnWhere> {
  typedef OMatrix<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<OScalar<T1>, OMatrix<T2>, OMatrix<T3>, FnWhere>::Type_t
where(const OScalar<T1>& a, const OMatrix<T2>& b, const OMatrix<T3>& c)
{
  typename TrinaryReturn<OScalar<T1>, OMatrix<T2>, OMatrix<T3>, FnWhere>::Type_t  d;

  // Not optimal - want to have where outside assignment
  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = where(a.elem(), b.elem(i,j), c.elem(i,j));

  return d;
}

/*! @} */  // end of group obsmatrix

} // namespace ENSEM

