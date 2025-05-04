// -*- C++ -*-
// $Id: ensem_obsscalar.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! \file
 * \brief Observable Scalar
 */

namespace ENSEM
{

//-------------------------------------------------------------------------------------
/*! \addtogroup obsscalar Scalar observable
 * \ingroup fiber
 *
 * Observable Scalar is a placeholder for no observable structure
 *
 * @{
 */

//! Observable Scalar
/*! Placeholder for no observable structure */
template<class T> class OScalar
{
public:
  OScalar() {}
  ~OScalar() {}

  //---------------------------------------------------------
  //! construct dest = const
  OScalar(const typename WordType<T>::Type_t& rhs) : F(rhs) {}

  //! construct dest = rhs
  template<class T1>
  OScalar(const OScalar<T1>& rhs) : F(rhs.elem()) {}

  //! construct dest = rhs
  template<class T1>
  OScalar(const T1& rhs) : F(rhs) {}

  //! construct dest = 0
  OScalar(const Zero& rhs) {zero_rep(F);}

  //---------------------------------------------------------
#if 1
  //! dest = const
  /*! Fill with an integer constant. Will be promoted to underlying word type */
  inline
  OScalar& operator=(const typename WordType<T>::Type_t& rhs)
    {
      elem() = rhs;
      return *this;
    }
#endif

  //! OScalar = zero
  inline
  OScalar& operator=(const Zero& rhs)
    {
      zero_rep(elem());
      return *this;
    }

  //! OScalar = OScalar
  /*! Set equal to another OScalar */
  inline
  OScalar& operator=(const OScalar& rhs) 
    {
      elem() = rhs.elem();
      return *this;
    }

  //! OScalar = OScalar
  /*! Set equal to another OScalar */
  template<class T1>
  inline
  OScalar& operator=(const OScalar<T1>& rhs) 
    {
      elem() = rhs.elem();
      return *this;
    }

  //! OScalar += OScalar
  template<class T1>
  inline
  OScalar& operator+=(const OScalar<T1>& rhs) 
    {
      elem() += rhs.elem();
      return *this;
    }

  //! OScalar -= OScalar
  template<class T1>
  inline
  OScalar& operator-=(const OScalar<T1>& rhs) 
    {
      elem() -= rhs.elem();
      return *this;
    }

  //! OScalar *= OScalar
  template<class T1>
  inline
  OScalar& operator*=(const OScalar<T1>& rhs) 
    {
      elem() *= rhs.elem();
      return *this;
    }

  //! OScalar /= OScalar
  template<class T1>
  inline
  OScalar& operator/=(const OScalar<T1>& rhs) 
    {
      elem() /= rhs.elem();
      return *this;
    }

  //! OScalar %= OScalar
  template<class T1>
  inline
  OScalar& operator%=(const OScalar<T1>& rhs) 
    {
      elem() %= rhs.elem();
      return *this;
    }

  //! OScalar |= OScalar
  template<class T1>
  inline
  OScalar& operator|=(const OScalar<T1>& rhs) 
    {
      elem() |= rhs.elem();
      return *this;
    }

  //! OScalar &= OScalar
  template<class T1>
  inline
  OScalar& operator&=(const OScalar<T1>& rhs) 
    {
      elem() &= rhs.elem();
      return *this;
    }

  //! OScalar ^= OScalar
  template<class T1>
  inline
  OScalar& operator^=(const OScalar<T1>& rhs) 
    {
      elem() ^= rhs.elem();
      return *this;
    }

  //! OScalar <<= OScalar
  template<class T1>
  inline
  OScalar& operator<<=(const OScalar<T1>& rhs) 
    {
      elem() <<= rhs.elem();
      return *this;
    }

  //! OScalar >>= OScalar
  template<class T1>
  inline
  OScalar& operator>>=(const OScalar<T1>& rhs) 
    {
      elem() >>= rhs.elem();
      return *this;
    }

  //! Deep copies here
  OScalar(const OScalar& a): F(a.F) {/* fprintf(stderr,"copy OScalar\n"); */}

  //---------------------------------------------------------
public:
  //! The backdoor
  inline const T* getF() const {return &F;}
  inline T* getF() {return &F;}

public:
  //! Has this object been initialized (resized and such)
  bool initP() const {return true;}

  inline int size() const {return 1;}
  inline int numElem() const {return 1;}

  inline void resize(const OScalar&) {}
  inline void resize(int n1) {}

public:
  inline T& elem() {return F;}
  inline const T& elem() const {return F;}

private:
  T F;
};


// I/O
//! Binary input
template<class T>  
inline
void read(ADATIO::BinaryReader& bin, OScalar<T>& d)
{
  read(bin, d.elem());
}

//! Binary output
template<class T>
inline
void write(ADATIO::BinaryWriter& bin, const OScalar<T>& d)
{
  write(bin, d.elem());
}

//! Ascii input
template<class T>
inline
std::istream& operator>>(std::istream& s, OScalar<T>& d)
{
  return s >> d.elem();
}

// Output
//! Ascii output
template<class T>
inline
std::ostream& operator<<(std::ostream& s, const OScalar<T>& d)
{
  return s << d.elem();
}

//! Text input
template<class T>
inline
ENSEM::TextReader& operator>>(ENSEM::TextReader& txt, OScalar<T>& d)
{
  int k;
  txt >> k >> d.elem();
  if (k != 0)
  {
    std::cerr << "error reading OScalar" << std::endl;
    exit(1);
  }
  return txt;
}

//! Text output
template<class T>
inline
ENSEM::TextWriter& operator<<(ENSEM::TextWriter& txt, const OScalar<T>& d)
{
  return txt << " 0 " << d.elem();
}

//! XML output
template<class T>
inline
ADATXML::XMLWriter& operator<<(ADATXML::XMLWriter& xml, const OScalar<T>& d)
{
  return xml << d.elem();
}

//! XML input
template<class T>
inline
void read(ADATXML::XMLReader& xml, const std::string& path, OScalar<T>& d)
{
  read(xml, path, d.elem());
}


/*! @} */  // end of group obsscalar


//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T>
struct WordType<OScalar<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};

// Internally used scalars
template<class T>
struct InternalScalar<OScalar<T> > {
  typedef OScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Internally used real scalars
template<class T>
struct RealScalar<OScalar<T> > {
  typedef OScalar<typename RealScalar<T>::Type_t>  Type_t;
};

// Makes a observable scalar leaving grid alone
template<class T>
struct PrimitiveScalar<OScalar<T> > {
  typedef OScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving observable indices alone
template<class T>
struct EnsemScalar<OScalar<T> > {
  typedef OScalar<typename EnsemScalar<T>::Type_t>  Type_t;
};

// Traits class to label IO types
template<class T> 
struct EnsbcIO<OScalar<T> > {
  enum {type = EnsbcIO<T>::type};
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(OScalar) -> OScalar
template<class T1, class Op>
struct UnaryReturn<OScalar<T1>, Op> {
  typedef OScalar<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};

// Default binary(OScalar,OScalar) -> OScalar
template<class T1, class T2, class Op>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, Op> {
  typedef OScalar<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};


// Assignment is different
template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpAssign > {
  typedef OScalar<T1> &Type_t;
};

template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpAddAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpSubtractAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpMultiplyAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpDivideAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpModAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpBitwiseOrAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpBitwiseAndAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpBitwiseXorAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpLeftShiftAssign > {
  typedef OScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpRightShiftAssign > {
  typedef OScalar<T1> &Type_t;
};
 



//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

/*! \addtogroup obsscalar */
/*! @{ */

// Observable Scalars

// ! OScalar
template<class T>
struct UnaryReturn<OScalar<T>, OpNot > {
  typedef OScalar<typename UnaryReturn<T, OpNot>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OScalar<T1>, OpNot>::Type_t
operator!(const OScalar<T1>& l)
{
  return ! l.elem();
}

// + OScalar
template<class T1>
inline typename UnaryReturn<OScalar<T1>, OpUnaryPlus>::Type_t
operator+(const OScalar<T1>& l)
{
  return +l.elem();
}

// - OScalar
template<class T1>
inline typename UnaryReturn<OScalar<T1>, OpUnaryMinus>::Type_t
operator-(const OScalar<T1>& l)
{
  return -l.elem();
}

// OScalar + OScalar
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpAdd>::Type_t
operator+(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() + r.elem();
}

// OScalar - OScalar
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpSubtract>::Type_t
operator-(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() - r.elem();
}

// OScalar * OScalar
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpMultiply>::Type_t
operator*(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() * r.elem();
}

// Optimized  adj(PMatrix)*PMatrix
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpAdjMultiply>::Type_t
adjMultiply(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return adjMultiply(l.elem(), r.elem());
}

// Optimized  PMatrix*adj(PMatrix)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return multiplyAdj(l.elem(), r.elem());
}

// Optimized  adj(PMatrix)*adj(PMatrix)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return adjMultiplyAdj(l.elem(), r.elem());
}

// OScalar / OScalar
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpDivide>::Type_t
operator/(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() / r.elem();
}


// OScalar << OScalar
template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpLeftShift > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpLeftShift>::Type_t
operator<<(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() << r.elem();
}

// OScalar >> OScalar
template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpRightShift > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpRightShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpRightShift>::Type_t
operator>>(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() >> r.elem();
}

// OScalar % OScalar
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpMod>::Type_t
operator%(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() % r.elem();
}

// OScalar ^ OScalar
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpBitwiseXor>::Type_t
operator^(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() ^ r.elem();
}

// OScalar & OScalar
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpBitwiseAnd>::Type_t
operator&(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() & r.elem();
}

// OScalar | OScalar
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpBitwiseOr>::Type_t
operator|(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() | r.elem();
}


// Comparisons
template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpLT > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpLT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpLT>::Type_t
operator<(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() < r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpLE > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpLE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpLE>::Type_t
operator<=(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() <= r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpGT > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpGT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpGT>::Type_t
operator>(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() > r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpGE > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpGE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpGE>::Type_t
operator>=(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() >= r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpEQ > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpEQ>::Type_t
operator==(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() == r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpNE > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpNE>::Type_t
operator!=(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() != r.elem();
}


template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpAnd > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpAnd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpAnd>::Type_t
operator&&(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() && r.elem();
}


template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, OpOr > {
  typedef OScalar<typename BinaryReturn<T1, T2, OpOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, OpOr>::Type_t
operator||(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return l.elem() || r.elem();
}


//-----------------------------------------------------------------------------
// Functions

// Adjoint
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnAdjoint>::Type_t
adj(const OScalar<T1>& s1)
{
  return adj(s1.elem());
}


// Conjugate
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnConjugate>::Type_t
conj(const OScalar<T1>& s1)
{
  return conj(s1.elem());
}


// Transpose
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnTranspose>::Type_t
transpose(const OScalar<T1>& s1)
{
  return transpose(s1.elem());
}


// TRACE
// trace = Trace(source1)
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnTrace>::Type_t
trace(const OScalar<T1>& s1)
{
  return trace(s1.elem());
}


// trace = Re(Trace(source1))
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnRealTrace>::Type_t
realTrace(const OScalar<T1>& s1)
{
  return realTrace(s1.elem());
}


// trace = Im(Trace(source1))
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnImagTrace>::Type_t
imagTrace(const OScalar<T1>& s1)
{
  return imagTrace(s1.elem());
}


// trace = colorTrace(source1)
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnTraceColor>::Type_t
traceColor(const OScalar<T1>& s1)
{
  return traceColor(s1.elem());
}


//! OScalar = traceSpin(OScalar)
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnTraceSpin>::Type_t
traceSpin(const OScalar<T1>& s1)
{
  return traceSpin(s1.elem());
}

//! OScalar = trace(OScalar * OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnTraceMultiply>::Type_t
traceMultiply(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return traceMultiply(l.elem(), r.elem());
}

//! OScalar = traceColor(OScalar * OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnTraceColorMultiply>::Type_t
traceColorMultiply(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return traceMultiply(l.elem(), r.elem());
}

//! OScalar = traceSpin(OScalar * OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return traceMultiply(l.elem(), r.elem());
}

//! OScalar = outerProduct(OScalar, OScalar)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const OScalar<T1>& l, const OScalar<T2>& r)
{
  return outerProduct(l.elem(),r.elem());
}


//! OScalar = Re(OScalar)
template<class T>
inline typename UnaryReturn<OScalar<T>, FnReal>::Type_t
real(const OScalar<T>& s1)
{
  return real(s1.elem());
}


// OScalar = Im(OScalar)
template<class T>
inline typename UnaryReturn<OScalar<T>, FnImag>::Type_t
imag(const OScalar<T>& s1)
{
  return imag(s1.elem());
}


// ArcCos
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnArcCos>::Type_t
acos(const OScalar<T1>& s1)
{
  return acos(s1.elem());
}

// ArcSin
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnArcSin>::Type_t
asin(const OScalar<T1>& s1)
{
  return asin(s1.elem());
}

// ArcTan
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnArcTan>::Type_t
atan(const OScalar<T1>& s1)
{
  return atan(s1.elem());
}

// Cos
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnCos>::Type_t
cos(const OScalar<T1>& s1)
{
  return cos(s1.elem());
}

// Exp
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnExp>::Type_t
exp(const OScalar<T1>& s1)
{
  return exp(s1.elem());
}

// Fabs
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnFabs>::Type_t
fabs(const OScalar<T1>& s1)
{
  return fabs(s1.elem());
}

// Log
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnLog>::Type_t
log(const OScalar<T1>& s1)
{
  return log(s1.elem());
}

// Sin
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnSin>::Type_t
sin(const OScalar<T1>& s1)
{
  return sin(s1.elem());
}

// Sqrt
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnSqrt>::Type_t
sqrt(const OScalar<T1>& s1)
{
  return sqrt(s1.elem());
}

// Tan
template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnTan>::Type_t
tan(const OScalar<T1>& s1)
{
  return tan(s1.elem());
}



//! OScalar<T> = pow(OScalar<T> , OScalar<T>)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnPow>::Type_t
pow(const OScalar<T1>& s1, const OScalar<T2>& s2)
{
  return pow(s1.elem(), s2.elem());
}

//! OScalar<T> = atan2(OScalar<T> , OScalar<T>)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnArcTan2>::Type_t
atan2(const OScalar<T1>& s1, const OScalar<T2>& s2)
{
  return atan2(s1.elem(), s2.elem());
}


//! OScalar<T> = (OScalar<T> , OScalar<T>)
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnCmplx>::Type_t
cmplx(const OScalar<T1>& s1, const OScalar<T2>& s2)
{
  return cmplx(s1.elem(), s2.elem());
}



// Global Functions
// OScalar = i * OScalar
template<class T>
inline typename UnaryReturn<OScalar<T>, FnTimesI>::Type_t
timesI(const OScalar<T>& s1)
{
  return timesI(s1.elem());
}

// OScalar = -i * OScalar
template<class T>
inline typename UnaryReturn<OScalar<T>, FnTimesMinusI>::Type_t
timesMinusI(const OScalar<T>& s1)
{
  return timesMinusI(s1.elem());
}


//! bool = isZero(OScalar)
template<class T>
bool
isZero(const OScalar<T>& s1)
{
  return isZero(s1.elem());
}


//! bool = isNaN(OScalar)
template<class T>
bool
isNaN(const OScalar<T>& s1)
{
  return isNaN(s1.elem());
}


//! bool = isInf(OScalar)
template<class T>
bool
isInf(const OScalar<T>& s1)
{
  return isInf(s1.elem());
}


//! bool = isFinite(OScalar)
template<class T>
bool
isFinite(const OScalar<T>& s1)
{
  return isFinite(s1.elem());
}


//! dest [float type] = source [seed type]
template<class T>
inline typename UnaryReturn<OScalar<T>, FnSeedToFloat>::Type_t
seedToFloat(const OScalar<T>& s1)
{
  return seedToFloat(s1.elem());
}


// OScalar = shift(OScalar, int)
template<class T1>
struct UnaryReturn<OScalar<T1>, FnShift> {
  typedef OScalar<typename UnaryReturn<T1, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnShift>::Type_t
shift(const OScalar<T1>& s1, int sh)
{
  return s1.elem();
}


// OScalar = cshift(OScalar, int)
template<class T1>
struct UnaryReturn<OScalar<T1>, FnCshift> {
  typedef OScalar<typename UnaryReturn<T1, FnCshift>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<OScalar<T1>, FnCshift>::Type_t
cshift(const OScalar<T1>& s1, int sh)
{
  return s1.elem();
}


//! Extract observable vector components 
/*! Generically, this is an identity operation. Defined differently under observable */
template<class T1>
struct UnaryReturn<OScalar<T1>, FnPeekObsVector> {
  typedef OScalar<typename UnaryReturn<T1, FnPeekObsVector>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OScalar<T>, FnPeekObsVector>::Type_t
peekObs(const OScalar<T>& l, int row)
{
  return l.elem();
}

//! Extract color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T>
inline typename UnaryReturn<OScalar<T>, FnPeekColorVector>::Type_t
peekColor(const OScalar<T>& l, int row)
{
  return peekColor(l.elem(),row);
}

//! Extract color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T>
inline typename UnaryReturn<OScalar<T>, FnPeekColorMatrix>::Type_t
peekColor(const OScalar<T>& l, int row, int col)
{
  return peekColor(l.elem(),row,col);
}

//! Extract spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1>
struct UnaryReturn<OScalar<T1>, FnPeekSpinVector> {
  typedef OScalar<typename UnaryReturn<T1, FnPeekSpinVector>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OScalar<T>, FnPeekSpinVector>::Type_t
peekSpin(const OScalar<T>& l, int row)
{
  return peekSpin(l.elem(),row);
}

//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<OScalar<T>, FnPeekSpinMatrix>::Type_t
peekSpin(const OScalar<T>& l, int row, int col)
{
  return peekSpin(l.elem(),row,col);
}


//! Insert color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T1, class T2>
inline OScalar<T1>&
pokeColor(OScalar<T1>& l, const OScalar<T2>& r, int row)
{
  pokeColor(l.elem(),r.elem(),row);
  return l;
}

//! Insert color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T1, class T2>
inline OScalar<T1>&
pokeColor(OScalar<T1>& l, const OScalar<T2>& r, int row, int col)
{
  pokeColor(l.elem(),r.elem(),row,col);
  return l;
}

//! Insert spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline OScalar<T1>&
pokeSpin(OScalar<T1>& l, const OScalar<T2>& r, int row)
{
  pokeSpin(l.elem(),r.elem(),row);
  return l;
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline OScalar<T1>&
pokeSpin(OScalar<T1>& l, const OScalar<T2>& r, int row, int col)
{
  pokeSpin(l.elem(),r.elem(),row,col);
  return l;
}


//-----------------------------------------------------------------------------
// Gamma algebra
template<int N, int m, class T, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, OScalar<T>, OpGammaConstMultiply> {
  typedef OScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! OScalar = Gamma<N,m> * OScalar
template<class T, int N, int m>
inline typename BinaryReturn<GammaConst<N,m>, OScalar<T>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<N,m>& l, const OScalar<T>& r)
{
  return l * r.elem();
}

template<class T, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<OScalar<T>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef OScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! OScalar = OScalar * Gamma<N,m>
template<class T, int N, int m>
inline typename BinaryReturn<OScalar<T>, GammaConst<N,m>, OpMultiplyGammaConst>::Type_t
operator*(const OScalar<T>& l, const GammaConst<N,m>& r)
{
  return l.elem() * r;
}


//-----------------------------------------------------------------------------
// Gamma algebra
template<int N, int m, class T, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, OScalar<T>, OpGammaConstDPMultiply> {
  typedef OScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! OScalar = Gamma<N,m> * OScalar
template<class T, int N, int m>
inline typename BinaryReturn<GammaConstDP<N,m>, OScalar<T>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<N,m>& l, const OScalar<T>& r)
{
  return l * r.elem();
}

template<class T, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<OScalar<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef OScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! OScalar = OScalar * Gamma<N,m>
template<class T, int N, int m>
inline typename BinaryReturn<OScalar<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP>::Type_t
operator*(const OScalar<T>& l, const GammaConstDP<N,m>& r)
{
  return l.elem() * r;
}

//-----------------------------------------------------------------------------
//! OScalar = chiralProjectPlus(OScalar)
template<class T>
inline typename UnaryReturn<OScalar<T>, FnChiralProjectPlus>::Type_t
chiralProjectPlus(const OScalar<T>& s1)
{
  return chiralProjectPlus(s1.elem());
}

//! OScalar = chiralProjectMinus(OScalar)
template<class T>
inline typename UnaryReturn<OScalar<T>, FnChiralProjectMinus>::Type_t
chiralProjectMinus(const OScalar<T>& s1)
{
  return chiralProjectMinus(s1.elem());
}


//-----------------------------------------------------------------------------
//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline void 
copymask(OScalar<T>& d, const OScalar<T1>& mask, const OScalar<T>& s1) 
{
  copymask(d.elem(),mask.elem(),s1.elem());
}

//! dest  = random  
template<class T>
inline void
fill_random(OScalar<T>& d)
{
  fill_random(d.elem());
}


//! dest  = gaussian  
template<class T>
inline void
fill_gaussian(OScalar<T>& d, OScalar<T>& r1, OScalar<T>& r2)
{
  fill_gaussian(d.elem(), r1.elem(), r2.elem());
}


#if 1
// Global sum over site indices only
template<class T>
struct UnaryReturn<OScalar<T>, FnSum > {
  typedef OScalar<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OScalar<T>, FnSum>::Type_t
sum(const OScalar<T>& s1)
{
  return sum(s1.elem());
}
#endif


// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T>
struct UnaryReturn<OScalar<T>, FnNorm2 > {
  typedef OScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<OScalar<T>, FnLocalNorm2 > {
  typedef OScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<OScalar<T>, FnLocalNorm2>::Type_t
localNorm2(const OScalar<T>& s1)
{
  return localNorm2(s1.elem());
}



//! OScalar<T> = InnerProduct(adj(OScalar<T1>)*OScalar<T2>)
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProduct > {
  typedef OScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, FnLocalInnerProduct > {
  typedef OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const OScalar<T1>& s1, const OScalar<T2>& s2)
{
  return localInnerProduct(s1.elem(), s2.elem());
}


//! OScalar<T> = InnerProductReal(adj(PMatrix<T1>)*PMatrix<T1>)
template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProductReal > {
  typedef OScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<OScalar<T1>, OScalar<T2>, FnLocalInnerProductReal > {
  typedef OScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnLocalInnerProductReal>::Type_t
localInnerProductReal(const OScalar<T1>& s1, const OScalar<T2>& s2)
{
  return localInnerProductReal(s1.elem(), s2.elem());
}


//! OScalar<T> = where(OScalar, OScalar, OScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<OScalar<T1>, OScalar<T2>, OScalar<T3>, FnWhere> {
  typedef OScalar<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<OScalar<T1>, OScalar<T2>, OScalar<T3>, FnWhere>::Type_t
where(const OScalar<T1>& a, const OScalar<T2>& b, const OScalar<T3>& c)
{
  return where(a.elem(), b.elem(), c.elem());
}


//-----------------------------------------------------------------------------
//! ENSEM Int to int observable in conversion routine
template<class T> 
inline int 
toInt(const OScalar<T>& s) 
{
  return toInt(s.elem());
}

//! ENSEM Real to float observable in conversion routine
template<class T> 
inline float
toFloat(const OScalar<T>& s) 
{
  return toFloat(s.elem());
}

//! ENSEM Double to double observable in conversion routine
template<class T> 
inline double
toDouble(const OScalar<T>& s) 
{
  return toDouble(s.elem());
}

//! ENSEM Boolean to bool observable in conversion routine
template<class T> 
inline bool
toBool(const OScalar<T>& s) 
{
  return toBool(s.elem());
}


//-----------------------------------------------------------------------------
// Other operations
//! dest = 0
template<class T> 
inline void 
zero_rep(OScalar<T>& dest) 
{
  zero_rep(dest.elem());
}

//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
cast_rep(T& d, const OScalar<T1>& s1)
{
  cast_rep(d, s1.elem());
}

//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
cast_rep(OScalar<T>& d, const OScalar<T1>& s1)
{
  cast_rep(d.elem(), s1.elem());
}

/*! @} */  // end of group obsscalar

}  // namespace ENSEM

