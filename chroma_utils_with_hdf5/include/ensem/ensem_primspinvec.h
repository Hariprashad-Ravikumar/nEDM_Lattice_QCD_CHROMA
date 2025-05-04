// -*- C++ -*-
// $Id: ensem_primspinvec.h,v 2.0 2008/12/05 04:43:34 edwards Exp $

/*! \file
 * \brief Primitive Spin Vector
 */


namespace ENSEM {

//-------------------------------------------------------------------------------------
/*! \addtogroup primspinvector Spin vector primitive
 * \ingroup primvector
 *
 * Primitive type that transforms like a Spin vector
 *
 * @{
 */

//! Primitive spin Vector class
/*! 
 * Spin vector class supports gamma matrix algebra 
 *
 * NOTE: the class is mostly empty - it is the specialized versions below
 * that know for a fixed size how gamma matrices (constants) should act
 * on the spin vectors.
 */
template <class T, int N> class PSpinVector : public PVector<T, N, PSpinVector>
{
public:
  //! PSpinVector = zero
  inline
  PSpinVector& operator=(const Zero& rhs)
    {
      this->assign(rhs);
      return *this;
    }

  //! PVector = PVector
  /*! Set equal to another PVector */
  template<class T1>
  inline
  PSpinVector& operator=(const PSpinVector<T1,N>& rhs) 
    {
      assign(rhs);
      return *this;
    }

};


//! Specialization of primitive spin Vector class for 4 spin components
/*! 
 * Spin vector class supports gamma matrix algebra for 4 spin components
 */
template<class T> class PSpinVector<T,4> : public PVector<T, 4, PSpinVector>
{
};


//! Specialization of primitive spin Vector class for 2 spin components
/*! 
 * Spin vector class supports gamma matrix algebra for 2 spin components
 * NOTE: this can be used for spin projection tricks of a 4 component spinor
 * to 2 spin components, or a 2 spin component Dirac fermion in 2 dimensions
 */
template<class T> class PSpinVector<T,2> : public PVector<T, 2, PSpinVector>
{
public:
};


/*! @} */   // end of group primspinvec

//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1, int N>
struct WordType<PSpinVector<T1,N> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

// Internally used scalars
template<class T, int N>
struct InternalScalar<PSpinVector<T,N> > {
  typedef PScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive into a scalar leaving grid alone
template<class T, int N>
struct PrimitiveScalar<PSpinVector<T,N> > {
  typedef PScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T, int N>
struct EnsemScalar<PSpinVector<T,N> > {
  typedef PSpinVector<typename EnsemScalar<T>::Type_t, N>  Type_t;
};

// Traits class to label IO types
template<class T, int N>
struct EnsbcIO<PSpinVector<T,N> > {
  enum {type = EnsbcIO<T>::type};
};

//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(PSpinVector) -> PSpinVector
template<class T1, int N, class Op>
struct UnaryReturn<PSpinVector<T1,N>, Op> {
  typedef PSpinVector<typename UnaryReturn<T1, Op>::Type_t, N>  Type_t;
};
// Default binary(PScalar,PSpinVector) -> PSpinVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinMatrix,PSpinVector) -> PSpinVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinVector<T2,N>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinVector,PScalar) -> PSpinVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinVector,PSpinVector) -> PSpinVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};


#if 0
template<class T1, class T2>
struct UnaryReturn<PScalar<T2>, OpCast<T1> > {
  typedef PScalar<typename UnaryReturn<T, OpCast>::Type_t>  Type_t;
//  typedef T1 Type_t;
};
#endif


// Assignment is different
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpAddAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpSubtractAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpMultiplyAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpDivideAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};
 


// SpinVector
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnLocalNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnLocalInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnLocalInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};



// Gamma algebra
template<int m, class T2, int N>
struct BinaryReturn<GammaConst<N,m>, PSpinVector<T2,N>, OpGammaConstMultiply> {
  typedef PSpinVector<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N>
struct BinaryReturn<GammaType<N>, PSpinVector<T2,N>, OpGammaTypeMultiply> {
  typedef PSpinVector<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};


// Gamma algebra
template<int m, class T2, int N>
struct BinaryReturn<GammaConstDP<N,m>, PSpinVector<T2,N>, OpGammaConstDPMultiply> {
  typedef PSpinVector<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N>
struct BinaryReturn<GammaTypeDP<N>, PSpinVector<T2,N>, OpGammaTypeDPMultiply> {
  typedef PSpinVector<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};


//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

/*! \addtogroup primspinvector */
/*! @{ */

// Peeking and poking
//! Extract spin vector components 
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnPeekSpinVector > {
  typedef PScalar<typename UnaryReturn<T, FnPeekSpinVector>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnPeekSpinVector>::Type_t
peekSpin(const PSpinVector<T,N>& l, int row)
{
  typename UnaryReturn<PSpinVector<T,N>, FnPeekSpinVector>::Type_t  d;

  // Note, do not need to propagate down since the function is eaten at this level
  d.elem() = l.elem(row);
  return d;
}

//! Insert spin vector components
template<class T1, class T2, int N>
inline PSpinVector<T1,N>&
pokeSpin(PSpinVector<T1,N>& l, const PScalar<T2>& r, int row)
{
  // Note, do not need to propagate down since the function is eaten at this level
  l.elem(row) = r.elem();
  return l;
}


//-----------------------------------------------

// SpinVector<4> = Gamma<4,m> * SpinVector<4>
// There are 16 cases here for Nd=4
template<class T2>
inline typename BinaryReturn<GammaConst<4,0>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,0>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,0>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;
  
  d.elem(0) =  r.elem(0);
  d.elem(1) =  r.elem(1);
  d.elem(2) =  r.elem(2);
  d.elem(3) =  r.elem(3);

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,1>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,1>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,1>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;
  
  d.elem(0) = timesI(r.elem(3));
  d.elem(1) = timesI(r.elem(2));
  d.elem(2) = timesMinusI(r.elem(1));
  d.elem(3) = timesMinusI(r.elem(0));

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,2>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,2>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,2>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = -r.elem(3);
  d.elem(1) =  r.elem(2);
  d.elem(2) =  r.elem(1);
  d.elem(3) = -r.elem(0);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,3>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,3>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,3>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(0));
  d.elem(1) = timesI(r.elem(1));
  d.elem(2) = timesMinusI(r.elem(2));
  d.elem(3) = timesI(r.elem(3));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,4>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,4>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,4>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(2));
  d.elem(1) = timesMinusI(r.elem(3));
  d.elem(2) = timesMinusI(r.elem(0));
  d.elem(3) = timesI(r.elem(1));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,5>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,5>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,5>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = -r.elem(1);
  d.elem(1) =  r.elem(0);
  d.elem(2) = -r.elem(3);
  d.elem(3) =  r.elem(2);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,6>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,6>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,6>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(1));
  d.elem(1) = timesMinusI(r.elem(0));
  d.elem(2) = timesMinusI(r.elem(3));
  d.elem(3) = timesMinusI(r.elem(2));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,7>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,7>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,7>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) =  r.elem(2);
  d.elem(1) =  r.elem(3);
  d.elem(2) = -r.elem(0);
  d.elem(3) = -r.elem(1);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,8>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,8>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,8>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) =  r.elem(2);
  d.elem(1) =  r.elem(3);
  d.elem(2) =  r.elem(0);
  d.elem(3) =  r.elem(1);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,9>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,9>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,9>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(1));
  d.elem(1) = timesI(r.elem(0));
  d.elem(2) = timesMinusI(r.elem(3));
  d.elem(3) = timesMinusI(r.elem(2));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,10>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,10>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,10>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = -r.elem(1);
  d.elem(1) =  r.elem(0);
  d.elem(2) =  r.elem(3);
  d.elem(3) = -r.elem(2);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,11>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,11>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,11>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(2));
  d.elem(1) = timesI(r.elem(3));
  d.elem(2) = timesMinusI(r.elem(0));
  d.elem(3) = timesI(r.elem(1));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,12>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,12>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,12>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(0));
  d.elem(1) = timesMinusI(r.elem(1));
  d.elem(2) = timesMinusI(r.elem(2));
  d.elem(3) = timesI(r.elem(3));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,13>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,13>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,13>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = -r.elem(3);
  d.elem(1) =  r.elem(2);
  d.elem(2) = -r.elem(1);
  d.elem(3) =  r.elem(0);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,14>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,14>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,14>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(3));
  d.elem(1) = timesMinusI(r.elem(2));
  d.elem(2) = timesMinusI(r.elem(1));
  d.elem(3) = timesMinusI(r.elem(0));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,15>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,15>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,15>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) =  r.elem(0);
  d.elem(1) =  r.elem(1);
  d.elem(2) = -r.elem(2);
  d.elem(3) = -r.elem(3);
  
  return d;
}


//-----------------------------------------------

// SpinVector<4> = GammaDP<4,m> * SpinVector<4>
// There are 16 cases here for Nd=4
template<class T2>
inline typename BinaryReturn<GammaConstDP<4,0>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,0>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,0>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;
  
  d.elem(0) =  r.elem(0);
  d.elem(1) =  r.elem(1);
  d.elem(2) =  r.elem(2);
  d.elem(3) =  r.elem(3);

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,1>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,1>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,1>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;
  
  d.elem(0) = timesMinusI(r.elem(3));
  d.elem(1) = timesMinusI(r.elem(2));
  d.elem(2) = timesI(r.elem(1));
  d.elem(3) = timesI(r.elem(0));

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,2>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,2>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,2>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = -r.elem(3);
  d.elem(1) =  r.elem(2);
  d.elem(2) =  r.elem(1);
  d.elem(3) = -r.elem(0);

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,3>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,3>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,3>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(0));
  d.elem(1) = timesMinusI(r.elem(1));
  d.elem(2) = timesI(r.elem(2));
  d.elem(3) = timesMinusI(r.elem(3));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,4>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,4>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,4>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(2));
  d.elem(1) = timesI(r.elem(3));
  d.elem(2) = timesI(r.elem(0));
  d.elem(3) = timesMinusI(r.elem(1));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,5>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,5>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,5>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = -r.elem(1);
  d.elem(1) =  r.elem(0);
  d.elem(2) = -r.elem(3);
  d.elem(3) =  r.elem(2);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,6>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,6>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,6>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(1));
  d.elem(1) = timesI(r.elem(0));
  d.elem(2) = timesI(r.elem(3));
  d.elem(3) = timesI(r.elem(2));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,7>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,7>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,7>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) =  r.elem(2);
  d.elem(1) =  r.elem(3);
  d.elem(2) = -r.elem(0);
  d.elem(3) = -r.elem(1);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,8>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,8>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,8>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) =  r.elem(0);
  d.elem(1) =  r.elem(1);
  d.elem(2) = -r.elem(2);
  d.elem(3) = -r.elem(3);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,9>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,9>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,9>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(3));
  d.elem(1) = timesI(r.elem(2));
  d.elem(2) = timesI(r.elem(1));
  d.elem(3) = timesI(r.elem(0));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,10>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,10>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,10>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) =  r.elem(3);
  d.elem(1) = -r.elem(2);
  d.elem(2) = -r.elem(1);
  d.elem(3) =  r.elem(0);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,11>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,11>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,11>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(0));
  d.elem(1) = timesMinusI(r.elem(1));
  d.elem(2) = timesMinusI(r.elem(2));
  d.elem(3) = timesI(r.elem(3));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,12>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,12>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,12>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(2));
  d.elem(1) = timesMinusI(r.elem(3));
  d.elem(2) = timesI(r.elem(0));
  d.elem(3) = timesMinusI(r.elem(1));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,13>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,13>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,13>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = -r.elem(1);
  d.elem(1) =  r.elem(0);
  d.elem(2) =  r.elem(3);
  d.elem(3) = -r.elem(2);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,14>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,14>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,14>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(1));
  d.elem(1) = timesI(r.elem(0));
  d.elem(2) = timesMinusI(r.elem(3));
  d.elem(3) = timesMinusI(r.elem(2));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,15>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,15>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,15>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = -r.elem(2);
  d.elem(1) = -r.elem(3);
  d.elem(2) = -r.elem(0);
  d.elem(3) = -r.elem(1);
  
  return d;
}


//-----------------------------------------------------------------------------
//! PSpinVector<T,4> = P_+ * PSpinVector<T,4>
template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnChiralProjectPlus>::Type_t
chiralProjectPlus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnChiralProjectPlus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  zero_rep(d.elem(2));
  zero_rep(d.elem(3));

  return d;
}

//! PSpinVector<T,4> = P_- * PSpinVector<T,4>
template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnChiralProjectMinus>::Type_t
chiralProjectMinus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnChiralProjectMinus>::Type_t  d;

  zero_rep(d.elem(0));
  zero_rep(d.elem(1));
  d.elem(2) = s1.elem(2);
  d.elem(3) = s1.elem(3);

  return d;
}


/*! @} */   // end of group primspinvector

} // namespace ENSEM

