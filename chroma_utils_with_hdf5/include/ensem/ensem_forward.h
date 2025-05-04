// -*- C++ -*-
// $Id: ensem_forward.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! @file
 * @brief Forward declarations for Ensem
 */

namespace ENSEM {

// Reality
template<class T> class RScalar;
template<class T> class RComplex;

// Primitives
template<class T> class PScalar;
template <class T, int N, template<class,int> class C> class PMatrix;
template <class T, int N, template<class,int> class C> class PVector;
template <class T, int N> class PColorVector;
template <class T, int N> class PSpinVector;
template <class T, int N> class PColorMatrix;
template <class T, int N> class PSpinMatrix;
template <class T> class PSeed;

template<int N> class GammaType;
template<int N, int m> class GammaConst;

template<int N> class GammaTypeDP;
template<int N, int m> class GammaConstDP;

// Observables
template<class T> class OScalar;
template<class T> class OVector;
template<class T> class OMatrix;

// Ensemble
template<class T> class EScalar;
template<class T> class Ensem;

// Simple scalar trait class
template<class T> struct SimpleScalar;
template<class T> struct InternalScalar;
template<class T> struct EnsemScalar;
template<class T> struct PrimitiveScalar;
template<class T> struct RealScalar;
template<class T> struct WordType;

} // namespace ENSEM
  
