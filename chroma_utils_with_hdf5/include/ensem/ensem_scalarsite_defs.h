// -*- C++ -*-
// $Id: ensem_scalarsite_defs.h,v 2.0 2008/12/05 04:43:34 edwards Exp $

/*! \file
 * \brief Type definitions
 */

namespace ENSEM {

/*! \addtogroup defs Type definitions
 *
 * User constructed types made from ENSEM type compositional nesting.
 * The layout is suitable for a scalar-like implementation. Namely,
 * configurations are the slowest varying index.
 *
 * @{
 */

#include "ensem_precision.h"

//----------------------------------------------------------------------
//! Gamma matrices are conveniently defined for this Ns
typedef GammaType<Ns> Gamma;

//! Gamma matrices are conveniently defined for this Ns
typedef GammaTypeDP<Ns> GammaDP;


// Aliases for a scalar architecture

// Floating aliases
typedef Ensem< OScalar< PSpinVector< PScalar< RComplex<REAL> >, Ns> > > EnsemSpinVector;
typedef Ensem< OScalar< PSpinMatrix< PScalar< RComplex<REAL> >, Ns> > > EnsemSpinMatrix;
typedef Ensem< OScalar< PScalar< PScalar< RComplex<REAL> > > > > EnsemComplex;

typedef Ensem< OVector< PScalar< PScalar< RComplex<REAL> > > > > EnsemVectorComplex;
typedef Ensem< OVector< PScalar< PScalar< RScalar<REAL> > > > > EnsemVectorReal;

typedef Ensem< OTensor< PScalar< PScalar< RComplex<REAL> > > > > EnsemTensorComplex;
typedef Ensem< OTensor< PScalar< PScalar< RScalar<REAL> > > > > EnsemTensorReal;

typedef Ensem< OScalar< PScalar< PSeed < RScalar<INTEGER32> > > > > EnsemSeed;
typedef Ensem< OScalar< PScalar< PScalar< RScalar<INTEGER32> > > > > EnsemInteger;
typedef Ensem< OScalar< PScalar< PScalar< RScalar<REAL> > > > > EnsemReal;
typedef Ensem< OScalar< PScalar< PScalar< RScalar<DOUBLE> > > > > EnsemDouble;
typedef Ensem< OScalar< PScalar< PScalar< RScalar<LOGICAL> > > > > EnsemBoolean;

typedef EScalar< OScalar< PSpinVector< PScalar< RComplex<REAL> >, Ns> > > SpinVector;
typedef EScalar< OScalar< PSpinMatrix< PScalar< RComplex<REAL> >, Ns> > > SpinMatrix;
typedef EScalar< OScalar< PScalar< PScalar< RComplex<REAL> > > > > Complex;

typedef EScalar< OVector< PScalar< PScalar< RComplex<REAL> > > > > VectorComplex;
typedef EScalar< OVector< PScalar< PScalar< RScalar<REAL> > > > > VectorReal;

typedef EScalar< OTensor< PScalar< PScalar< RComplex<REAL> > > > > TensorComplex;
typedef EScalar< OTensor< PScalar< PScalar< RScalar<REAL> > > > > TensorReal;

typedef EScalar< OScalar< PScalar< PSeed< RScalar<INTEGER32> > > > > Seed;
typedef EScalar< OScalar< PScalar< PScalar< RScalar<INTEGER32> > > > > Integer;
typedef EScalar< OScalar< PScalar< PScalar< RScalar<REAL> > > > > Real;
typedef EScalar< OScalar< PScalar< PScalar< RScalar<DOUBLE> > > > > Double;
typedef EScalar< OScalar< PScalar< PScalar< RScalar<LOGICAL> > > > > Boolean;

typedef EScalar< OScalar< PSpinVector< PScalar< RComplex<DOUBLE> >, Ns> > > DSpinVector;
typedef EScalar< OScalar< PSpinMatrix< PScalar< RComplex<DOUBLE> >, Ns> > > DSpinMatrix;
typedef EScalar< OScalar< PScalar< PScalar< RComplex<DOUBLE> > > > > DComplex;

typedef EScalar< OScalar< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, Ns> > > DPropagator;

// Level below ensemble for internal convenience
//typedef OScalar< PScalar< PScalar< RScalar<REAL> > > > IntReal;
//typedef OScalar< PScalar< PScalar< RScalar<REAL32> > > > IntReal32;
//typedef OScalar< PScalar< PScalar< RScalar<INTEGER32> > > > IntInteger;
//typedef OScalar< PScalar< PScalar< RScalar<REAL64> > > > IntReal64;
//typedef OScalar< PScalar< PScalar< RScalar<DOUBLE> > > > IntDouble;
//typedef OScalar< PScalar< PScalar< RScalar<LOGICAL> > > > IntBoolean;

// Odd-ball to support random numbers
//typedef Real IEnsemReal;
//typedef Seed IEnsemSeed;

//
// Fixed precision
//
// REAL32 types
typedef Ensem< OScalar< PSpinMatrix< PScalar< RComplex<REAL32> >, Ns> > > EnsemSpinMatrixF;
typedef Ensem< OScalar< PSpinVector< PScalar< RComplex<REAL32> >, Ns> > > EnsemSpinVectorF;

typedef Ensem< OVector< PScalar< PScalar< RComplex<REAL32> > > > > EnsemVectorComplexF;
typedef Ensem< OVector< PScalar< PScalar< RScalar<REAL32> > > > > EnsemVectorRealF;

typedef Ensem< OTensor< PScalar< PScalar< RComplex<REAL32> > > > > EnsemTensorComplexF;
typedef Ensem< OTensor< PScalar< PScalar< RScalar<REAL32> > > > > EnsemTensorRealF;

typedef Ensem< OScalar< PScalar< PScalar< RComplex<REAL32> > > > > EnsemComplexF;
typedef Ensem< OScalar< PScalar< PScalar< RScalar<REAL32> > > > > EnsemRealF;

typedef EScalar< OScalar< PSpinMatrix< PScalar< RComplex<REAL32> >, Ns> > > SpinMatrixF;
typedef EScalar< OScalar< PSpinVector< PScalar< RComplex<REAL32> >, Ns> > > SpinVectorF;

typedef EScalar< OVector< PScalar< PScalar< RComplex<REAL32> > > > > VectorComplexF;
typedef EScalar< OVector< PScalar< PScalar< RScalar<REAL32> > > > > VectorRealF;

typedef EScalar< OTensor< PScalar< PScalar< RComplex<REAL32> > > > > TensorComplexF;
typedef EScalar< OTensor< PScalar< PScalar< RScalar<REAL32> > > > > TensorRealF;

typedef EScalar< OScalar< PScalar< PScalar< RComplex<REAL32> > > > > ComplexF;
typedef EScalar< OScalar< PScalar< PScalar< RScalar<REAL32> > > > > RealF;

// REAL64 types
typedef Ensem< OScalar< PSpinMatrix< PScalar< RComplex<REAL64> >, Ns> > > EnsemSpinMatrixD;
typedef Ensem< OScalar< PSpinVector< PScalar< RComplex<REAL64> >, Ns> > > EnsemSpinVectorD;

typedef Ensem< OVector< PScalar< PScalar< RComplex<REAL64> > > > > EnsemVectorComplexD;
typedef Ensem< OVector< PScalar< PScalar< RScalar<REAL64> > > > > EnsemVectorRealD;

typedef Ensem< OTensor< PScalar< PScalar< RComplex<REAL64> > > > > EnsemTensorComplexD;
typedef Ensem< OTensor< PScalar< PScalar< RScalar<REAL64> > > > > EnsemTensorRealD;

typedef Ensem< OScalar< PScalar< PScalar< RComplex<REAL64> > > > > EnsemComplexD;
typedef Ensem< OScalar< PScalar< PScalar< RScalar<REAL64> > > > > EnsemRealD;

typedef EScalar< OScalar< PSpinMatrix< PScalar< RComplex<REAL64> >, Ns> > > SpinMatrixD;
typedef EScalar< OScalar< PSpinVector< PScalar< RComplex<REAL64> >, Ns> > > SpinVectorD;

typedef EScalar< OTensor< PScalar< PScalar< RComplex<REAL64> > > > > TensorComplexD;
typedef EScalar< OTensor< PScalar< PScalar< RScalar<REAL64> > > > > TensorRealD;

typedef EScalar< OScalar< PScalar< PScalar< RComplex<REAL64> > > > > ComplexD;
typedef EScalar< OScalar< PScalar< PScalar< RScalar<REAL64> > > > > RealD;

// Equivalent names
typedef Integer  Int;

typedef RealF  Real32;
typedef RealD  Real64;
typedef ComplexF  Complex32;
typedef ComplexD  Complex64;

typedef EnsemInteger  EnsemInt;


/*! @} */   // end of group defs

} // namespace ENSEM

