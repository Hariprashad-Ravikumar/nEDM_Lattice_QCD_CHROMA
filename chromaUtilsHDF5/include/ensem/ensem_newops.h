// -*- C++ -*-
// $Id: ensem_newops.h,v 2.1 2009/10/27 02:11:14 edwards Exp $

/*! @file
 * @brief Additional operations on EnsemTypes
 */

namespace ENSEM {


//-----------------------------------------------------------------------------
// Operator tags that are only used for type resolution
//-----------------------------------------------------------------------------

struct FnMean
{
  PETE_EMPTY_CONSTRUCTORS(FnMean)
};

struct FnVariance
{
  PETE_EMPTY_CONSTRUCTORS(FnVariance)
};

struct FnContract
{
  PETE_EMPTY_CONSTRUCTORS(FnContract)
};

struct FnQuarkContractXX
{
  PETE_EMPTY_CONSTRUCTORS(FnQuarkContractXX)
};

struct FnSum
{
  PETE_EMPTY_CONSTRUCTORS(FnSum)
};

struct FnNorm2
{
  PETE_EMPTY_CONSTRUCTORS(FnNorm2)
};

struct FnInnerProduct
{
  PETE_EMPTY_CONSTRUCTORS(FnInnerProduct)
};

struct FnInnerProductReal
{
  PETE_EMPTY_CONSTRUCTORS(FnInnerProductReal)
};

struct FnAdjoint
{
  PETE_EMPTY_CONSTRUCTORS(FnAdjoint)
};

struct FnConjugate
{
  PETE_EMPTY_CONSTRUCTORS(FnConjugate)
};

struct FnTranspose
{
  PETE_EMPTY_CONSTRUCTORS(FnTranspose)
};

struct FnTransposeColor
{
  PETE_EMPTY_CONSTRUCTORS(FnTransposeColor)
};

struct FnTransposeSpin
{
  PETE_EMPTY_CONSTRUCTORS(FnTransposeSpin)
};

struct FnTrace
{
  PETE_EMPTY_CONSTRUCTORS(FnTrace)
};

struct FnRealTrace
{
  PETE_EMPTY_CONSTRUCTORS(FnRealTrace)
};

struct FnImagTrace
{
  PETE_EMPTY_CONSTRUCTORS(FnImagTrace)
};

struct FnTraceColor
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceColor)
};

struct FnTraceSpin
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceSpin)
};

struct FnReal
{
  PETE_EMPTY_CONSTRUCTORS(FnReal)
};

struct FnImag
{
  PETE_EMPTY_CONSTRUCTORS(FnImag)
};

struct FnLocalNorm2
{
  PETE_EMPTY_CONSTRUCTORS(FnLocalNorm2)
};

struct FnTimesI
{
  PETE_EMPTY_CONSTRUCTORS(FnTimesI)
};

struct FnTimesMinusI
{
  PETE_EMPTY_CONSTRUCTORS(FnTimesMinusI)
};

struct FnSeedToFloat
{
  PETE_EMPTY_CONSTRUCTORS(FnSeedToFloat)
};

struct FnChiralProjectPlus
{
  PETE_EMPTY_CONSTRUCTORS(FnChiralProjectPlus)
};

struct FnChiralProjectMinus
{
  PETE_EMPTY_CONSTRUCTORS(FnChiralProjectMinus)
};

struct FnCmplx
{
  PETE_EMPTY_CONSTRUCTORS(FnCmplx)
};

struct FnOuterProduct
{
  PETE_EMPTY_CONSTRUCTORS(FnOuterProduct)
};

struct FnColorVectorContract
{
  PETE_EMPTY_CONSTRUCTORS(FnColorVectorContract)
};

struct FnColorCrossProduct
{
  PETE_EMPTY_CONSTRUCTORS(FnColorCrossProduct)
};

struct FnLocalInnerProduct
{
  PETE_EMPTY_CONSTRUCTORS(FnLocalInnerProduct)
};

struct FnLocalInnerProductReal
{
  PETE_EMPTY_CONSTRUCTORS(FnLocalInnerProductReal)
};

struct FnShift
{
  PETE_EMPTY_CONSTRUCTORS(FnShift)
};

struct FnCshift
{
  PETE_EMPTY_CONSTRUCTORS(FnCshift)
};

struct FnQuarkContract13
{
  PETE_EMPTY_CONSTRUCTORS(FnQuarkContract13)
};

struct FnQuarkContract14
{
  PETE_EMPTY_CONSTRUCTORS(FnQuarkContract14)
};

struct FnQuarkContract23
{
  PETE_EMPTY_CONSTRUCTORS(FnQuarkContract23)
};

struct FnQuarkContract24
{
  PETE_EMPTY_CONSTRUCTORS(FnQuarkContract24)
};

struct FnQuarkContract12
{
  PETE_EMPTY_CONSTRUCTORS(FnQuarkContract12)
};

struct FnQuarkContract34
{
  PETE_EMPTY_CONSTRUCTORS(FnQuarkContract34)
};

struct FnColorContract
{
  PETE_EMPTY_CONSTRUCTORS(FnColorContract)
};


//-----------------------------------------------------------------------------
// Optimization hooks
//-----------------------------------------------------------------------------

struct OpAdjMultiply
{
  PETE_EMPTY_CONSTRUCTORS(OpAdjMultiply)
};

struct OpMultiplyAdj
{
  PETE_EMPTY_CONSTRUCTORS(OpMultiplyAdj)
};

struct OpAdjMultiplyAdj
{
  PETE_EMPTY_CONSTRUCTORS(OpAdjMultiplyAdj)
};

struct FnTraceMultiply
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceMultiply)
};

struct FnTraceColorMultiply
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceColorMultiply)
};

struct FnTraceSpinMultiply
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceSpinMultiply)
};

//-----------------------------------------------------------------------------
// Operators and tags for accessing elements of a Ensem object
//-----------------------------------------------------------------------------

//! Structure for extracting observable vector components
struct FnPeekEnsem
{
  PETE_EMPTY_CONSTRUCTORS(FnPeekEnsem)
};

//! Structure for extracting observable vector components
struct FnPeekObsVector
{
  PETE_EMPTY_CONSTRUCTORS(FnPeekObsVector)
};

//! Structure for extracting observable tensor components
struct FnPeekObsTensor
{
  PETE_EMPTY_CONSTRUCTORS(FnPeekObsTensor)
};

//! Structure for extracting color matrix components
struct FnPeekColorMatrix
{
  PETE_EMPTY_CONSTRUCTORS(FnPeekColorMatrix)
};

//! Structure for extracting color vector components
struct FnPeekColorVector
{
  PETE_EMPTY_CONSTRUCTORS(FnPeekColorVector)
};

//! Structure for extracting spin matrix components
struct FnPeekSpinMatrix
{
  PETE_EMPTY_CONSTRUCTORS(FnPeekSpinMatrix)
};

//! Structure for extracting spin vector components
struct FnPeekSpinVector
{
  PETE_EMPTY_CONSTRUCTORS(FnPeekSpinVector)
};

//! Structure for inserting observable vector components
struct FnPokeEnsem
{
  PETE_EMPTY_CONSTRUCTORS(FnPokeEnsem)
};

//! Structure for inserting observable vector components
struct FnPokeObsVector
{
  PETE_EMPTY_CONSTRUCTORS(FnPokeObsVector)
};

//! Structure for inserting observable tensor components
struct FnPokeObsTensor
{
  PETE_EMPTY_CONSTRUCTORS(FnPokeObsTensor)
};

//! Structure for inserting color matrix components
struct FnPokeColorMatrix
{
  PETE_EMPTY_CONSTRUCTORS(FnPokeColorMatrix)
};


//! Structure for inserting color vector components
struct FnPokeColorVector
{
  PETE_EMPTY_CONSTRUCTORS(FnPokeColorVector)
};

//! Structure for inserting spin matrix components
struct FnPokeSpinMatrix
{
  PETE_EMPTY_CONSTRUCTORS(FnPokeSpinMatrix)
};

//! Structure for inserting spin vector components
struct FnPokeSpinVector
{
  PETE_EMPTY_CONSTRUCTORS(FnPokeSpinVector)
};


//-----------------------------------------------------------------------------
// Additional operator tags 
//-----------------------------------------------------------------------------

struct OpGammaConstMultiply
{
  PETE_EMPTY_CONSTRUCTORS(OpGammaConstMultiply)
};


struct OpMultiplyGammaConst
{
  PETE_EMPTY_CONSTRUCTORS(OpMultiplyGammaConst)
};


// Member function definition in primgamma.h
struct OpGammaTypeMultiply
{
  PETE_EMPTY_CONSTRUCTORS(OpGammaTypeMultiply)
};


// Member function definition in primgamma.h
struct OpMultiplyGammaType
{
  PETE_EMPTY_CONSTRUCTORS(OpMultiplyGammaType)
};


//-----------------------------------------------------------------------------
// Additional operator tags 
//-----------------------------------------------------------------------------

struct OpGammaConstDPMultiply
{
  PETE_EMPTY_CONSTRUCTORS(OpGammaConstDPMultiply)
};


struct OpMultiplyGammaConstDP
{
  PETE_EMPTY_CONSTRUCTORS(OpMultiplyGammaConstDP)
};


// Member function definition in primgamma.h
struct OpGammaTypeDPMultiply
{
  PETE_EMPTY_CONSTRUCTORS(OpGammaTypeDPMultiply)
};


// Member function definition in primgamma.h
struct OpMultiplyGammaTypeDP
{
  PETE_EMPTY_CONSTRUCTORS(OpMultiplyGammaTypeDP)
};


} // namespace ENSEM

