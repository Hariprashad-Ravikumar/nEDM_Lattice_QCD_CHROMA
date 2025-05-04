// -*- C++ -*-
// $Id: ensem_precision.h,v 2.0 2008/12/05 04:43:34 edwards Exp $

/*! \file
 * \brief PRECISION ISSUES
 */

#ifndef ENSEM_PRECISION_H
#define ENSEM_PRECISION_H

// Fix default precision
#if ! defined(BASE_PRECISION)
#define BASE_PRECISION 64
#endif

// These are fixed precision versions
typedef int       INTEGER32;
typedef float     REAL32;
typedef double    REAL64;
typedef bool      LOGICAL;

// Set the base floating precision
#if BASE_PRECISION == 32
// Use single precision for base precision
typedef REAL32    REAL;
typedef REAL64    DOUBLE;

#elif BASE_PRECISION == 64
// Use double precision for base precision
typedef REAL64    REAL;
typedef REAL64    DOUBLE;

#else
#error "Unknown BASE_PRECISION"
#endif

#endif
