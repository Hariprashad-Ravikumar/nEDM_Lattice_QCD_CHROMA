// -*- C++ -*-
// $Id: ensem.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! \file
 * \brief Primary include file for ENSEM ensbc
 *
 * No other file should be included by the user
 */


// #define ENSEM_DEBUG 3


#ifndef ENSEM_INCLUDE
#define ENSEM_INCLUDE

#include <xml_array.h>
#include <xml_array2d.h>
#include "io/adat_xmlio.h"
#include "io/adat_io.h"

namespace ENSEM
{

  // Namespace composition
  using namespace ADATXML;
  using namespace ADATIO;
  using XMLArray::Array;
  using XMLArray::Array2d;

  // Get local precision info
#include "ensem_precision.h"

  // GNU specific stuff
#if defined(__GNUC__)
// Under g++, enforce using V3 or greater
#if __GNUC__ < 3
#error "ENSEM++ requires g++ 3.0 or higher. This version of the g++ compiler is not supported"
#endif
#endif

  // Under gcc, set some attributes
#if defined(__GNUC__)
// gcc
#define ENSEM_ALIGN8   __attribute__ ((aligned (8)))
#define ENSEM_ALIGN16  __attribute__ ((aligned (16)))
#define ENSEM_INLINE   __attribute__ ((always_inline))
// The attributes in ENSEM_CONST is buggering g++-3.4 
// #define ENSEM_CONST    __attribute__ ((const,pure))
#define ENSEM_CONST
#define ENSEM_CINLINE  __attribute__ ((always_inline,const,pure))
#else
// default
#define ENSEM_ALIGN8
#define ENSEM_ALIGN16
#define ENSEM_INLINE
#define ENSEM_CONST
#define ENSEM_CINLINE
#endif

#define ENSEM_ALIGNMENT_SIZE  16

//---------------------------------------------------------------
// Snarfed bits from PETE
#if defined(PETE_MAKE_EMPTY_CONSTRUCTORS)

#define PETE_EMPTY_CONSTRUCTORS(CLASS)  \
  CLASS() { }   \
  CLASS(const CLASS &) { } \
  CLASS &operator=(const CLASS &) { return *this; }

#define PETE_EMPTY_CONSTRUCTORS_TEMPLATE(CLASS, ARG)  \
  CLASS() { }   \
  CLASS(const CLASS<ARG> &) { } \
  CLASS &operator=(const CLASS<ARG> &) { return *this; }

#else

#define PETE_EMPTY_CONSTRUCTORS(CLASS)
#define PETE_EMPTY_CONSTRUCTORS_TEMPLATE(CLASS, ARG)

#endif

#include "ensem_type_computations.h"
#include "ensem_operator_tags.h"
//---------------------------------------------------------------
} // namespace ENSEM

#include "ensem_forward.h"

#include "ensem_params.h"
#include "ensem_io.h"

#include "ensem_traits.h"

#include "ensem_newops.h"
#include "ensem_simpleword.h"
#include "ensem_reality.h"
#include "ensem_primitive.h"
#include "ensem_observable.h"
#include "ensem_ensem.h"

#include "ensem_scalarsite_defs.h"
#include "ensem_specializations.h"

#include "ensem_random.h"

#include "ensem_scalar_specific.h"

#endif  // ENSEM_INCLUDE
