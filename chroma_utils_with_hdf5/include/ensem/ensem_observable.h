// -*- C++ -*-
// $Id: ensem_observable.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! \file
 * \brief Observable classes
 *
 * Observables are the various types on the fibers at the lattice sites
 */


/*! \defgroup fiber Fiber only types and operations
 * \ingroup fiberbundle
 *
 * Observables are the various types on the fibers at the lattice sites.
 *
 * The observable indices, including Reality (also known as complex or real),
 * is represented as a tensor product over various vector spaces. Different
 * kinds of object can transform in those vector spaces, like Scalar, Vector, and
 * Matrix.
 */

#include "ensem_obsscalar.h"
//#include "ensem_obsmatrix.h"
#include "ensem_obsvector.h"
#include "ensem_obstensor.h"

