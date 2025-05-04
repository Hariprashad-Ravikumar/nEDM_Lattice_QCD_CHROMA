// -*- C++ -*-
// $Id: ensem_primitive.h,v 2.1 2009/10/27 02:11:14 edwards Exp $

/*! \file
 * \brief Primitive classes
 *
 * Primitives are the various types on the fibers at the lattice sites
 */


/*! \defgroup fiber Fiber only types and operations
 * \ingroup fiberbundle
 *
 * Primitives are the various types on the fibers at the lattice sites.
 *
 * The primitive indices, including Reality (also known as complex or real),
 * is represented as a tensor product over various vector spaces. Different
 * kinds of object can transform in those vector spaces, like Scalar, Vector, and
 * Matrix.
 */

#include "ensem_primscalar.h"
#include "ensem_primmatrix.h"
#include "ensem_primvector.h"
#include "ensem_primseed.h"
#include "ensem_primcolormat.h"
#include "ensem_primcolorvec.h"
#include "ensem_primgamma.h"
#include "ensem_primspinmat.h"
#include "ensem_primspinvec.h"

