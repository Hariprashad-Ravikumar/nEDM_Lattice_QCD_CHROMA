// -*- C++ -*-
// $Id: formfac_chisqq.h,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Chisq
 */

#ifndef __formfac_chisqq_h__
#define __formfac_chisqq_h__

namespace FF
{
  // ********************** FUNCTION CHISQQ *****************************
  //! Returns (logically)  gammq(ndof/2,chisq/2)
  double chisqq(int n, double x);

  // Same as "chisqq", but will check if ndof is <=0 and return 0.
  double chisqq_zero(int ndof, double x);

} // end namespace

#endif
