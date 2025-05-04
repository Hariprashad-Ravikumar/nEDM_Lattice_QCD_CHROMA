// -*- C++ -*-
// $Id: formfac_util.h,v 2.0 2008/12/05 04:43:37 edwards Exp $

/*! \file
 * \brief Utilities
 */

#ifndef __formfac_util_h__
#define __formfac_util_h__

#include "ensem/ensem.h"
#include "formfac/formfac_qsq.h"

namespace FF
{

  //! Compute energy via dispersion relation
  EnsemReal energyDisp(const EnsemReal& m, const ArrayInt& p, const LatticeParam& lattice);

  //! special photon energy routine
  EnsemReal photonenergyDisp(const EnsemReal& Qsq, const ArrayInt& p, const LatticeParam& lattice);

  //! Decompose a lexicographic site into coordinates
  Array<int> crtesn(int ipos, const Array<int>& latt_size);

  //! Calculates the lexicographic site index from the coordinate of a site
  /*! 
   * Nothing specific about the actual lattice size, can be used for 
   * any kind of latt size 
   */
  int local_site(const Array<int>& coord, const Array<int>& latt_size);

} // namespace FF

#endif
