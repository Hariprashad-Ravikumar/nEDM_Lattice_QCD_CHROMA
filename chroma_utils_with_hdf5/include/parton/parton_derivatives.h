// -*- C++ -*-
// $Id: parton_derivatives.h,v 2.0 2008/12/05 04:43:38 edwards Exp $
/*! \file
 * \brief Derivatives
 */

#ifndef __derivatives_h__
#define __derivatives_h__

#include "parton/parton_correlation_functions.h"

namespace FF 
{

/*#####################################################################################*/
/* generic quark derivatives                                                                 */
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_D_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu1);
  EnsemVectorComplexF QBar_Gamma_DD_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu2, int mu1);
  EnsemVectorComplexF QBar_Gamma_DDD_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu3, int mu2, int mu1);

} // namespace FF

#endif

