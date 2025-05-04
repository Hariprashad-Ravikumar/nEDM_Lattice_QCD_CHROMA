// -*- C++ -*-
// $Id: parton_distribution_moments.h,v 2.0 2008/12/05 04:43:38 edwards Exp $
/*! \file
 * \brief Parton distribution moments
 */

#ifndef __parton_distribution_moments_h__
#define __parton_distribution_moments_h__

#include "formfac/formfac_manage_3pt.h"

namespace FF 
{

  EnsemVectorRealF q_d1(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_d2(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_lpol_0(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_lpol_1_a(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_lpol_1_b(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_lpol_2(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_tpol_1(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_upol_0(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_upol_1_a(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_upol_1_b(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_upol_2(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

  EnsemVectorRealF q_upol_3(Manage3PtFuncReduced& threept, const PiPf& pi_pf);

} // namespace FF

#endif

