// -*- C++ -*-
// $Id: parton_correlation_functions.h,v 2.0 2008/12/05 04:43:37 edwards Exp $
/*! \file
 * \brief Correlation functions
 */

#ifndef __correlation_functions_h__
#define __correlation_functions_h__

#include "formfac/formfac_manage_3pt.h"

namespace FF 
{

/*#####################################################################################*/
/* index types                                                                         */
/*#####################################################################################*/

//  typedef unsigned short int dirac;      /* 0, 1, ..., 15             */
//  typedef   signed short int time;       /* ..., -2, -1, 0, 1, 2, ... */
//  typedef unsigned short int direction;  /* 0, 1, ..., 7              */

/*#####################################################################################*/
/* correlation functions                                                               */
/*#####################################################################################*/

/* momentum is always inserted in the quark field */

/* right to left indexing */

  //! \f$\bar{q}\Gamma q\f$
  EnsemVectorComplexF QBar_Gamma_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma);

  //! \f$\bar{q}\Gamma U_{\mu_1}^\dag q\f$
  EnsemVectorComplexF QBar_Gamma_UDag_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu1);

  //! \f$\bar{q}\Gamma U_{\mu_2}^\dag U^\dag_{\mu_1} q\f$
  EnsemVectorComplexF QBar_Gamma_UDagUDag_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu2, int mu1);

  //! \f$\bar{q}\Gamma U_{\mu_3}^\dag U_{\mu_2}^\dag U^\dag_{\mu_1} q\f$
  EnsemVectorComplexF QBar_Gamma_UDagUDagUDag_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu3, int mu2, int mu1);

  //! \f$\bar{q}\Gamma U_{\mu_4}^\dag U_{\mu_3}^\dag U_{\mu_2}^\dag U^\dag_{\mu_1} q\f$
  EnsemVectorComplexF QBar_Gamma_UDagUDagUDagUDagUDag_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu5, int mu4, int mu3, int mu2, int mu1);

  EnsemVectorComplexF QBar_Gamma_UDagUDagUDagUDagUDag_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu5, int mu4, int mu3, int mu2, int mu1);


/* left to right indexing */

  //! \f$\bar{q} U_{\mu_1} \Gamma q\f$
  EnsemVectorComplexF QBar_U_Gamma_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu1);

  //! \f$\bar{q} U_{\mu_1} U_{\mu_2} \Gamma q\f$
  EnsemVectorComplexF QBar_UU_Gamma_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu1, int mu2);

  //! \f$\bar{q} U_{\mu_1} U_{\mu_2} U_{\mu_3} \Gamma q\f$
  EnsemVectorComplexF QBar_UUU_Gamma_Q(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int gamma, int mu1, int mu2, int mu3);
} // namespace FF  

#endif

