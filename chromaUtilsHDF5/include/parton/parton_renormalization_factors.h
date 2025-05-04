// -*- C++ -*-
// $Id: parton_renormalization_factors.h,v 2.0 2008/12/05 04:43:38 edwards Exp $
/*! \file
 * \brief Renormalization factors
 */

#ifndef __renormalization_factors_h__
#define __renormalization_factors_h__

#include "formfac/formfac_manage_3pt.h"

namespace FF
{
  extern double   Z_x_q_a;
  extern double   Z_x_q_b;
  extern double  Z_xx_q;
  extern double Z_xxx_q;

  extern double  Z_1_Delta_q;
  extern double  Z_x_Delta_q_a;
  extern double  Z_x_Delta_q_b;
  extern double Z_xx_Delta_q;

  extern double Z_1_delta_q;
  extern double Z_x_delta_q;

  extern double Z_d1_q;
  extern double Z_d2_q;

  void SetRenormalizationFactors( const float beta, const float a ); // [a] = GeV

  void SetRenormalizationFactors_DD( const float beta, const float a ); // [a] = GeV

} // namespace FF

#endif

