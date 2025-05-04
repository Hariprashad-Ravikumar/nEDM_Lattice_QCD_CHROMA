// -*- C++ -*-
// $Id: formfac_ensem_linfit.h,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Chi-sq fit to a sum of linear params with generic function weights
 */

#ifndef __formfac_ensem_linfit_h__
#define __formfac_ensem_linfit_h__

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"

using namespace ENSEM;

namespace FF
{

  //-----------------------------------------------------------------------------------
  //! Data is a single-valued function dependent on a vector array
  struct EnsemCorrLinGenFitData_t
  {
    Array<EnsemReal>  X;
    EnsemReal         f;
  };


  //! Data is a single-valued function dependent on a vector array
  struct EnsemCorrLinGenFitResults_t
  {
    Array<EnsemReal>  coeffs;
    
    Array2d<EnsemReal> coeff_cov;

    int  ndof;
    int  count_singular_val;

    Real Q;
    Real chisq;
  };


  //  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  EnsemCorrLinGenFitResults_t corrLinGenFit(const Array<EnsemCorrLinGenFitData_t>& data,
					    const Array2d<Real>& inv_covar,
					    double tol);

  //  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  EnsemCorrLinGenFitResults_t corrLinGenFit(const Array<EnsemCorrLinGenFitData_t>& data,
					    const Array2d<Real>& inv_covar,
					    double tol,
					    int Tmin, int Tmax);

} // end namespace

#endif
