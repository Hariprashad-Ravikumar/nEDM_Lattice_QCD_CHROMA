// -*- C++ -*-
// $Id: formfac_linear_covar_genfit.h,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Chi-sq fit to a sum of linear params with generic function weights
 */

#ifndef __formfac_linear_covar_genfit_h__
#define __formfac_linear_covar_genfit_h__

#include "ensem/ensem.h"

using namespace ENSEM;

namespace FF
{

  //! Data is a single-valued function dependent on a vector array
  struct CorrLinGenFitData_t
  {
    Array<Real>  X;
    Real         f;
  };


  //! Data is a single-valued function dependent on a vector array
  struct CorrLinGenFitResults_t
  {
    struct Coeff_t
    {
      Real a;
      Real a_err;
    };

    Array<Coeff_t>  coeffs;
    Array2d<Real> coeff_cov;

    int  Tmin;
    int  Tmax;

    int  ndof;
    int  count_singular_val;

    Real Q;
    Real chisq;
  };

  //!  Write the results
  void write(XMLWriter& xml_out, const std::string& xml_group, const CorrLinGenFitResults_t& p);


  //!  Fit to the single value functions  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  CorrLinGenFitResults_t corrLinGenFit(const Array<CorrLinGenFitData_t>& data,
				       const Array2d<Real>& inv_covar,
				       double tol,
				       int Tmin, int Tmax);

  //!  Fit to the single value functions  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  CorrLinGenFitResults_t corrLinGenFit(const Array<CorrLinGenFitData_t>& data,
				       const Array2d<Real>& inv_covar,
				       double tol);

  //!  Fit to the single value functions  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  CorrLinGenFitResults_t corrLinGenFit(XMLWriter& xml_out, const std::string& xml_group,
				       const Array<CorrLinGenFitData_t>& data,
				       const Array2d<Real>& inv_covar,
				       double tol);

  //!  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  CorrLinGenFitResults_t corrLinGenFit(XMLWriter& xml_out, const std::string& xml_group,
				       const Array<CorrLinGenFitData_t>& data,
				       const Array2d<Real>& inv_covar,
				       double tol,
				       int Tmin, int Tmax);

} // end namespace

#endif
