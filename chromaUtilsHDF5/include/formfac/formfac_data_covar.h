// -*- C++ -*-
// $Id: formfac_data_covar.h,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Covariance and inverse of covariance
 */

#ifndef __formfac_data_covar_h__
#define __formfac_data_covar_h__

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include <vector>

namespace FF
{
  using namespace ENSEM;

  //! Covariance on real data
  Array2d<Real> data_covar(const std::vector<ENSEM::EnsemReal>& data);

  //struct to hold the inverse covariance matrix and the number of reset sing values
  //also holds the info about the inverted range
  struct data_covar_inv_t
  {
    Array2d<Real> inv_cov;
    
    int n_reset_sv;
    int Tmin;
    int Tmax;
  };

  //! Invariance of the covariance matrix
  //Array2d<Real> data_covar_inv(const Array2d<Real>& covar, double tol);
  data_covar_inv_t data_covar_inv(const Array2d<Real>& covar, double tol);
  
  //! Invariance of a subset of the covariance matrix
  /*! 
   * Only a subset of the matrix is inverted. However, the original size matrix
   * is returned to save on data type manipulation.
   */
  // Array2d<Real> data_covar_inv(const Array2d<Real>& covar, double tol, int Tmin, int Tmax);
  data_covar_inv_t data_covar_inv(const Array2d<Real>& covar, double tol, int Tmin, int Tmax);





  //! Write correlation matrix
  void writeCorrMat(XMLWriter& xml_out, const std::string& path, const Array2d<Real>& covar);

  //! Write correlation matrix times its inverse as a test
  void writeCorrMatTest(XMLWriter& xml_out, const std::string& path, 
			const Array2d<Real>& covar,
			const Array2d<Real>& inv_covar,
			int Tmin, int Tmax);
} // end namespace

#endif
