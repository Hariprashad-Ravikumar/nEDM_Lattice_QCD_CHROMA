// -*- C++ -*-
// $Id: formfac_solver_driver.h,v 2.0 2008/12/05 04:43:37 edwards Exp $
/*! \file
 * \brief Driver for call the LLSq solver
 */

#ifndef __formfac_solver_driver_h__
#define __formfac_solver_driver_h__

#include <list> 
#include "formfac/formfac_solver.h"
#include "formfac/formfac_qsq.h"
#include "io/adat_xmlio.h"

namespace FF
{
  //! Driver
  void llsq_driver(XMLWriter& xml_out, LLSqComponent<EnsemReal>& sys, const std::list<QsqVal>& qsq_val,
		   double cov_tol, double solve_tol);

} // namespace FF

#endif
