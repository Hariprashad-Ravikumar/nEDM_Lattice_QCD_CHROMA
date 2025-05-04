// -*- C++ -*-
// $Id: formfac_solver_row.h,v 2.0 2008/12/05 04:43:37 edwards Exp $
/*! \file
 * \brief Keys used for solver
 */

#ifndef __formfac_solver_row_h__
#define __formfac_solver_row_h__

#include "adat/handle.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_ensemble.h"
#include <string>

namespace FF
{
  // Namespace composition
  using namespace ADAT;
  using namespace ENSEM;

  //-----------------------------------------------------------------------------------
  //! The row entries for the linear least squares system
  struct TimeFitRange_t
  {
    int t_i;   /*!< Starting time-slice */
    int t_f;   /*!< Ending time-slice */
  };

  //! Reader
  void read(XMLReader& xml, const std::string& path, TimeFitRange_t& param);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const TimeFitRange_t& param);


  //-----------------------------------------------------------------------------------
  //! The row entries for the linear least squares system
  struct LLSqRow_t
  {
    int          t;               /*!< Time-slice */
    PiPf         pi_pf;           /*!< Source and sink momentum */
    Array<int>   insert_lorentz;  /*!< Insertion lorentz indices */

    int          dt;              /*!< Source-sink separation */
    int          quark;           /*!< Some number indicating which quark line */

    std::string  src_name;        /*!< Source operator name */ 
    std::string  src_smear;       /*!< Source smearing */ 
    Array<int>   src_lorentz;     /*!< Source lorentz indices */
    int          src_spin;        /*!< Source Dirac spin indices */

    std::string  snk_name;        /*!< Sink operator name */
    std::string  snk_smear;       /*!< Sink smearing */
    Array<int>   snk_lorentz;     /*!< Sink lorentz indices */
    int          snk_spin;        /*!< Sink Dirac spin indices */

    std::string  mass;            /*!< Some string label for the mass(es) in the corr */
    std::string  ensemble;        /*!< Label for the ensemble */
  };


}  // namespace FF

#endif
