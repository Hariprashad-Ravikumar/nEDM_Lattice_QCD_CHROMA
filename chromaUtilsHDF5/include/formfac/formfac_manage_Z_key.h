// -*- C++ -*-
// $Id: formfac_manage_Z_key.h,v 2.1 2009/03/21 21:33:45 edwards Exp $
/*! \file
 * \brief Manage amplitudes
 */

#ifndef __formfac_manage_Z_key_h__
#define __formfac_manage_Z_key_h__

#include "ensem/ensem.h"
#include "io/adat_xmlio.h"
#include "io/adat_io.h"

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  //----------------------------------------------------------------------------
  //! Amp key
  struct AmpArg
  {
    int          level;       /*!< Amp level */
    Array<int>   mom;         /*!< Momentum */

    std::string  name;        /*!< Some string label for the operator */
    std::string  smear;       /*!< Some string label for the smearing of this operator */

    std::string  mass;        /*!< Some string label for the mass(es) in the corr */
    std::string  ensemble;    /*!< Label for the ensemble */
  };

  //----------------------------------------------------------------------------
  //! Reader of a amp arg
  void read(XMLReader& xml, const std::string& path, AmpArg& param);

  //! Writer of a amp arg
  void write(XMLWriter& xml, const std::string& path, const AmpArg& param);

  //----------------------------------------------------------------------------
  //! Reader of a amp arg
  void read(BinaryReader& xml, AmpArg& param);

  //! Writer of a amp arg
  void write(BinaryWriter& xml, const AmpArg& param);

} // namespace FF

#endif
