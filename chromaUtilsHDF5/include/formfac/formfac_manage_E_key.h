// -*- C++ -*-
// $Id: formfac_manage_E_key.h,v 2.1 2009/03/21 21:33:45 edwards Exp $
/*! \file
 * \brief Keys for manage energy funcs
 */

#ifndef __formfac_manage_E_key_h__
#define __formfac_manage_E_key_h__

#include "ensem/ensem.h"
#include "io/adat_xmlio.h"
#include "io/adat_io.h"

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  //----------------------------------------------------------------------------
  //! Energy key
  struct EnergyArg
  {
    int          level;       /*!< Energy level */
    Array<int>   mom;         /*!< Momentum */

    std::string  mass;        /*!< Some string label for the mass(es) in the corr */
    std::string  ensemble;    /*!< Label for the ensemble */
  };

  //----------------------------------------------------------------------------
  //! Reader of a energy arg
  void read(XMLReader& xml, const std::string& path, EnergyArg& param);

  //! Writer of a energy arg
  void write(XMLWriter& xml, const std::string& path, const EnergyArg& param);

  //----------------------------------------------------------------------------
  //! Reader of a energy arg
  void read(BinaryReader& bin, EnergyArg& param);

  //! Writer of a energy arg
  void write(BinaryWriter& bin, const EnergyArg& param);

} // namespace FF

#endif
