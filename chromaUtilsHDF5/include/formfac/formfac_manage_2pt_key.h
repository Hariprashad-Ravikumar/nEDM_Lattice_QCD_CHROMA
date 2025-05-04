// -*- C++ -*-
// $Id: formfac_manage_2pt_key.h,v 2.1 2009/03/21 21:33:45 edwards Exp $

/*! \file
 * \brief Keys for manage 2-pt funcs
 */

#ifndef __formfac_manage_2pt_key_h__
#define __formfac_manage_2pt_key_h__

#include "ensem/ensem.h"
#include "io/adat_xmlio.h"
#include "io/adat_io.h"

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  //----------------------------------------------------------------------------
  //! Two pt
  struct TwoPtArg
  {
    Array<int>   mom;         /*!< Momentum */

    std::string  src_name;    /*!< Some string label for the operator */
    std::string  src_smear;   /*!< Some string label for the smearing of this operator */
    Array<int>   src_lorentz; /*!< Source Lorentz indices */
    int          src_spin;    /*!< Source Dirac spin indices */

    std::string  snk_name;    /*!< Some string label for the operator */
    std::string  snk_smear;   /*!< Some string label for the smearing of this operator */
    Array<int>   snk_lorentz; /*!< Sink Lorentz indices */
    int          snk_spin;    /*!< Sink Dirac spin indices */

    std::string  mass;        /*!< Some string label for the mass(es) in the corr */
    std::string  ensemble;    /*!< Label for the ensemble */
  };

  //----------------------------------------------------------------------------
  //! Reader of a twopt arg
  void read(XMLReader& xml, const std::string& path, TwoPtArg& param);

  //! Writer of a twopt arg
  void write(XMLWriter& xml, const std::string& path, const TwoPtArg& param);

  //----------------------------------------------------------------------------
  //! Reader of a twopt arg
  void read(BinaryReader& bin, TwoPtArg& param);

  //! Writer of a twopt arg
  void write(BinaryWriter& xml, const TwoPtArg& param);

} // namespace FF

#endif
