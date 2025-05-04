// -*- C++ -*-
// $Id: formfac_manage_3pt_key.h,v 2.1 2009/03/21 21:33:45 edwards Exp $
/*! \file
 * \brief Key for manage 3-pt funcs
 */

#ifndef __formfac_manage_3pt_key_h__
#define __formfac_manage_3pt_key_h__

#include "ensem/ensem.h"
#include "formfac/formfac_ensemble.h"
#include "io/adat_xmlio.h"
#include "io/adat_io.h"

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  //----------------------------------------------------------------------------
  //! Three pt
  struct ThreePtArg
  {
    ThreePtArg() {gamma=-1;dt=-1;quark=-1;src_spin=-1;snk_spin=-1;}

    PiPf         pi_pf;       /*!< Source and sink momenta */
    int          gamma;       /*!< Gamma matrix index (0 .. 15). In DP basis */
    Array<int>   links;       /*!< Gauge link insertions */

    int          dt;          /*!< Source-sink separation */
    int          quark;       /*!< Some number indicating which quark line */

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
  //! Reader of a threept arg
  void read(XMLReader& xml, const std::string& path, ThreePtArg& param);

  //! Writer of a threept arg
  void write(XMLWriter& xml, const std::string& path, const ThreePtArg& param);

  //----------------------------------------------------------------------------
  //! Reader of a threept arg
  void read(BinaryReader& bin, ThreePtArg& param);

  //! Writer of a threept arg
  void write(BinaryWriter& bin, const ThreePtArg& param);

  //----------------------------------------------------------------------------
  //! Reduced version of a three pt key
  struct ThreePtArgReduced
  {
    ThreePtArgReduced() {gamma=-1;}

    PiPf         pi_pf;      /*!< Source and sink momenta */
    int          gamma;      /*!< Gamma matrix index (0 .. 15). In DP basis */
    Array<int>   links;      /*!< Gauge link insertions */
  };


  //----------------------------------------------------------------------------
  //! Merge this key with another key to produce a final key
  void mergeKeys(ThreePtArg& out, const ThreePtArgReduced& in, const ThreePtArg& prototype);

} // namespace FF

#endif
