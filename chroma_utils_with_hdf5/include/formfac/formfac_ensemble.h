// -*- C++ -*-
// $Id: formfac_ensemble.h,v 2.0 2008/12/05 04:43:35 edwards Exp $

/*! \file
 * \brief Ensemble and lattice info
 */

#ifndef __formfac_ensemble_h__
#define __formfac_ensemble_h__

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  typedef Array<int>  ArrayInt;   // save some typing
  typedef Array<double> ArrayDouble;  // save some typing

  //! Hold momenta
  struct PiPf
  {
    ArrayInt  p_i;
    ArrayInt  p_f;
  };

  //! Support for maps
  bool operator<(const PiPf& a, const PiPf& b);

  //! Reader
  void read(XMLReader& xml, const std::string& path, PiPf& val);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const PiPf& val);


  //! Anisotropy params
  struct AnisoParam
  {
    AnisoParam() {anisoP = false; c_sq = xi = xi_0 = nu = 1.0; m0_at = 0.0;}
    bool           anisoP;
    double         c_sq;
    double         xi;
    double         xi_0;
    double         nu;
    double         m0_at;
  };

  //! Read params
  void read(XMLReader& xml, const std::string& path, AnisoParam& param);

  //! Write params
  void write(XMLWriter& xml, const std::string& path, const AnisoParam& param);


  //! Lattice params
  struct LatticeParam
  {
    ArrayInt     latt_size;
    int          decay_dir;
    AnisoParam   aniso;
    double       a_fm;  // spatial lattice spacing in fm
  };

  //! Read params
  void read(XMLReader& xml, const std::string& path, LatticeParam& param);

  //! Write params
  void write(XMLWriter& xml, const std::string& path, const LatticeParam& param);


  //! Ensemble Info
  struct EnsembleInfo
  {
    EnsembleInfo() {nbin=0;}
      
    int                 nbin;
    LatticeParam        lattice;
  };

  //! Reader
  void read(XMLReader& xml, const std::string& path, EnsembleInfo& param);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const EnsembleInfo& param);


} // namespace FF

#endif
