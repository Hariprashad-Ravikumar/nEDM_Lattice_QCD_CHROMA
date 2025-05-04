// -*- C++ -*-
// $Id: hadron_3pt_corr.h,v 2.4 2009/03/21 21:33:45 edwards Exp $
/*! \file
 * \brief Hadron 3pt correlators
 */

#ifndef __hadron_3pt_corr_h__
#define __hadron_3pt_corr_h__

#include "ensem/ensem.h"
#include "adat/map_obj.h"
#include "io/key_val_db.h"
#include "ConfDataStoreDB.h"

#include "formfac/formfac_ensemble.h"
#include "formfac/formfac_manage_3pt_key.h"

#include <iostream>

namespace FF
{
  using namespace ENSEM;
  using namespace ADATIO;
  using namespace ADAT;

  //----------------------------------------------------------------------------
  //! Key for Hadron 3pt corr
  struct KeyHadron3PtCorr_t
  {
    int          num_vecs;    /*!< Number of vectors used in this corr */

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
  //! Used for error output
  std::ostream& operator<<(std::ostream& os, const KeyHadron3PtCorr_t& d);
	
  //----------------------------------------------------------------------------
  //! KeyHadron3PtCorr reader
  void read(XMLReader& xml, const std::string& path, KeyHadron3PtCorr_t& param);

  //! KeyHadron3PtCorr writer
  void write(XMLWriter& xml, const std::string& path, const KeyHadron3PtCorr_t& param);

  //----------------------------------------------------------------------------
  //! KeyHadron3PtCorr reader
  void read(BinaryReader& bin, KeyHadron3PtCorr_t& param);

  //! KeyHadron3PtCorr writer
  void write(BinaryWriter& bin, const KeyHadron3PtCorr_t& param);

  //----------------------------------------------------------------------------
  //! Merge this key with another key to produce a final key
  void mergeKeys(KeyHadron3PtCorr_t& out, const ThreePtArg& in, const KeyHadron3PtCorr_t& prototype);

  //----------------------------------------------------------------------------
  //! Use an in-memory map
  typedef MapObject<KeyHadron3PtCorr_t,EnsemVectorComplex>  MapHadron3PtCorr_t;

  //! Use a DB
  typedef FILEDB::ConfDataStoreDB< SerialDBKey<KeyHadron3PtCorr_t>, SerialDBData<VectorComplex> >  DBHadron3PtSingleCorr_t;

  //! Use a DB
  typedef FILEDB::ConfDataStoreDB< SerialDBKey<KeyHadron3PtCorr_t>, SerialDBData<EnsemVectorComplex> >  DBHadron3PtCorr_t;

} // namespace ColorVec

#endif
