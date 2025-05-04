// -*- C++ -*-
/*! \file
 * \brief Hadron 1pt correlators
 */

#ifndef __hadron_1pt_corr_h__
#define __hadron_1pt_corr_h__

#include "ensem/ensem.h"
#include "adat/map_obj.h"
#include "io/key_val_db.h"
#include "ConfDataStoreDB.h"

#include <iostream>

namespace FF
{
  using namespace ENSEM;
  using namespace ADATIO;
  using namespace ADAT;

  //----------------------------------------------------------------------------
  //! Key for Hadron 1pt corr
  struct KeyHadron1PtCorr_t
  {
    int          num_vecs;    /*!< Number of vectors used in this corr */

    int          creation_op; /*!< Is this a creation ops?  */

    std::string  name;        /*!< Some string label for the operator */
    std::string  smear;       /*!< Some string label for the smearing of this operator */
    Array<int>   lorentz;     /*!< Source Lorentz indices */
    int          spin;        /*!< Source Dirac spin indices */

    Array<int>   mom;         /*!< D-1 momentum of the sink operator */
    std::string  mass;        /*!< Some string label for the mass(es) in the corr */
    std::string  ensemble;    /*!< Label for the ensemble */
  };

  //----------------------------------------------------------------------------
  //! Used for error output
  std::ostream& operator<<(std::ostream& os, const KeyHadron1PtCorr_t& d);
	
  //----------------------------------------------------------------------------
  //! KeyHadron1PtCorr reader
  void read(XMLReader& xml, const std::string& path, KeyHadron1PtCorr_t& param);

  //! KeyHadron1PtCorr writer
  void write(XMLWriter& xml, const std::string& path, const KeyHadron1PtCorr_t& param);

  //----------------------------------------------------------------------------
  //! KeyHadron1PtCorr reader
  void read(BinaryReader& bin, KeyHadron1PtCorr_t& param);

  //! KeyHadron1PtCorr writer
  void write(BinaryWriter& bin, const KeyHadron1PtCorr_t& param);


  //----------------------------------------------------------------------------
  //! Use a DB
  typedef FILEDB::ConfDataStoreDB< SerialDBKey<KeyHadron1PtCorr_t>, SerialDBData<VectorComplex> >  DBHadron1PtSingleCorr_t;

} // namespace ColorVec

#endif
