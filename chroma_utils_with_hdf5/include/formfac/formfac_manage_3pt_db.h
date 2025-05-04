// -*- C++ -*-
// $Id: formfac_manage_3pt_db.h,v 2.3 2009/03/09 02:41:17 edwards Exp $

/*! \file
 * \brief Manage 3-pt funcs
 */

#ifndef __formfac_manage_3pt_db_h__
#define __formfac_manage_3pt_db_h__

#include "formfac/formfac_manage_3pt.h"
#include "adat/handle.h"
#include "io/key_val_db.h"
#include "formfac/hadron_3pt_corr.h"

#include "AllConfStoreDB.h"


namespace FF
{
  //! Manager factory
  namespace FormfacManage3PtFuncDBEnv
  { 
    // Register the callbacks
    bool registerAll(void);
  }


  //----------------------------------------------------------------------------
  //! BB parameters
  struct Manage3PtFuncDBParams_t
  {
    Manage3PtFuncDBParams_t(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_in, const std::string& path) const;

    std::string          db_file;
    KeyHadron3PtCorr_t   proto_key;
    LatticeParam         lattice;
  };

  // Parameters
  void read(XMLReader& xml_in, const std::string& path, Manage3PtFuncDBParams_t& param);

  // Parameters
  void write(XMLWriter& xml_out, const std::string& path, const Manage3PtFuncDBParams_t& param);


  //----------------------------------------------------------------------------
  //! Class to hold 3-pt functions specific to use database
  class Manage3PtFuncDB : public Manage3PtFunc
  {
  public:
    //! Constructor
    Manage3PtFuncDB(const std::string& db_file_,
		    const KeyHadron3PtCorr_t& proto_key_,
		    const LatticeParam& lattice);

    //! Destructor
    ~Manage3PtFuncDB();

    //! Return a value given a key
    EnsemVectorComplexF operator[](const ThreePtArg& arg);

    //! Number of bins
    int size() const;

    //! Time extent
    int timeLen() const;

    //! Decay direction
    int decayDir() const;

    //! Lattice size
    const ArrayInt& lattSize() const;

  private:
    // Save some typing
    typedef KeyHadron3PtCorr_t       K;
    typedef EnsemVectorComplexF      V;
    typedef EnsemScalar<V>::Type_t   SV;

    FILEDB::AllConfStoreDB< SerialDBKey<K>,  SerialDBData<SV> > database;
    K               proto_key;
    LatticeParam    lattice;
  };


} // namespace FF

#endif
