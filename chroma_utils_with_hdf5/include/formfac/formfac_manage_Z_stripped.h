// -*- C++ -*-
// $Id: formfac_manage_Z_stripped.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 * \brief Manage amplitudes
 */

#ifndef __formfac_manage_Z_stripped_h__
#define __formfac_manage_Z_stripped_h__

#include "adat/map_traits.h"
#include "formfac/formfac_manage_Z.h"
#include "adat/handle.h"
#include "io/adat_xml_group_reader.h"

namespace FF
{
  //! Manager factory
  namespace FormfacManageAmpFuncStrippedEnv
  { 
    // Register the callbacks
    bool registerAll(void);
  }


  //----------------------------------------------------------------------------
  //! Parameters
  struct ManageAmpFuncStrippedParams_t
  {
    ManageAmpFuncStrippedParams_t(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_in, const std::string& path) const;

    GroupXML_t     state;            /*!< xml holding state group */
    bool           avg_2pt_func;
    LatticeParam   lattice;
  };

  // Parameters
  void read(XMLReader& xml_in, const std::string& path, ManageAmpFuncStrippedParams_t& param);

  // Parameters
  void write(XMLWriter& xml_out, const std::string& path, const ManageAmpFuncStrippedParams_t& param);


  //----------------------------------------------------------------------------
  //! Class to hold Z amplitudes
  class ManageAmpFuncStripped : public ManageAmpFuncMap
  {
  public:
    //! Empty constructor
    ManageAmpFuncStripped() {}

    //! Constructor
    ManageAmpFuncStripped(ADAT::Handle<StateAmpFunc> state_, 
			  bool avg_2pt_func_,
			  const LatticeParam& lattice);

    //! destructor to help with cleanup;
    ~ManageAmpFuncStripped() {}

    //! Virtual constructor
    void create(ADAT::Handle<StateAmpFunc> state_, 
		bool avg_2pt_func_,
		const LatticeParam& lattice);

    //! Return a Z ref. Fetch if not in memory
    EnsemReal operator[](const AmpArg& arg);

    //! Always insert a key,value
    void insert(const AmpArg& k, const EnsemReal& v);

    //! Does key exist?
    bool exist(const AmpArg& k);

    //! Number of bins
    int size() const {return ensemble_info.nbin;}

    //! Time extent
    int timeLen() const {return 1;}

    //! Decay direction
    int decayDir() const {return ensemble_info.lattice.decay_dir;}

    //! Lattice Size
    const ArrayInt& lattSize() const {return ensemble_info.lattice.latt_size;}

  protected:
    //! Read files
    void do_readZ(const AmpArg& arg);

  private:
    ADAT::Handle<StateAmpFunc>  state;
    bool                  avg_2pt_func;
    EnsembleInfo          ensemble_info;

    //! TR1 unordered_map version (a hash table)
    typedef std::unordered_map<AmpArg,
			       EnsemReal,
			       ADAT::UnorderedMapTraits<AmpArg>,
			       ADAT::UnorderedMapTraits<AmpArg> > AmpType_t;

    //! Memory version of a data-base
    AmpType_t         Z;
  };


} // namespace FF

#endif
