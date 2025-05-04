// -*- C++ -*-
// $Id: formfac_manage_E_stripped.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 * \brief Manage energy funcs
 */

#ifndef __formfac_manage_E_stripped_h__
#define __formfac_manage_E_stripped_h__

#include "adat/map_traits.h"
#include "formfac/formfac_manage_E.h"
#include "adat/handle.h"
#include "io/adat_xml_group_reader.h"

namespace FF
{
  //! Manager factory
  namespace FormfacManageEnergyFuncStrippedEnv
  { 
    // Register the callbacks
    bool registerAll(void);
  }


  //----------------------------------------------------------------------------
  //! Parameters
  struct ManageEnergyFuncStrippedParams_t
  {
    ManageEnergyFuncStrippedParams_t(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_in, const std::string& path) const;

    GroupXML_t     state;            /*!< xml holding state group */
    LatticeParam   lattice;
  };

  // Parameters
  void read(XMLReader& xml_in, const std::string& path, ManageEnergyFuncStrippedParams_t& param);

  // Parameters
  void write(XMLWriter& xml_out, const std::string& path, const ManageEnergyFuncStrippedParams_t& param);


  //----------------------------------------------------------------------------
  //! Class to hold energies
  class ManageEnergyFuncStripped : public ManageEnergyFuncMap
  {
  public:
    //! Empty Constructor
    ManageEnergyFuncStripped() {}

    //! Constructor
    ManageEnergyFuncStripped(ADAT::Handle<StateEnergyFunc> state_,
			     const LatticeParam& lattice);

    //! destructor to help with cleanup;
    ~ManageEnergyFuncStripped() {}

    //! Virtual constructor
    void create(ADAT::Handle<StateEnergyFunc> state_,
		const LatticeParam& lattice);

    //! Return an energy
    EnsemReal operator[](const EnergyArg& p);

    //! Insert a key,value
    void insert(const EnergyArg& k, const EnsemReal& v);

    //! Does key exist?
    bool exist(const EnergyArg& k);

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
    void do_readE(const EnergyArg& p);

    //! Compute via dispersion relation
    void do_computeE(const EnergyArg& p);

  private:
    ADAT::Handle<StateEnergyFunc>  state;
    EnsembleInfo          ensemble_info;

    //! TR1 unordered_map version (a hash table)
    typedef std::unordered_map<EnergyArg,
			       EnsemReal,
			       ADAT::UnorderedMapTraits<EnergyArg>,
			       ADAT::UnorderedMapTraits<EnergyArg> > EnergyType_t;

    //! memory version of a data-base
    EnergyType_t         E;
  };

} // namespace FF

#endif
