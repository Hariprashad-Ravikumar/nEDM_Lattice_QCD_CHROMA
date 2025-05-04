// -*- C++ -*-
// $Id: formfac_manage_2pt_stripped.h,v 2.1 2009/03/26 21:21:44 edwards Exp $
/*! \file
 * \brief Manage 2-pt funcs using stripped ensemble files
 */

#ifndef __formfac_manage_2pt_stripped_h__
#define __formfac_manage_2pt_stripped_h__

#include "adat/map_traits.h"
#include "adat/handle.h"
#include "formfac/formfac_manage_2pt.h"
#include "io/adat_xml_group_reader.h"

namespace FF
{
  //! Manager factory
  namespace FormfacManage2PtFuncStrippedEnv
  { 
    // Register the callbacks
    bool registerAll(void);
  }


  //----------------------------------------------------------------------------
  //! Parameters
  struct Manage2PtFuncStrippedParams_t
  {
    Manage2PtFuncStrippedParams_t(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_in, const std::string& path) const;

    GroupXML_t     state;            /*!< xml holding state group */
    bool           avg_2pt_func;
    LatticeParam   lattice;
  };

  // Parameters
  void read(XMLReader& xml_in, const std::string& path, Manage2PtFuncStrippedParams_t& param);

  // Parameters
  void write(XMLWriter& xml_out, const std::string& path, const Manage2PtFuncStrippedParams_t& param);


  //----------------------------------------------------------------------------
  //! Class to hold 2-pt functions
  class Manage2PtFuncStripped : public Manage2PtFuncMap
  {
  public:
    //! Empty Constructor
    Manage2PtFuncStripped() {}

    //! Constructor
    Manage2PtFuncStripped(ADAT::Handle<State2PtFunc> state_, 
			  bool avg_2pt_func_,
			  const LatticeParam& lattice);

    //! destructor to help with cleanup;
    ~Manage2PtFuncStripped() {}

    //! Virtual constructor
    void create(ADAT::Handle<State2PtFunc> state_, 
		bool avg_2pt_func_,
		const LatticeParam& lattice);

    //! Return a 2-pt ref. Fetch if not in memory
    EnsemVectorReal operator[](const TwoPtArg& arg);

    //! Always insert a key,value
//    void insert(const TwoPtArg& k, const EnsemVectorReal& v);

    //! Does key exist?
    bool exist(const TwoPtArg& k);

    //! Number of bins
    int size() const {return ensemble_info.nbin;}

    //! Time extent
    int timeLen() const {return ensemble_info.lattice.latt_size[ensemble_info.lattice.decay_dir];}

    //! Decay direction
    int decayDir() const {return ensemble_info.lattice.decay_dir;}

    //! Lattice Size
    const ArrayInt& lattSize() const {return ensemble_info.lattice.latt_size;}

  protected:
    //! Read files
    void do_read2pt(const TwoPtArg& arg);

  private:
    ADAT::Handle<State2PtFunc>  state;
    bool                  avg_2pt_func;
    EnsembleInfo          ensemble_info;

    //! TR1 unordered_map version (a hash table)
    typedef std::unordered_map<TwoPtArg,
			       EnsemVectorReal,
			       ADAT::UnorderedMapTraits<TwoPtArg>,
			       ADAT::UnorderedMapTraits<TwoPtArg> > TwoPtType_t;

    //! Memory version of a data-base
    TwoPtType_t         twopt;
  };


} // namespace FF

#endif
