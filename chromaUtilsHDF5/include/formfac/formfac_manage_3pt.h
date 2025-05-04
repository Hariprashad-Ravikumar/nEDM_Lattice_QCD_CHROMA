// -*- C++ -*-
// $Id: formfac_manage_3pt.h,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Manage 3-pt funcs
 */

#ifndef __formfac_manage_3pt_h__
#define __formfac_manage_3pt_h__

#include "adat/handle.h"
#include "ensem/ensem.h"
#include "formfac/formfac_manage.h"
#include "formfac/formfac_ensemble.h"
#include "formfac/formfac_manage_3pt_key.h"
#include "formfac/formfac_manage_state.h"

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  //----------------------------------------------------------------------------
  //! State function
  class State3PtFunc : public StateFunc<ThreePtArg>
  {
  public:
    //! Help with cleanup
    virtual ~State3PtFunc() {}

    //! Return 
    virtual std::string operator()(int cfg, const ThreePtArg& arg) const = 0;
  };


  //----------------------------------------------------------------------------
  //! Interface for 3-pt functions using a reduced size key
  class Manage3PtFuncReduced;


  //----------------------------------------------------------------------------
  //! Abstract interface for 3-pt functions
  class Manage3PtFunc : virtual public ManageEnsem<ThreePtArg, EnsemVectorComplexF>
  {
  public:
    //! Destructor
    virtual ~Manage3PtFunc() {}

    //! Return a "view" to a managed object using a prototype key "arg"
    virtual Manage3PtFuncReduced* createView(const ThreePtArg& arg);
  };


  //----------------------------------------------------------------------------
  //! Class to hold 3-pt functions with map semantics
  class Manage3PtFuncMap : public Manage3PtFunc, public ManageEnsemMap<ThreePtArg, EnsemVectorComplexF>
  {
  public:
    //! Destructor
    virtual ~Manage3PtFuncMap() {}
  };


  //----------------------------------------------------------------------------
  //! Interface for 3-pt functions using a reduced size key. Could be abstract if needed.
  class Manage3PtFuncReduced : virtual public ManageEnsem<ThreePtArgReduced, EnsemVectorComplexF>
  {
  public:
    //! Constructor
    Manage3PtFuncReduced(ADAT::Handle<Manage3PtFunc> obj, const ThreePtArg& proto_key);

    //! Destructor
    ~Manage3PtFuncReduced();

    //! Return the prototype key
    ThreePtArg getProtoKey() const;

    //! Return a value given a key
    EnsemVectorComplexF operator[](const ThreePtArgReduced& arg);

    //! Number of bins
    int size() const;

    //! Time extent
    int timeLen() const;

    //! Decay direction
    int decayDir() const;

    //! Lattice size
    const ArrayInt& lattSize() const;

  private:
    ADAT::Handle<Manage3PtFunc>  obj;
    ThreePtArg                   proto_key;
  };


} // namespace FF

#endif
