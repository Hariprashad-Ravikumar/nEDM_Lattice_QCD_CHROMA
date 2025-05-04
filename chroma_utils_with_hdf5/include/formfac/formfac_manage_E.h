// -*- C++ -*-
// $Id: formfac_manage_E.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 * \brief Manage energy funcs
 */

#ifndef __formfac_manage_E_h__
#define __formfac_manage_E_h__

#include "formfac/formfac_manage.h"
#include "formfac/formfac_manage_E_key.h"
#include "formfac/formfac_manage_state.h"

namespace FF
{
  //----------------------------------------------------------------------------
  //! State function
  class StateEnergyFunc : public StateFunc<EnergyArg>
  {
  public:
    //! Help with cleanup
    virtual ~StateEnergyFunc() {}

    //! Return 
    virtual std::string operator()(int cfg, const EnergyArg& arg) const = 0;
  };


  //----------------------------------------------------------------------------
  //! Class to hold energies
  class ManageEnergyFunc : virtual public ManageEnsem<EnergyArg, EnsemReal>
  {
  public:
    //! Virtual destructor
    virtual ~ManageEnergyFunc() {}
  };


  //----------------------------------------------------------------------------
  //! Class to hold energies
  class ManageEnergyFuncMap : public ManageEnergyFunc, public ManageEnsemMap<EnergyArg, EnsemReal>
  {
  public:
    //! Virtual destructor
    virtual ~ManageEnergyFuncMap() {}
  };

} // namespace FF

#endif
