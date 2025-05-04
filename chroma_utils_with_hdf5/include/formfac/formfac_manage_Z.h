// -*- C++ -*-
// $Id: formfac_manage_Z.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 * \brief Manage amplitudes
 */

#ifndef __formfac_manage_Z_h__
#define __formfac_manage_Z_h__

#include "formfac/formfac_manage.h"
#include "formfac/formfac_manage_Z_key.h"
#include "formfac/formfac_manage_state.h"

namespace FF
{
  //----------------------------------------------------------------------------
  //! State function
  class StateAmpFunc : public StateFunc<AmpArg>
  {
  public:
    //! Help with cleanup
    virtual ~StateAmpFunc() {}

    //! Return 
    virtual std::string operator()(int cfg, const AmpArg& arg) const = 0;
  };


  //----------------------------------------------------------------------------
  //! Class to hold 2-pt functions
  class ManageAmpFunc : virtual public ManageEnsem<AmpArg, EnsemReal>
  {
  public:
    //! Virtual destructor
    virtual ~ManageAmpFunc() {}
  };


  //----------------------------------------------------------------------------
  //! Class to hold Z amplitudes
  class ManageAmpFuncMap : public ManageAmpFunc, public ManageEnsemMap<AmpArg, EnsemReal>
  {
  public:
    //! Virtual destructor
    virtual ~ManageAmpFuncMap() {}
  };


} // namespace FF

#endif
