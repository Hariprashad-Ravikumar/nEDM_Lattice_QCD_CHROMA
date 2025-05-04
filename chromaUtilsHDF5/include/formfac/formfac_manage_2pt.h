// -*- C++ -*-
// $Id: formfac_manage_2pt.h,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Manage 2-pt funcs
 */

#ifndef __formfac_manage_2pt_h__
#define __formfac_manage_2pt_h__

#include "formfac/formfac_manage_2pt_key.h"
#include "formfac/formfac_manage_state.h"
#include "formfac/formfac_manage.h"

namespace FF
{
  //----------------------------------------------------------------------------
  //! State function
  class State2PtFunc : public StateFunc<TwoPtArg>
  {
  public:
    //! Help with cleanup
    virtual ~State2PtFunc() {}

    //! Return 
    virtual std::string operator()(int cfg, const TwoPtArg& arg) const = 0;
  };


  //----------------------------------------------------------------------------
  //! Class to hold 2-pt functions
  class Manage2PtFunc : virtual public ManageEnsem<TwoPtArg, EnsemVectorReal>
  {
  public:
    //! Virtual destructor
    virtual ~Manage2PtFunc() {}
  };


  //----------------------------------------------------------------------------
  //! Class to hold 2-pt functions
  class Manage2PtFuncMap : public Manage2PtFunc, public ManageEnsemMap<TwoPtArg, EnsemVectorReal>
  {
  public:
    //! Destructor
    virtual ~Manage2PtFuncMap() {}
  };

} // namespace FF

#endif
