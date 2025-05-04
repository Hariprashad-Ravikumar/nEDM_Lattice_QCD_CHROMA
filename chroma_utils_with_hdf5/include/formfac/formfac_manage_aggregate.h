// -*- C++ -*-
// $Id: formfac_manage_aggregate.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 *  \brief All formfactor factories
 */

#ifndef __formfac_manage_aggregate_h__
#define __formfac_manage_aggregate_h__

namespace FF
{
  //! Registration aggregator
  namespace FormfacManage2PtFuncEnv
  {
    bool registerAll();
  }

  //! Registration aggregator
  namespace FormfacManageEnergyFuncEnv
  {
    bool registerAll();
  }

  //! Registration aggregator
  namespace FormfacManageAmpFuncEnv
  {
    bool registerAll();
  }

  //! Registration aggregator
  namespace FormfacManage3PtFuncEnv
  {
    bool registerAll();
  }
}

#endif
