// -*- C++ -*-
// $Id: formfac_manage_factory.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 *  \brief Factory for form-factor managers
 */

#ifndef __formfac_manage_factory_h__
#define __formfac_manage_factory_h__

#include "formfac/formfac_manage_2pt.h"
#include "formfac/formfac_manage_E.h"
#include "formfac/formfac_manage_Z.h"
#include "formfac/formfac_manage_3pt.h"

#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "adat/handle.h"

namespace FF
{
  //----------------------------------------------------------------------------
  using namespace Util;

  //! 2-pt state factory
  typedef SingletonHolder< 
    ObjectFactory<State2PtFunc,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  State2PtFunc* (*)(XMLReader&,
				    const std::string&),
		  StringFactoryError> >
  TheState2PtFuncFactory;


  //! Energy state factory
  typedef SingletonHolder< 
    ObjectFactory<StateEnergyFunc,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  StateEnergyFunc* (*)(XMLReader&,
				       const std::string&),
		  StringFactoryError> >
  TheStateEnergyFuncFactory;


  //! Amp state factory
  typedef SingletonHolder< 
    ObjectFactory<StateAmpFunc,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  StateAmpFunc* (*)(XMLReader&,
				    const std::string&),
		  StringFactoryError> >
  TheStateAmpFuncFactory;


  //! 3-pt state factory
  typedef SingletonHolder< 
    ObjectFactory<State3PtFunc,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  State3PtFunc* (*)(XMLReader&,
				    const std::string&),
		  StringFactoryError> >
  TheState3PtFuncFactory;


  //! 2-pt management factory
  typedef SingletonHolder< 
    ObjectFactory<Manage2PtFunc,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  Manage2PtFunc* (*)(XMLReader&,
				     const std::string&), 
		  StringFactoryError> >
  TheManage2PtFuncFactory;

  //! 2-pt management factory
  typedef SingletonHolder< 
    ObjectFactory<Manage2PtFuncMap,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  Manage2PtFuncMap* (*)(XMLReader&,
					const std::string&), 
		  StringFactoryError> >
  TheManage2PtFuncMapFactory;


  //! Energy management factory
  typedef SingletonHolder< 
    ObjectFactory<ManageEnergyFunc,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  ManageEnergyFunc* (*)(XMLReader&,
					const std::string&), 
		  StringFactoryError> >
  TheManageEnergyFuncFactory;

  //! Energy management factory
  typedef SingletonHolder< 
    ObjectFactory<ManageEnergyFuncMap,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  ManageEnergyFuncMap* (*)(XMLReader&,
					   const std::string&), 
		  StringFactoryError> >
  TheManageEnergyFuncMapFactory;


  //! Amp management factory
  typedef SingletonHolder< 
    ObjectFactory<ManageAmpFunc,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  ManageAmpFunc* (*)(XMLReader&,
				     const std::string&), 
		  StringFactoryError> >
  TheManageAmpFuncFactory;

  //! 3-pt management factory
  typedef SingletonHolder< 
    ObjectFactory<ManageAmpFuncMap,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  ManageAmpFuncMap* (*)(XMLReader&,
					const std::string&), 
		  StringFactoryError> >
  TheManageAmpFuncMapFactory;


  //! 3-pt management factory
  typedef SingletonHolder< 
    ObjectFactory<Manage3PtFunc,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  Manage3PtFunc* (*)(XMLReader&,
				     const std::string&), 
		  StringFactoryError> >
  TheManage3PtFuncFactory;

  //! 3-pt management factory
  typedef SingletonHolder< 
    ObjectFactory<Manage3PtFuncMap,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  Manage3PtFuncMap* (*)(XMLReader&,
					const std::string&), 
		  StringFactoryError> >
  TheManage3PtFuncMapFactory;


  //! Reduced 3-pt management factory
  typedef SingletonHolder< 
    ObjectFactory<Manage3PtFuncReduced,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  Manage3PtFuncReduced* (*)(XMLReader&,
					    const std::string&), 
		  StringFactoryError> >
  TheManage3PtFuncReducedFactory;

}


#endif
