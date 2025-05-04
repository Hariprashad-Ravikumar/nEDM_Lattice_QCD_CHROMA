// -*- C++ -*-
// $Id: formfac_manage_npr.h,v 2.0 2008/12/05 04:43:36 edwards Exp $

/*! \file
 * \brief Manage NPR funcs
 */

#ifndef __formfac_manage_npr_h__
#define __formfac_manage_npr_h__

#include <map>
#include "ensem/ensem.h"
#include "adat/handle.h"
#include "io/adat_xml_group_reader.h"

#include "adat/singleton.h"
#include "adat/objfactory.h"

#include "formfac/formfac_manage.h"

namespace FF
{
  //----------------------------------------------------------------------------
  //! Three pt
  struct NprArg
  {
    NprArg() {links.resize(0); src.resize(1);src[0]=0; snk.resize(1);snk[0]=0; g=-1;}

    ArrayInt     mom;
    int          g;
    ArrayInt     links;
    std::string  src_name;
    ArrayInt     src;
    std::string  snk_name;
    ArrayInt     snk;
  };

  //! Support for maps
  bool operator<(const NprArg& a, const NprArg& b);


  //! State function
  class StateNprFunc
  {
  public:
    //! Help with cleanup
    virtual ~StateNprFunc() {}

    //! Return 
    virtual std::string operator()(int cfg, const NprArg& arg) const = 0;
  };

} // namespace FF

#endif
