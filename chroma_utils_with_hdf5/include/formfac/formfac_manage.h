// -*- C++ -*-
// $Id: formfac_manage.h,v 2.0 2008/12/05 04:43:35 edwards Exp $

/*! \file
 * \brief Manage ensembles
 */

#ifndef __formfac_manage_h__
#define __formfac_manage_h__

#include "formfac/formfac_ensemble.h"

namespace FF
{
  //----------------------------------------------------------------------------
  //! Abstract class of managers
  template<typename K, typename V>
  class ManageEnsem
  {
  public:
    //! Destructor
    virtual ~ManageEnsem() {}

    //! Return a value given a key
    virtual V operator[](const K& arg) = 0;

    //! Number of bins
    virtual int size() const = 0;

    //! Time extent
    virtual int timeLen() const = 0;

    //! Decay direction
    virtual int decayDir() const = 0;

    //! Lattice size
    virtual const ArrayInt& lattSize() const = 0;
  };


  //----------------------------------------------------------------------------
  //! Abstract class of managers that resembles a map
  template<typename K, typename V>
  class ManageEnsemMap : virtual public ManageEnsem<K,V>
  {
  public:
    //! Destructor
    virtual ~ManageEnsemMap() {}

    //! Does key exist?
    virtual bool exist(const K& k) = 0;
  };

} // namespace FF

#endif
