// -*- C++ -*-
// $Id: proginfo.h,v 1.1 2009/04/23 03:27:48 edwards Exp $
/*! \file
 *  \brief Print out basic info about this program
 */

#ifndef __proginfo_h__
#define __proginfo_h__

namespace ADAT 
{

  //! Print out basic info about this program
  /*!
   * Arguments:
   *
   *  \param xml          The xml stream to write the info
   */
  
  void proginfo(ADATXML::XMLWriter& xml);

}  // end namespace ADAT

#endif
