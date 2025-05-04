// -*- C++ -*-
/*! \file
 * \brief Support for unordered maps
 */

#ifndef __map_traits_h__
#define __map_traits_h__

#include <unordered_map>
//#include "xml_array.h"
#include "io/adat_io.h"

namespace ADAT
{
  using namespace XMLArray;

  //! Hash function call
  size_t hashfunc(const std::string& that);

  //! Equivalence
  bool hashcmp(const std::string& that1, const std::string& that2);


  //! Support for unordered_map from the TR1 lib.
  template<typename K>
  struct UnorderedMapTraits
  {
    //! Hash function call
    size_t operator()(const K& that) const 
    {
       ADATIO::BinaryBufferWriter bin;
       write(bin, that);
       return hashfunc(bin.str());
    }

    //! Equivalence
    bool operator()(const K& that1, const K& that2) const
    {
      std::string a, b;
      {
	ADATIO::BinaryBufferWriter bin;
	write(bin, that1);
	a = bin.str();
      }
      {
	ADATIO::BinaryBufferWriter bin;
	write(bin, that2);
	b = bin.str();
      }
      return hashcmp(a,b);
    }
  };


} // namespace ADAT

#endif
