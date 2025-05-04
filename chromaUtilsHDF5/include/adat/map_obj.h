// -*- C++ -*-
/*! \file
 * \brief Wrapper over maps
 */

#ifndef __map_obj_h__
#define __map_obj_h__

#include <unordered_map>
#include <vector>
#include <string>
#include <map>

#include "io/adat_io.h"

namespace ADAT
{

  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObject
  {
  public:
    //! Default constructor
    MapObject() {}

    //! Destructor
    ~MapObject() {}

    //! Exists?
    bool exist(const K& key) const {
      ADATIO::BinaryBufferWriter bin;
      write(bin, key);
      return (src_map.find(bin.str()) == src_map.end()) ? false : true;
    }
			
    //! Clear the object
    void clear() {src_map.clear();}

    //! Erase a key-value
    void erase(const K& key) {
      ADATIO::BinaryBufferWriter bin;
      write(bin, key);

      typename MapType_t::const_iterator iter = src_map.find(bin.str());
      if (iter != src_map.end())
      {
	src_map.erase(iter);
      }
    }
			
    //! Insert
    void insert(const K& key, const V& val) {
      ADATIO::BinaryBufferWriter bin;
      write(bin, key);

      const std::string bin_key(bin.str());
      typename MapType_t::iterator iter = src_map.find(bin_key);
      if (iter != src_map.end())
      {
	iter->second = val;
      }
      else
      {
	src_map.insert(std::make_pair(bin_key,val));
      }
    }
			
    //! Getter
    const V& operator[](const K& key) const {
      ADATIO::BinaryBufferWriter bin;
      write(bin, key);

      typename MapType_t::const_iterator iter = src_map.find(bin.str());
      if (iter == src_map.end())
      {
	std::cerr << "MapObject: key not found" << std::endl;
	exit(1);
      }

      return iter->second;
    }
			
    //! Setter
    V& operator[](const K& key) {
      ADATIO::BinaryBufferWriter bin;
      write(bin, key);

      typename MapType_t::iterator iter = src_map.find(bin.str());
      if (iter == src_map.end())
      {
	std::cerr << "MapObject: key not found" << std::endl;
	exit(1);
      }

      return iter->second;
    }
			
    //! The number of elements
    size_t size() const {return src_map.size();}

    //! Dump keys
    std::vector<K> keys() const {
      std::vector<K> _keys;
      typename MapType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter)
      {
	ADATIO::BinaryBufferReader bin(iter->first);
	K key;
	read(bin, key);
	_keys.push_back(key);
      }
      return _keys;
    }

    //! Dump keys and values
    virtual void keysAndValues(std::vector<K>& _keys, std::vector<V>& _vals) const {
      _keys.resize(0);
      _vals.resize(0);
      typename MapType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter)
      {
	ADATIO::BinaryBufferReader bin(iter->first);
	K key;
	read(bin, key);
	_keys.push_back(key);
	_vals.push_back(iter->second);
      }
    }

  protected:  
    //! Map type convenience
    typedef std::unordered_map<std::string, V> MapType_t;
    
    //! Map of objects
    mutable MapType_t src_map;
  };

} // namespace ColorVec

#endif
