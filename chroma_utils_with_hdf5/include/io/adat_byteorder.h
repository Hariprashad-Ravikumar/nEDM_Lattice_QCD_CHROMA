// -*- C++ -*-
// $Id: adat_byteorder.h,v 2.0 2008/12/05 04:43:37 edwards Exp $

/*! \file
 * \brief Byte manipulation
 */

#ifndef __adat_byteorder_h__
#define __adat_byteorder_h__

#include <cstdio>

namespace ADATUtil
{
#if 0
  // #include <adat_stdint.h>
#else
  // NEED PROPER CONSTRUCTION OF TYPES
  //! See adat_crc32.cc for defs and discussion
  typedef unsigned short int  n_uint16_t;
  typedef unsigned int        n_uint32_t;
#endif


  //! crc32
  n_uint32_t crc32(n_uint32_t crc, const unsigned char *buf, size_t len);

  //! crc32
  n_uint32_t crc32(n_uint32_t crc, const char *buf, size_t len);

  //! Is the native byte order big endian?
  bool big_endian();

  //! Byte-swap an array of data each of size nmemb
  void byte_swap(void *ptr, size_t size, size_t nmemb);

  //! fread on a binary file written in big-endian order
  size_t bfread(void *ptr, size_t size, size_t nmemb, FILE *stream);

  //! fwrite to a binary file in big-endian order
  size_t bfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
}

#endif
