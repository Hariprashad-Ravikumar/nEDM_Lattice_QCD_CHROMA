// -*- C++ -*-
// $Id: formfac_bb.h,v 2.0 2008/12/05 04:43:35 edwards Exp $

/*! \file
 * \brief Utilities
 */

#ifndef __formfac_bb_h__
#define __formfac_bb_h__

#include "ensem/ensem.h"
#include <vector>

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  typedef Array<int>  ArrayInt;   // save some typing
  typedef Array<double> ArrayDouble;  // save some typing

  //! Image of building blocks file
  struct BuildingBlocks_t
  {
    struct Param_t
    {
      std::string        SeqSourceType;
      signed short int   Flavor;  // currently assumes u and d are given as f=0 and f=1
      unsigned short int Contraction;
      signed short int   GammaInsertion;
      unsigned short int DecayDir;
      unsigned short int NMomPerms;
      
      Array<int>   sink_mom;
      int          links_max;
      Array<int>   latt_size;
      int          t_srce;
      int          t_sink;
    };

    struct Links_t
    {
      struct Gamma_t
      {
	struct Momenta_t
	{
	  Array<int>       inser_mom;
	  Array<ComplexF>  current;
	};

	int               gamma_value;
	Array<Momenta_t>  momenta;
      };

      Array<int>        link_value;
      Array<Gamma_t>    gamma;
    };

    int                  version;
    Param_t              param;
    Array<Links_t>       links;
  };


  //! Mother of all readers
  void read(const std::string& BuildingBlocksFileName, BuildingBlocks_t& bar);

  //! Mother of all readers
  void read(const std::string& BuildingBlocksFileName, std::vector<BuildingBlocks_t> &bar);

  //! Writer
  void write(const std::string& BuildingBlocksFileName, const BuildingBlocks_t& bar);

  //! Thin out a BB structure
  void thinBB(BuildingBlocks_t& new_bar, 
	      int new_links_max,
	      const BuildingBlocks_t& old_bar);

} // namespace FF

#endif
