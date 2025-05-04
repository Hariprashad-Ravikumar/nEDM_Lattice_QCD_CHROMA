// -*- C++ -*-
// $Id: ensem_random.h,v 2.0 2008/12/05 04:43:34 edwards Exp $
//
// ENSEM interface
//
// Random number support

namespace ENSEM {

//! Random number generator namespace
/*!
 * A collection of routines and data for supporting random numbers
 * 
 * It is a linear congruential with modulus m = 2**47, increment c = 0,
 * and multiplier a = (2**36)*m3 + (2**24)*m2 + (2**12)*m1 + m0.  
 */

namespace RNG
{
  //! Default initialization of the RNG
  /*! Uses arbitrary internal seed to initialize the RNG */
  bool initDefaultRNG(void);

  //! Initialize the internals of the RNG
  void initRNG(void);

  //! Hook to get code to link
  extern const bool registered;

  //! Initialize the RNG seed
  /*!
   * Seeds are big-ints
   */
  void setrn(const Seed& lseed);

  //! Recover the current RNG seed
  /*!
   * Seeds are big-ints
   */
  void savern(Seed& lseed);


  //! Internal RNG
  float sranf();
}


//! dest  = random
inline void
fill_random(float& d)
{
  d = float(RNG::sranf());
}

//! dest  = random
inline void
fill_random(double& d)
{
  d = double(RNG::sranf());
}

} // namespace ENSEM
