// -*- C++ -*-
// $Id: formfac_state.h,v 2.0 2008/12/05 04:43:37 edwards Exp $

/*! \file
 * \brief State management functions
 */

#ifndef __formfac_state_h__
#define __formfac_state_h__

#warning "THIS FILE NEEDS TO BE REARRANGED"
#if 0


#include "list"
#include "ensem/ensem.h"
#include "formfac/formfac_manage.h"
#include "formfac/formfac_qsq.h"

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  //! Driver
//  void ff_driver(XMLWriter& xml_out, LLSqComponent<EnsemReal>& sys, const std::list<QsqVal>& qsq_val);

  //-----------------------------------------------------------------------------------

  //! Make a state
  class Local3PtFunc : public State3PtFunc
  {
  public:
    //! Constructor
    Local3PtFunc(const std::string& pattern_) : pattern(pattern_) {}

    //! Construct 3pt names
    std::string operator()(int cfg, const ThreePtArg& arg) const
      {
	ArrayInt p_i = arg.pi_pf.p_i;
	ArrayInt p_f = arg.pi_pf.p_f;
	ArrayInt q   = p_f - p_i;      // note sign convention
  
	// cout << __func__ << ": pattern=" << pattern.c_str() << std::endl;

	char qxsign = (q[0] < 0 ) ? '-' : '+';
	char qysign = (q[1] < 0 ) ? '-' : '+';
	char qzsign = (q[2] < 0 ) ? '-' : '+';

	char pfxsign = (p_f[0] < 0 ) ? '-' : '+';
	char pfysign = (p_f[1] < 0 ) ? '-' : '+';
	char pfzsign = (p_f[2] < 0 ) ? '-' : '+';

	char filen[1024];
	sprintf(filen, pattern.c_str(), 
		pfzsign, abs(p_f[2]), pfysign, abs(p_f[1]), pfxsign, abs(p_f[0]),
		cfg,
		qzsign, abs(q[2]), qysign, abs(q[1]), qxsign, abs(q[0]));

	return std::string(filen);
      }

  private:
    std::string pattern;
  };

  
  //! State function
  class Local2PtFunc : public State2PtFunc
  {
  public:
    //! Constructor
    Local2PtFunc(const std::string& pattern_, const std::string& particle_) :
      pattern(pattern_), particle(particle_) {}

    //! Construct 2pt names
    std::string operator()(const TwoPtArg& arg) const
      {
	ArrayInt zero(Nd-1);
	zero = 0;
	std::ostringstream pion;

	pion << particle;

	if (arg.p != zero)
	  pion << "_px" << arg.p[0] << "_py" << arg.p[1] << "_pz" << arg.p[2];

	char filen[1024];
	sprintf(filen, pattern.c_str(), pion.str().c_str());

	return std::string(filen);
      }

  private:
    std::string pattern;
    std::string particle;
  };


  //! State function
  class LocalZFunc : public State2PtFunc
  {
  public:
    //! Constructor
    LocalZFunc(const std::string& pattern_, const std::string& particle_) :
      pattern(pattern_), particle(particle_) {}

    //! Construct Z names
    std::string operator()(const TwoPtArg& arg) const
      {
	ArrayInt zero(Nd-1);
	zero = 0;
	std::ostringstream pion;

	pion << particle;

	if (arg.p != zero)
	  pion << "_px" << arg.p[0] << "_py" << arg.p[1] << "_pz" << arg.p[2];

	char filen[1024];
	sprintf(filen, pattern.c_str(), pion.str().c_str());

	return std::string(filen);
      }

  private:
    std::string pattern;
    std::string particle;
  };


  //! State function
  class LocalEnergyFunc : public StateEnergyFunc
  {
  public:
    //! Constructor
    LocalEnergyFunc(const std::string& pattern_, const std::string& particle_) : 
      pattern(pattern_), particle(particle_) {}

    // Construct energy names
    std::string operator()(const ArrayInt& mom) const
      {
	ArrayInt zero(Nd-1);
	zero = 0;
	std::ostringstream pion;

	pion << particle;

	if (mom != zero)
	  pion << "_px" << mom[0] << "_py" << mom[1] << "_pz" << mom[2];

	char filen[1024];
	sprintf(filen, pattern.c_str(), pion.str().c_str());

	return std::string(filen);
      }

  private:
    std::string pattern;
    std::string particle;
  };


} // namespace FF

#endif

#endif
