// -*- C++ -*-
// $Id: formfac_qsq.h,v 2.0 2008/12/05 04:43:36 edwards Exp $

/*! \file
 * \brief Construct Qsq list
 */

#ifndef __formfac_qsq_h__
#define __formfac_qsq_h__

#include "formfac/formfac_ensemble.h"
#include <list>

namespace FF
{
  //! 4D dot product
  EnsemReal dot(const EnsemReal& E1, const ArrayDouble& p1, const EnsemReal& E2, const ArrayDouble& p2);

  //! 4D Minkowski q^2
  EnsemReal minkQsq(const EnsemReal& E1, const ArrayDouble& p1, const EnsemReal& E2, const ArrayDouble& p2);

  //! Structure holding all info at a single Qsq point
  struct QsqVal
  {
    double           Qsq_c;
    double           Qsq_l;
    int              qsq;
    std::list<PiPf>  pi_pf;
  };

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const QsqVal& val);


  //! Info need for Qsq map construction
  struct QsqParam
  {
    int                 nq_norm_sq_max;
    int                 ni_norm_sq_max;
    LatticeParam        lattice;
    double              mass_i;
    double              mass_f;
    std::list<ArrayInt> nf_list;
  };

  //! 4D Minkowski q^2 using continuum disp. relation in physical units
  double qsqContDispPhys(double m1, const ArrayInt& p1, double m2, const ArrayInt& p2,
			 const LatticeParam& lattice);

  //! 4D Minkowski q^2 using lattice disp. relation in physical units
  double qsqLattDispPhys(double m1, const ArrayInt& p1, double m2, const ArrayInt& p2,
			 const LatticeParam& lattice);

  //! Make a continuum momentum from integers
  ArrayDouble contMom(const ArrayInt& nf, int ns);

  //! Make a lattice momentum from integers
  ArrayDouble lattMom(const ArrayInt& nf, int ns);

  //! Canonically order an array of momenta
  /*! \return abs(mom[0]) >= abs(mom[1]) >= ... >= abs(mom[mu]) >= ... >= 0 */
  ArrayInt canonicalOrder(const ArrayInt& mom);

  //! Decompose a lexicographic site into coordinates
  Array<int> crtesn(int ipos, const Array<int>& latt_size);

  //! List pretty printer
  std::ostream& operator<<(std::ostream& s, const std::list<QsqVal>& qsq_val);

  //! Construct the map of all Qsq points
  std::list<QsqVal> constructQsq(const QsqParam& param);

  //! Construct the map of all Qsq points
  /*! Fixes norm2(p_f), norm2(p_i) and norm2(q) */
  std::list<QsqVal> constructQsqFixedMom(const QsqParam& param);

} // namespace FF

#endif
