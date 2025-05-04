// -*- C++ -*-
// $Id: formfac_solver.h,v 2.0 2008/12/05 04:43:37 edwards Exp $

/*! \file
 * \brief Solver interface
 */

#ifndef __formfac_solver_h__
#define __formfac_solver_h__

#include "formfac/formfac_solver_row.h"
#include "formfac/formfac_qsq.h"
#include <string>
#include <vector>

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  //-----------------------------------------------------------------------------------
  //! Generate components for linear least squares
  /*! @ingroup formfac
   *
   *  Construct components for SVD solver
   */
  template<typename T>
  class LLSqComponent
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~LLSqComponent() {}

    //! checks if an entire row is zero and if so ignores the row and the rhs
    virtual bool isRowZero(const LLSqRow_t& row) const
      {
	bool isz = true;
	std::vector<std::string> FF_names(unknowns());
	for(int nFF = 0; nFF < FF_names.size(); nFF++)
	{
	  isz &= isZero(mat(row,FF_names[nFF]));
	  if (! isz)
	    break;
	}
	return isz;
      }

    //! The names of the unknowns. These are the possible column indices
    virtual std::vector<std::string> unknowns() const = 0;

    //! Number of FF
    virtual int nUnknowns() const {return unknowns().size();}

    //! Set qsq to use
    virtual void setQsq(const QsqVal& qsq) = 0;

    //! Reset time fit interval
    virtual void setTimeFit(const TimeFitRange_t& range) = 0;

    //! Get time fit interval
    virtual TimeFitRange_t getTimeFit() const = 0;

    //! Get time fit interval
    virtual std::vector<TimeFitRange_t> getFitRanges() const = 0;

    //! Time extent
    virtual int timeLen() const = 0;

    //! Return the list of the row indices
    virtual std::vector<LLSqRow_t> rowList() const = 0;

    //! Given row and column indices, return the corresponding matrix element
    virtual T mat(const LLSqRow_t& row, const std::string& FF_name) const = 0;

    //! Given row index, return the corresponding rhs element
    virtual T rhs(const LLSqRow_t& row) const = 0;

    //! Number of bins
    virtual int nbins() const = 0;

    //! Possibly save the solutions
    virtual void save(const Array<EnsemReal>& F, 
		      const Real& qsq, int t_i, int t_f) const
      {
	// default behavior is to do nothing
      }
  };


  //-----------------------------------------------------------------------------------
  //! Hold results of linear least squares
  template<typename T>
  struct LLSqResults
  {
    LLSqResults() {count_singular_val = ndof = 0;}

    int  ndof;
    int  count_singular_val;
    
    Array<T>  x;      /*!< Solutions */

    Array2d<T> x_cov;

    Real      Q;
    Real      chisq;  /*!< Chi-sq */
  };


  // Prototypes
  //! Linear least squares for real systems
  LLSqResults<EnsemReal> llsq(const LLSqComponent<EnsemReal>& sys, double cov_tol, double solve_tol);

  //! Linear least squares for complex systems
//  LLSqResults llsq(EnsemVectorComplex& x, EnsemReal& chisq, const LLSqComponent<EnsemComplex>& sys);

}  // namespace FF

#endif
