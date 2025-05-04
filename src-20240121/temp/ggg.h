/*! \file
 *  \brief Calculate Weinberg three-gluon operator term
 *
 *  (1/3) f_{ABC} G_{\mu\nu}^A \tilde{G}^{B\nu\lambda} G_{\lambda}^{C\mu}
 *
 */

#ifndef __ggg_h__
#define __ggg_h__

namespace Chroma 
{

  //! Compute Weinberg tree-gluon operator
  /*!
   * \ingroup glue
   *
   * \param u          gauge field (Read)
   * \param ggg_ts     ggg for each timeslice (Write)
   * \param ggg        ggg sum (Write)
   * \param qtop_ts    topological charge for each timeslice (Write)
   * \param qtop       topological charge sum (Write)
   */

  void weinberg_ggg(const multi1d<LatticeColorMatrix>& u, 
     multi1d<DComplex> &ggg_ts, DComplex &ggg,
     multi1d<DComplex> &qtop_ts, DComplex &qtop);

}  // end namespace Chroma

#endif
