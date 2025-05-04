/*! \file
 *  \brief Calculate Weinberg three-gluon operator term
 *
 *  (1/3) f_{ABC} G_{\mu\nu}^A \tilde{G}^{B\nu\lambda} G_{\lambda}^{C\mu}
 *
 */

#include "chromabase.h"
#include "meas/glue/ggg.h"

namespace Chroma 
{
  class TimeSliceFunc : public SetFunc
  {
  public:
    TimeSliceFunc(int dir): dir_decay(dir) {}
  
    int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
    int numSubsets() const {return Layout::lattSize()[dir_decay];}
  
    int dir_decay;
  
  private:
    TimeSliceFunc() {}  // hide default constructor
  };

  // Calculate lattice version of the continuum field-strength tensor following
  // Highly-improved lattice field-strength tensor 
  // (Bilson-Thompson, et al., https://arxiv.org/pdf/hep-lat/0203008.pdf)
  // Code is taken from qnaive.cc
  void calculate_FieldStrength(const multi1d<LatticeColorMatrix>& u, int mu1, int nu1, LatticeColorMatrix &u_clov_1)
  {
    // k5 is a tunable free paramter.
    // A 4-loop improved field-strength can be obtained by setting k5 among
    // the following three choices: k5 = 0, k5 = 19/495, or k5 = 1/576
    // Here we set k5 = 0.0 but the code (taken from qnaive.cc) can handle any k5.
    Real k5 = 0.0;

    /* Local Variables */
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;

    Real k1,k2,k3,k4,kk5;

    if( Nd != 4 )
      QDP_error_exit("Nd for the ggg has to be 4 but: ", Nd);

    k1 = 19.0/9.0 - 55.0 * k5;
    k1 *= 2.0;
    k2 = 1.0/36.0 - 16.0 * k5;
    k2 *= 2.0;
    k3 = 64.0 * k5 - 32.0/45.0;
    k4 = 1.0/15.0 - 6.0 * k5;
    kk5 = k5;
    kk5 *= 2.0;

    if( toBool(k1 != 0) ) {
      /* First "plus-plus" 1x1 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 = k1 * tmp_1;

      /* First "plus-minus" 1x1 */
     // tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
     // tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
     // tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
     // u_clov_1 -= k1 * tmp_1;
      /* First "plus-minus" 1x1 */
     tmp_3 = shift(adj(u[nu1]), FORWARD, mu1);
     tmp_1 = u[mu1] * shift(tmp_3, BACKWARD, nu1);

     tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
     tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
     u_clov_1 -= k1 * tmp_1;


      /* First "minus-minus" 1x1 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 += k1 * tmp_1;
      /* First "minus-minus" 1x1 */
      tmp_1 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_2 = shift(tmp_1, BACKWARD, nu1);
      tmp_3 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_1 = tmp_3 * tmp_2;

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = tmp_1 * tmp_2;

      tmp_1 = tmp_3 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k1 * tmp_1;


      ///* First "minus-plus" 1x1 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 -= k1 * tmp_1;
      /* First "minus-plus" 1x1 */
      tmp_2 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_3 = shift(u[nu1], BACKWARD, mu1);
      tmp_1 = tmp_2 * tmp_3;

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k1 * tmp_1;


    }

    if( toBool(k2!=0) ) {
      ///* First "plus-plus" 2x2 */
      //tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);

      //tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 += k2 * tmp_1;
      /* First "plus-plus" 2x2 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_3 = shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, mu1);

      tmp_3 = shift(u[nu1], FORWARD, mu1);
      tmp_1 = shift(tmp_3, FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * tmp_3;

      tmp_2 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k2 * tmp_1;


      ///* First "plus-minus" 2x2 */
      //tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 -= k2 * tmp_1;
      /* First "plus-minus" 2x2 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);

      tmp_2 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_1 = shift(tmp_3, FORWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[nu1], BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k2 * tmp_1;


      ///* First "minus-minus" 2x2 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 += k2 * tmp_1;
      /* First "minus-minus" 2x2 */
      tmp_2 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_3;

      tmp_2 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = tmp_2 * shift(tmp_1, BACKWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_3 * shift(tmp_1, BACKWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[nu1], BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k2 * tmp_1;


      ///* First "minus-plus" 2x2 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 -= k2 * tmp_1;
      /* First "minus-plus" 2x2 */
      tmp_3 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_2;

      tmp_3 = shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, mu1);

      tmp_1 = shift(u[nu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k2 * tmp_1;
    }

    if( toBool(k3!=0) ) {
      ///* First "plus-plus" 2x1 */
      //tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 += k3 * tmp_1;
      /* First "plus-plus" 2x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_3 = shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, mu1);

      tmp_3 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k3 * tmp_1;


      ///* First "plus-minus" 2x1 */
      //tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 -= k3 * tmp_1;
      /* First "plus-minus" 2x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k3 * tmp_1;


      ///* First "minus-minus" 2x1 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]),BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[mu1],BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 += k3 * tmp_1;
      /* First "minus-minus" 2x1 */
      tmp_2 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_3;

      tmp_2 = shift(adj(u[nu1]),BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(u[mu1],BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k3 * tmp_1;


      ///* First "minus-plus" 2x1 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 -= k3 * tmp_1;
      /* First "minus-plus" 2x1 */
      tmp_2 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_3;

      tmp_3 =shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, mu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k3 * tmp_1;



      ///* First "plus-plus" 1x2 */
      //tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 += k3 * tmp_1;
      /* First "plus-plus" 1x2 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_3 = shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k3 * tmp_1;


      ///* First "plus-minus" 1x2 */
      //tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 -= k3 * tmp_1;
      /* First "plus-minus" 1x2 */
      tmp_2 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_1 = u[mu1] * shift(tmp_2, BACKWARD, nu1);

      tmp_2 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[nu1], BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k3 * tmp_1;


      ///* First "minus-minus" 1x2 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 += k3 * tmp_1;
      /* First "minus-minus" 1x2 */
      tmp_2 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_3;

      tmp_2 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[nu1], BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k3 * tmp_1;


      ///* First "minus-plus" 1x2 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 -= k3 * tmp_1;
      /* First "minus-plus" 1x2 */
      tmp_2 = shift(u[nu1], BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_2;

      tmp_3 = shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k3 * tmp_1;


    }

    if( toBool(k4!=0) ) {
      ///* First "plus-plus" 3x1 */
      //tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 += k4 * tmp_1;
      /* First "plus-plus" 3x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_3 = shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, mu1);

      tmp_1 = shift(u[nu1], FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, mu1);

      tmp_2 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k4 * tmp_1;


      ///* First "plus-minus" 3x1 */
      //tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 -= k4 * tmp_1;
      /* First "plus-minus" 3x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_3 = shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, mu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_1 = shift(tmp_3, FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k4 * tmp_1;


      ///* First "minus-minus" 3x1 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 += k4 * tmp_1;
      /* First "minus-minus" 3x1 */
      tmp_3 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_2;

      tmp_2 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, mu1);

      tmp_3 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k4 * tmp_1;


      ///* First "minus-plus" 3x1 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 -= k4 * tmp_1;
      /* First "minus-plus" 3x1 */
      tmp_3 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_2;

      tmp_2 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, mu1);

      tmp_1 = shift(u[nu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, mu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k4 * tmp_1;



      ///* First "plus-plus" 1x3 */
      //tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 += k4 * tmp_1;
      /* First "plus-plus" 1x3 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_3 = shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = shift(u[nu1], FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = shift(adj(u[mu1]), FORWARD, nu1);
      tmp_3 = shift(tmp_2, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k4 * tmp_1;


      ///* First "plus-minus" 1x3 */
      //tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 -= k4 * tmp_1;
      /* First "plus-minus" 1x3 */
      tmp_3 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_1 = u[mu1] * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(u[nu1], BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[nu1], BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k4 * tmp_1;


      ///* First "minus-minus" 1x3 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]),BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 += k4 * tmp_1;
      /* First "minus-minus" 1x3 */
      tmp_3 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, nu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_2;

      tmp_2 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[nu1]),BACKWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(u[nu1], BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[nu1], BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k4 * tmp_1;


      ///* First "minus-plus" 1x3 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 -= k4 * tmp_1;
      /* First "minus-plus" 1x3 */
      tmp_2 = shift(u[nu1], BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_2;

      tmp_3 = shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = shift(u[nu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, FORWARD, nu1);
      tmp_3 = shift(tmp_2, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k4 * tmp_1;


    }

    if( toBool(kk5!=0) ) {
      ///* First "plus-plus" 3x3 */
      //tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      //tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      //tmp_1 = tmp_2 * adj(u[nu1]);
      //u_clov_1 += kk5 * tmp_1;
      /* First "plus-plus" 3x3 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_3 = shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, mu1);

      tmp_1 = shift(u[nu1], FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, mu1);

      tmp_3 = shift(u[nu1], FORWARD, mu1);
      tmp_2 = shift(tmp_3, FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = shift(u[nu1], FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, mu1);
      tmp_1 = shift(tmp_3, FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, mu1);
      tmp_2 = shift(tmp_3, FORWARD, nu1);
      tmp_3 = shift(tmp_2, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_1 = shift(tmp_3, FORWARD, nu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = shift(adj(u[mu1]), FORWARD, nu1);
      tmp_3 = shift(tmp_2, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += kk5 * tmp_1;


      ///* First "plus-minus" 3x3 */
      //tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 -= kk5 * tmp_1;
      /* First "plus-minus" 3x3 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_3 = shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, mu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_1 = shift(tmp_3, FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, mu1);
      tmp_2 = shift(tmp_3, FORWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, mu1);
      tmp_1 = shift(tmp_3, FORWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_3 = shift(tmp_2, FORWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[mu1]), FORWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(u[nu1], BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[nu1], BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= kk5 * tmp_1;

      ///* First "minus-minus" 3x3 */
      //tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      //tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      //tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      //u_clov_1 += kk5 * tmp_1;
      /* First "minus-minus" 3x3 */
      tmp_3 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_2;

      tmp_2 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, mu1);

      tmp_3 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_2 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(adj(u[nu1]), BACKWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, nu1);
      tmp_3 = shift(tmp_2, BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);

      tmp_1 = shift(u[nu1], BACKWARD, nu1);
      tmp_3 = shift(tmp_1, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, nu1);

      tmp_3 = shift(u[nu1], BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += kk5 * tmp_1;


      /* First "minus-plus" 3x3 */
      tmp_3 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * tmp_2;

      tmp_2 = shift(adj(u[mu1]), BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, BACKWARD, mu1);

      tmp_1 = shift(u[nu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(tmp_3, BACKWARD, mu1);

      tmp_3 = shift(u[nu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = shift(u[nu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, BACKWARD, mu1);
      tmp_3 = shift(tmp_2, BACKWARD, mu1);
      tmp_2 = shift(tmp_3, FORWARD, nu1);
      tmp_3 = shift(tmp_2, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_1 = shift(u[mu1], BACKWARD, mu1);
      tmp_3 = shift(tmp_1, BACKWARD, mu1);
      tmp_1 = shift(tmp_3, FORWARD, nu1);
      tmp_3 = shift(tmp_1, FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(u[mu1], BACKWARD, mu1);
      tmp_2 = shift(tmp_3, FORWARD, nu1);
      tmp_3 = shift(tmp_2, FORWARD, nu1);
      tmp_2 = tmp_1 * shift(tmp_3, FORWARD, nu1);

      tmp_3 = shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(tmp_3, FORWARD, nu1);

      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= kk5 * tmp_1;

    }

   // Take traceless antihermitian projection
   // W' = (1/2) [W - W^dag - (1/Nc) Tr(W - W^dag)]
   tmp_1=1;
   tmp_2 = adj(u_clov_1);
   u_clov_1 -= tmp_2; //W - W^dag
   tmp_2 = tmp_1 * trace(u_clov_1)/Nc; // (1/3) Tr(W - W^dag)
   u_clov_1 -= tmp_2;
   u_clov_1 *= 0.5;

   // Multiply 1/8 to make i*F_\mu\nu from the plaquette sum
   u_clov_1 *= 0.125;

   // Multiply -i to make F
   u_clov_1 *= cmplx(Real(0), -Real(1));
  }



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
     multi1d<DComplex> &qtop_ts, DComplex &qtop)
  {
    START_CODE();

    //---------------------------
    // Calculate F_{\mu\nu}
    //---------------------------
    multi2d<LatticeColorMatrix> F(Nd,Nd);
    for(int mu=0; mu<Nd; ++mu)
    for(int nu=0; nu<Nd; ++nu)
      calculate_FieldStrength(u, mu, nu, F[mu][nu]);

    //---------------------------
    // Calculate F^A = Tr(F T^A)
    //---------------------------
    if( Nc != 3 )
      QDP_error_exit("Nc for the ggg has to be 3 but: ", Nc);

    // Build SU(3) generators
    // 1/2 are included in the definitions of belows
    Complex  one = cmplx( Real(0.5), Real(0.0));
    Complex mone = cmplx(-Real(0.5), Real(0.0));
    Complex  img = cmplx( Real(0.0), Real(0.5));
    Complex mimg = cmplx( Real(0.0),-Real(0.5));

    Complex invsqrt3     = cmplx( Real(0.5)/sqrt(Real(3.0)), Real(0.0));
    Complex mtwoinvsqrt3 = cmplx(-Real(1.0)/sqrt(Real(3.0)), Real(0.0));

    // SU(3) generator = (1/2) * Gell-Mann matrix
    multi1d<ColorMatrix> T(8);
    T = zero;

    pokeColor(T[0], one,  0, 1);
    pokeColor(T[0], one,  1, 0);

    pokeColor(T[1], mimg, 0, 1);
    pokeColor(T[1], img,  1, 0);

    pokeColor(T[2], one,  0, 0);
    pokeColor(T[2], mone, 1, 1);

    pokeColor(T[3], one,  0, 2);
    pokeColor(T[3], one,  2, 0);

    pokeColor(T[4], mimg, 0, 2);
    pokeColor(T[4], img,  2, 0);

    pokeColor(T[5], one,  1, 2);
    pokeColor(T[5], one,  2, 1);

    pokeColor(T[6], mimg, 1, 2);
    pokeColor(T[6],  img, 2, 1);

    pokeColor(T[7], invsqrt3,      0, 0);
    pokeColor(T[7], invsqrt3,      1, 1);
    pokeColor(T[7], mtwoinvsqrt3,  2, 2);

    // Calculate F^A = Tr(F T^A)
    multi3d<LatticeComplex> FA(Nd,Nd,8);
    for(int mu=0; mu<Nd; ++mu)
    for(int nu=0; nu<Nd; ++nu)
    for(int A=0;  A<8;   ++A)
      FA[mu][nu][A] = traceColor(F[mu][nu]*T[A])*Real(2);

    // Structure constant
    multi3d<Real> fABC(8,8,8);
    for(int A=0; A<8; ++A)
    for(int B=0; B<8; ++B)
    for(int C=0; C<8; ++C)
      fABC[A][B][C] = Real(0);

    fABC[0][1][2] =  Real(1);
    fABC[0][2][1] = -Real(1);
    fABC[0][3][6] =  Real(0.5);
    fABC[0][4][5] = -Real(0.5);
    fABC[0][5][4] =  Real(0.5);
    fABC[0][6][3] = -Real(0.5);
    fABC[1][0][2] = -Real(1);
    fABC[1][2][0] =  Real(1);
    fABC[1][3][5] =  Real(0.5);
    fABC[1][4][6] =  Real(0.5);
    fABC[1][5][3] = -Real(0.5);
    fABC[1][6][4] = -Real(0.5);
    fABC[2][0][1] =  Real(1);
    fABC[2][1][0] = -Real(1);
    fABC[2][3][4] =  Real(0.5);
    fABC[2][4][3] = -Real(0.5);
    fABC[2][5][6] = -Real(0.5);
    fABC[2][6][5] =  Real(0.5);
    fABC[3][0][6] = -Real(0.5);
    fABC[3][1][5] = -Real(0.5);
    fABC[3][2][4] = -Real(0.5);
    fABC[3][4][2] =  Real(0.5);
    fABC[3][4][7] =  sqrt(Real(3))/Real(2);
    fABC[3][5][1] =  Real(0.5);
    fABC[3][6][0] =  Real(0.5);
    fABC[3][7][4] = -sqrt(Real(3))/Real(2);
    fABC[4][0][5] =  Real(0.5);
    fABC[4][1][6] = -Real(0.5);
    fABC[4][2][3] =  Real(0.5);
    fABC[4][3][2] = -Real(0.5);
    fABC[4][3][7] = -sqrt(Real(3))/Real(2);
    fABC[4][5][0] = -Real(0.5);
    fABC[4][6][1] =  Real(0.5);
    fABC[4][7][3] =  sqrt(Real(3))/Real(2);
    fABC[5][0][4] = -Real(0.5);
    fABC[5][1][3] =  Real(0.5);
    fABC[5][2][6] =  Real(0.5);
    fABC[5][3][1] = -Real(0.5);
    fABC[5][4][0] =  Real(0.5);
    fABC[5][6][2] = -Real(0.5);
    fABC[5][6][7] =  sqrt(Real(3))/Real(2);
    fABC[5][7][6] = -sqrt(Real(3))/Real(2);
    fABC[6][0][3] =  Real(0.5);
    fABC[6][1][4] =  Real(0.5);
    fABC[6][2][5] = -Real(0.5);
    fABC[6][3][0] = -Real(0.5);
    fABC[6][4][1] = -Real(0.5);
    fABC[6][5][2] =  Real(0.5);
    fABC[6][5][7] = -sqrt(Real(3))/Real(2);
    fABC[6][7][5] =  sqrt(Real(3))/Real(2);
    fABC[7][3][4] =  sqrt(Real(3))/Real(2);
    fABC[7][4][3] = -sqrt(Real(3))/Real(2);
    fABC[7][5][6] =  sqrt(Real(3))/Real(2);
    fABC[7][6][5] = -sqrt(Real(3))/Real(2);
 
    // Epsilon tensor
    // Structure constant
    multi4d<Real> epsilon(4,4,4,4);
    for(int i=0; i<4; ++i)
    for(int j=0; j<4; ++j)
    for(int k=0; k<4; ++k)
    for(int l=0; l<4; ++l)
      epsilon[i][j][k][l] = Real(0);

    epsilon[0][1][2][3] =  Real(1);
    epsilon[0][1][3][2] = -Real(1);
    epsilon[0][2][1][3] = -Real(1);
    epsilon[0][2][3][1] =  Real(1);
    epsilon[0][3][1][2] =  Real(1);
    epsilon[0][3][2][1] = -Real(1);
    epsilon[1][0][2][3] = -Real(1);
    epsilon[1][0][3][2] =  Real(1);
    epsilon[1][2][0][3] =  Real(1);
    epsilon[1][2][3][0] = -Real(1);
    epsilon[1][3][0][2] = -Real(1);
    epsilon[1][3][2][0] =  Real(1);
    epsilon[2][0][1][3] =  Real(1);
    epsilon[2][0][3][1] = -Real(1);
    epsilon[2][1][0][3] = -Real(1);
    epsilon[2][1][3][0] =  Real(1);
    epsilon[2][3][0][1] =  Real(1);
    epsilon[2][3][1][0] = -Real(1);
    epsilon[3][0][1][2] = -Real(1);
    epsilon[3][0][2][1] =  Real(1);
    epsilon[3][1][0][2] =  Real(1);
    epsilon[3][1][2][0] = -Real(1);
    epsilon[3][2][0][1] = -Real(1);
    epsilon[3][2][1][0] =  Real(1);

    LatticeComplex ggg_tmp = zero;

    // ggg = (1/3) f_{ABC} G_{ij}^A \tilde{G}_{jk}^{B} G_{ki}^{C}
    //     = (1/3) f_{ABC} G_{ij}^A (1/2) epsilon_{jklm} G_{lm}^{B} G_{ki}^{C}
    // Note that (1/2) in front of epsilon is taken out by calculating only m > l
    for(int A=0; A<8; ++A)
    for(int B=0; B<8; ++B)
    for(int C=0; C<8; ++C)
    for(int i=0; i<4; ++i)
    for(int j=0; j<4; ++j)
    for(int k=0; k<4; ++k)
    for(int l=0; l<4; ++l)
    for(int m=l+1; m<4; ++m)
    {
      if(toBool(epsilon[j][k][l][m] == 0)) continue;
      if(toBool(fABC[A][B][C]) == 0)       continue;

      ggg_tmp += fABC[A][B][C] * epsilon[j][k][l][m] * ((FA[i][j][A] * FA[l][m][B]) * FA[k][i][C]);
    }
    ggg_tmp *= (Real(1)/Real(3));

    Set timeslice;
    timeslice.make(TimeSliceFunc(Nd-1));

    const multi1d<int>& latt_size = Layout::lattSize();
    ggg_ts.resize(latt_size[Nd-1]);
    ggg_ts = sumMulti(ggg_tmp, timeslice);
    ggg = sum(ggg_tmp);
    
    QDPIO::cout << "ggg = " << ggg << std::endl;


    //----------------------------------------------------------------------------------------------
    // Calculate topological charge 
    //----------------------------------------------------------------------------------------------
    // q = (1/32/pi^2) epsilon_ijkl Tr(G_ij Gkl)
    // Additional factor Real(4) is multiplied for calculating only j>i and l>k components
    // Note that this is not the optimal calculation; one can reduce number of operators further 
    // as in qnaive.cc
    LatticeComplex qtop_tmp = zero;
    for(int i=0; i<4; ++i)
    for(int j=i+1; j<4; ++j)
    for(int k=0; k<4; ++k)
    for(int l=k+1; l<4; ++l)
    {
      if(toBool(epsilon[i][j][k][l] == 0)) continue;
      qtop_tmp += epsilon[i][j][k][l] * trace(F[i][j] * F[k][l]);
    }

    qtop_tmp *= Real(4)/(Real(8)*twopi*twopi);

    qtop_ts.resize(latt_size[Nd-1]);
    qtop_ts = sumMulti(qtop_tmp, timeslice);
    qtop = sum(qtop_tmp);
    QDPIO::cout << "qtop = " << qtop << std::endl;


    END_CODE();
  }

}  // end namespace Chroma
