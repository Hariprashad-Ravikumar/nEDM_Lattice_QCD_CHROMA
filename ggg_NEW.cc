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
  
  
  multi1d<LatticeColorMatrix> Wline_down(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
      multi1d<LatticeColorMatrix> Wline_d(m);
      Wline_d[0] = u[mu]; 
      for (int i = 1; i < m; i++){
        Wline_d[i] = shift(Wline_d[i - 1], FORWARD, mu);  
      }
      return Wline_d;
  }
  
  multi1d<LatticeColorMatrix> Wline_left(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
      multi1d<LatticeColorMatrix> Wline_l(n);
      Wline_l[0] = u[nu]; 
      for (int i = 1; i < n; i++){
        Wline_l[i] = shift(Wline_l[i - 1], FORWARD, nu);  
      }
      return Wline_l;
  }
  
  multi1d<LatticeColorMatrix> Wline_right(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
     multi1d<LatticeColorMatrix> Wline_r(n);
     LatticeColorMatrix Wline_d_p1 = u[nu];  

     for (int sh = 0; sh < m; sh++) {
         Wline_d_p1 = shift(Wline_d_p1, FORWARD, mu);
     }
     Wline_r[0] = Wline_d_p1;
     for (int i = 1; i < n; i++) {
         Wline_r[i] = shift(Wline_r[i - 1], FORWARD, nu);
     }
     return Wline_r;
  }
  
  multi1d<LatticeColorMatrix> Wline_top(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
     multi1d<LatticeColorMatrix> Wline_t(m);
     LatticeColorMatrix Wline_l_p1 = u[mu];  

     for (int sh = 0; sh < n; sh++) {
         Wline_l_p1 = shift(Wline_l_p1, FORWARD, nu);
     }
     Wline_t[0] = Wline_l_p1;
     for (int i = 1; i < m; i++) {
         Wline_t[i] = shift(Wline_t[i - 1], FORWARD, mu);
     }
     return Wline_t;
  }
  
  LatticeColorMatrix Wilson_loop1(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
    multi1d<LatticeColorMatrix> Wline_d = Wline_down(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_r = Wline_right(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_l = Wline_left(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_t = Wline_top(u, m, n, mu, nu);

    // attach adj() to all elements of Wline_l
    for (int i = 0; i < Wline_l.size(); i++){
        Wline_l[i] = adj(Wline_l[i]);
    }
    //attach adj() to all elements of Wline_t
    for (int j = 0; j < Wline_t.size(); j++){
        Wline_t[j] = adj(Wline_t[j]);
    }
    LatticeColorMatrix Wl1_d = Wline_d[0];
    for (int a = 1; a < Wline_d.size(); a++){
        Wl1_d = Wl1_d * Wline_d[a];
    }
    LatticeColorMatrix Wl1_r = Wline_r[0];
    for (int b = 1; b < Wline_r.size(); b++){
        Wl1_r = Wl1_r * Wline_r[b];
    }
    LatticeColorMatrix Wl1_t = Wline_t[0];
    for (int c = 1; c < Wline_t.size(); c++){
        Wl1_t = Wline_t[c] * Wl1_t;
    }
    LatticeColorMatrix Wl1_l = Wline_l[0];
    for (int d = 1; d < Wline_l.size(); d++){
        Wl1_l = Wline_l[d] * Wl1_l;
    }
    LatticeColorMatrix Wl1 = Wl1_d * Wl1_r * Wl1_t * Wl1_l;
    return Wl1;
  }
  
  LatticeColorMatrix Wilson_loop2(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
    multi1d<LatticeColorMatrix> Wline_d = Wline_down(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_r = Wline_right(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_l = Wline_left(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_t = Wline_top(u, m, n, mu, nu);

    // attach adj() to all elements of Wline_l
    for (int i = 0; i < Wline_l.size(); i++){
        Wline_l[i] = adj(Wline_l[i]);
    }
    //attach adj() to all elements of Wline_t
    for (int j = 0; j < Wline_t.size(); j++){
        Wline_t[j] = adj(Wline_t[j]);
    }
    LatticeColorMatrix Wl1_d = Wline_d[0];
    for (int a = 1; a < Wline_d.size(); a++){
        Wl1_d = Wl1_d * Wline_d[a];
    }
    LatticeColorMatrix Wl1_r = Wline_r[0];
    for (int b = 1; b < Wline_r.size(); b++){
        Wl1_r = Wl1_r * Wline_r[b];
    }
    LatticeColorMatrix Wl1_t = Wline_t[0];
    for (int c = 1; c < Wline_t.size(); c++){
        Wl1_t = Wline_t[c] * Wl1_t;
    }
    LatticeColorMatrix Wl1_l = Wline_l[0];
    for (int d = 1; d < Wline_l.size(); d++){
        Wl1_l = Wline_l[d] * Wl1_l;
    }
    LatticeColorMatrix Wl2 =Wl1_l * Wl1_d * Wl1_r * Wl1_t;
    for (int sh = 0; sh < n; sh++) {
         Wl2 = shift(Wl2, BACKWARD, nu);
     }
    return Wl2;
  }
  
  LatticeColorMatrix Wilson_loop3(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
    multi1d<LatticeColorMatrix> Wline_d = Wline_down(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_r = Wline_right(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_l = Wline_left(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_t = Wline_top(u, m, n, mu, nu);

    // attach adj() to all elements of Wline_l
    for (int i = 0; i < Wline_l.size(); i++){
        Wline_l[i] = adj(Wline_l[i]);
    }
    //attach adj() to all elements of Wline_t
    for (int j = 0; j < Wline_t.size(); j++){
        Wline_t[j] = adj(Wline_t[j]);
    }
    LatticeColorMatrix Wl1_d = Wline_d[0];
    for (int a = 1; a < Wline_d.size(); a++){
        Wl1_d = Wl1_d * Wline_d[a];
    }
    LatticeColorMatrix Wl1_r = Wline_r[0];
    for (int b = 1; b < Wline_r.size(); b++){
        Wl1_r = Wl1_r * Wline_r[b];
    }
    LatticeColorMatrix Wl1_t = Wline_t[0];
    for (int c = 1; c < Wline_t.size(); c++){
        Wl1_t = Wline_t[c] * Wl1_t;
    }
    LatticeColorMatrix Wl1_l = Wline_l[0];
    for (int d = 1; d < Wline_l.size(); d++){
        Wl1_l = Wline_l[d] * Wl1_l;
    }
    LatticeColorMatrix Wl3 = Wl1_r * Wl1_t * Wl1_l * Wl1_d;
    for (int sh = 0; sh < m; sh++) {
         Wl3 = shift(Wl3, BACKWARD, mu);
     }
    return Wl3;
  }
  
  LatticeColorMatrix Wilson_loop4(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
    multi1d<LatticeColorMatrix> Wline_d = Wline_down(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_r = Wline_right(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_l = Wline_left(u, m, n, mu, nu);
    multi1d<LatticeColorMatrix> Wline_t = Wline_top(u, m, n, mu, nu);

    // attach adj() to all elements of Wline_l
    for (int i = 0; i < Wline_l.size(); i++){
        Wline_l[i] = adj(Wline_l[i]);
    }
    //attach adj() to all elements of Wline_t
    for (int j = 0; j < Wline_t.size(); j++){
        Wline_t[j] = adj(Wline_t[j]);
    }
    LatticeColorMatrix Wl1_d = Wline_d[0];
    for (int a = 1; a < Wline_d.size(); a++){
        Wl1_d = Wl1_d * Wline_d[a];
    }
    LatticeColorMatrix Wl1_r = Wline_r[0];
    for (int b = 1; b < Wline_r.size(); b++){
        Wl1_r = Wl1_r * Wline_r[b];
    }
    LatticeColorMatrix Wl1_t = Wline_t[0];
    for (int c = 1; c < Wline_t.size(); c++){
        Wl1_t = Wline_t[c] * Wl1_t;
    }
    LatticeColorMatrix Wl1_l = Wline_l[0];
    for (int d = 1; d < Wline_l.size(); d++){
        Wl1_l = Wline_l[d] * Wl1_l;
    }
    LatticeColorMatrix Wl4 = Wl1_t * Wl1_l * Wl1_d * Wl1_r;
    for (int sh = 0; sh < m; sh++) {
         Wl4 = shift(Wl4, BACKWARD, mu);
     }
    for (int sh = 0; sh < n; sh++) {
         Wl4 = shift(Wl4, BACKWARD, nu);
     }
    return Wl4;
  }
  
  LatticeColorMatrix clover_loop_mxn(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
    LatticeColorMatrix Wl1 = Wilson_loop1(u, m, n, mu, nu);  
    LatticeColorMatrix Wl2 = Wilson_loop2(u, m, n, mu, nu);  
    LatticeColorMatrix Wl3 = Wilson_loop3(u, m, n, mu, nu);  
    LatticeColorMatrix Wl4 = Wilson_loop4(u, m, n, mu, nu);  
    return (Wl1 + Wl2 + Wl3 + Wl4);
  }
  
  LatticeColorMatrix C_mxn(const multi1d<LatticeColorMatrix>& u, int m, int n, int mu, int nu)
  {
    LatticeColorMatrix cloop_mxn = clover_loop_mxn(u, m, n, mu, nu);  
    LatticeColorMatrix cloop_nxm = clover_loop_mxn(u, n, m, mu, nu);  
    return (1.0/8.0)*(cloop_mxn + cloop_nxm);
  }

  void calculate_FieldStrength(const multi1d<LatticeColorMatrix>& u, int mu1, int nu1, LatticeColorMatrix &Fmunu_lat)
  {
    Real k5 = 0.0;

    LatticeColorMatrix c1x1 = C_mxn(u, 1, 1, mu1, nu1);
    LatticeColorMatrix c2x2 = C_mxn(u, 2, 2, mu1, nu1);
    LatticeColorMatrix c1x2 = C_mxn(u, 1, 2, mu1, nu1);
    LatticeColorMatrix c1x3 = C_mxn(u, 1, 3, mu1, nu1);
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
     
    Fmunu_lat = (19.0/9.0-55.0*k5)*c1x1 + (1.0/36.0-16.0*k5)*c2x2 + (64.0*k5-32.0/45.0)*c1x2 + (1.0/15.0-6.0*k5)*c1x3;
    // Take traceless antihermitian projection
   // W' = (1/2) [W - W^dag - (1/Nc) Tr(W - W^dag)]
    tmp_1=1;
    tmp_2 = adj(Fmunu_lat);
    Fmunu_lat -= tmp_2; //W - W^dag
    tmp_2 = tmp_1 * trace(Fmunu_lat)/Nc; // (1/3) Tr(W - W^dag)
    Fmunu_lat -= tmp_2;
    Fmunu_lat *= 0.5;

   // Multiply -i to make F
    Fmunu_lat *= cmplx(Real(0), -Real(1));
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
    
    QDPIO::cout << "New ggg sum = " << ggg << std::endl;


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
    QDPIO::cout << "New qtop sum = " << qtop << std::endl;


    END_CODE();
  }

}  // end namespace Chroma
