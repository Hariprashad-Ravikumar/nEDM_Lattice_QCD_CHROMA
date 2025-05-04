// $Id: formfac_meson_fit.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
// Meson-Meson form-factor code using fitting method

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac.h"
#include "state.h"
#include <list>
#include <iostream>


namespace FF
{
  //! Make a state
  class Meson3PtFunc : public State3PtFunc
  {
  public:
    //! Constructor
    Meson3PtFunc(const std::string& pattern_) : pattern(pattern_) {}

    //! Construct 3pt names
    std::string operator()(int cfg, const ThreePtArg& arg) const
      {
	ArrayInt p_i = arg.pi_pf.p_i;
	ArrayInt p_f = arg.pi_pf.p_f;
	ArrayInt q   = p_f - p_i;
  
	// cout << __func__ << ": pattern=" << pattern.c_str() << std::endl;

	char qxsign = (q[0] < 0 ) ? '-' : '+';
	char qysign = (q[1] < 0 ) ? '-' : '+';
	char qzsign = (q[2] < 0 ) ? '-' : '+';

	char pfxsign = (p_f[0] < 0 ) ? '-' : '+';
	char pfysign = (p_f[1] < 0 ) ? '-' : '+';
	char pfzsign = (p_f[2] < 0 ) ? '-' : '+';

	if (arg.src[0] < 0 || arg.src[0] >= 3)
	{
	  std::cerr << __func__ << ": src out of bounds\n";
	  exit(1);
	}
	char xyz[3] = {'X', 'Y', 'Z'};

	char filen[1024];
	sprintf(filen, pattern.c_str(),
		pfzsign, abs(p_f[2]), pfysign, abs(p_f[1]), pfxsign, abs(p_f[0]),
		cfg,
		xyz[arg.src[0]],
		qzsign, abs(q[2]), qysign, abs(q[1]), qxsign, abs(q[0]));

	std::cout << __func__ << ": filen=" << filen << "\n";

	return std::string(filen);
      }

  private:
    std::string pattern;
  };




  //! Print an array
  std::ostream& operator<<(std::ostream& s, const ArrayInt& d)
  {
    s << d[0];
    for(int i=1; i < d.size(); ++i)
      s << " " << d[i];

    return s;
  }


  //! Print an array
  std::ostream& operator<<(std::ostream& s, const ArrayDouble& d)
  {
    s << d[0];
    for(int i=1; i < d.size(); ++i)
      s << " " << d[i];

    return s;
  }


  //! Params for linear system
  struct MesonFFSysParam
  {
    MesonFFSysParam();
    MesonFFSysParam(XMLReader& xml_in, const std::string& path);

    struct ThreePtParam_t
    {
      string           source_id;
      string           sink_id;
      Array<ArrayInt>  sink_mom;
      GroupXML_t       threept;     /*!< xml holding 3pt group */
    };
    list<ThreePtParam_t>  threept_params;

    /*! 
     * Being bold here - instead of an array of ground+excited states, just 
     * list them all as separate states. This allows only excited state fits with
     * no ground state!
     */
    struct SourceSinkOperator_t
    {
      string       operator_id;   /*!< id in operator factory */
      string       name;          /*!< local name for debugging */
      string       quantum_id;    /*!< quantum number label for overlap-state */
      string       mass_label;    /*!< additional label distinguishing masses */
      GroupXML_t   amp;           /*!< xml holding Z group */
      GroupXML_t   energy;        /*!< xml holding E group */
    };
    list<SourceSinkOperator_t>  source_operators;
    list<SourceSinkOperator_t>  sink_operators;

    struct InsertionOperator_t
    {
      string   type;    /*!< type of operator insertion for matrix element */
      string   renorm;  /*!< ensemble of renormalization */
    } insertion_operator;

    LatticeParam   lattice;

    string    qsq_file_name;

    int       t_source;
    int       t_sink;

    struct MomentumRange
    {
      MomentumRange() {mom_rangeP = false; mom2_max = -1;}

      bool             mom_rangeP;
      int              mom2_max;
      int              pi2_max;
      Array<ArrayInt>  sink_mom;
    } mom_range;

    struct FixedMomentum
    {
      FixedMomentum() {fixed_momP = false;}

      bool             fixed_momP;
      Array<PiPf>      pi_pf;
    } fixed_mom;

    Array<ArrayInt>  fit_range;
  };



  // Three-pt params
  void read(XMLReader& xml_in, const std::string& path, MesonFFSysParam::ThreePtParam_t& param) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "source_id", param.source_id);
    read(paramtop, "sink_id", param.sink_id);
    read(paramtop, "sink_mom", param.sink_mom);
    param.threept = ADATXML::readXMLGroup(paramtop, "ThreePtFuncs", "ThreePtType");
  }


  // Source/sink operators
  void read(XMLReader& xml_in, const std::string& path, MesonFFSysParam::SourceSinkOperator_t& param) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "operator_id", param.operator_id);
    read(paramtop, "name", param.name);
    read(paramtop, "quantum_id", param.quantum_id);
    read(paramtop, "mass_label", param.mass_label);
    read(paramtop, "amp", param.amp);
    read(paramtop, "energy", param.energy);
    param.amp = ADATXML::readXMLGroup(paramtop, "Amplitude", "AmpType");
  }


  // Insertion operator
  void read(XMLReader& xml_in, const std::string& path, MesonFFSysParam::InsertionOperator_t& param) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "type", param.type);
    read(paramtop, "renorm_factor", param.renorm);
  }


  // Meson-meson sys params
  MesonFFSysParam::MesonFFSysParam(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "ThreePts", threept_params);
    read(paramtop, "SourceOperators", source_operators);
    read(paramtop, "SinkOperators", sink_operators);
    read(paramtop, "InsertionOperator", insertion_operator);

    read(paramtop, "LatticeParam", lattice);
    read(paramtop, "qsq_file_name", qsq_file_name);
    read(paramtop, "t_source", t_source);
    read(paramtop, "t_sink", t_sink);
    read(paramtop, "fit_range", fit_range);

    if (paramtop.count("MomentumRange") == 1)
    {
      XMLReader top(paramtop, "MomentumRange");

      mom_range.mom_rangeP = true;
      read(top, "mom2_max", mom_range.mom2_max);
      read(top, "pi2_max", mom_range.pi2_max);    
      read(top, "sink_mom", mom_range.sink_mom);
    }

    if (paramtop.count("FixedMomentum") == 1)
    {
      XMLReader top(paramtop, "FixedMomentum");
      
      fixed_mom.fixed_momP = true;
      read(top, "PiPf", fixed_mom.pi_pf);
    }

    if (mom_range.mom_rangeP && fixed_mom.fixed_momP)
      throw string("both fixed and varying momenta ranges specified - only allow one");
    
    if (! mom_range.mom_rangeP && ! fixed_mom.fixed_momP)
      throw string("neither fixed nor varying momenta ranges specified - must specify one");
  }


  //! System for linear least squares
  /*! System here is rho-pi FF */
  class MesonFFSys: public LLSqComponent<EnsemReal>
  {
  public:
    //! destructor to help with cleanup;
    MesonFFSys(const MesonFFSysParam& param_);

    //! destructor to help with cleanup;
    ~MesonFFSys() {}

    //! Set qsq to use
    void setQsq(const QsqVal& qsq);

    //! Reset time fit interval
    void setTimeFit(int ti, int tf);

    //! Get time fit interval
    ArrayInt getTimeFit() const;

    //! Return the matrix element insertion id
    std::string matrixElementInsertionId() const;

    //! Return the matrix element key for a given source and sink key
    std::string matrixElementId(const std::string& source, const std::string& sink) const;

    //! Should we use this matrix element or ignore it
    bool useMatrixElementP(const std::string& source, const std::string& sink) const;

    //! Get time fit ranges
    const Array<ArrayInt>& getFitRange() const;

    //! Return the row index list
    const std::list<ArrayInt>& getRowList() const {return row_list;}

    //! Return the column index list
    const std::list<ArrayInt>& getColumnList() const {return col_list;}

    //! Given source coordinates, return the corresponding matrix element
    EnsemReal mat(std::list<ArrayInt>::const_iterator row, 
		  std::list<ArrayInt>::const_iterator nFF) const;

    //! Given source coordinates, return the corresponding rhs element
    EnsemReal rhs(std::list<ArrayInt>::const_iterator nFF) const;

    //! Number of bins
    int nbins() const {return nbin;}

    //! Time extent
    int timeLen() const {return param.lattice.latt_size[param.lattice.decay_dir];}

    //! Number of FF
    int nUnknowns() const {return nUnknown;}

    //! Set number of FF
    void setNUnknowns(int n);

    //! Return the dimensions of all the indices
    const ArrayInt& size() const {return sz;}

    //! Energy of initial state at momentum p
    EnsemReal& getMassI(int n) const;

    //! Energy of final state at momentum p
    EnsemReal& getMassF(int m) const;

  protected:
    //! Convenience - make a threept arg
    ThreePtArg threePtArg(const PiPf& pi_pf, int g, int nu) const;

    //! Convenience - make a twopt arg
    TwoPtArg twoPtArg(const ArrayInt& p) const;

    //! Energy of rho at momentum p
    EnsemReal& energy_i(int level, const ArrayInt& p) const;

    //! Energy of pi at momentum p
    EnsemReal& energy_f(int level, const ArrayInt& p) const;

    //! Compute matrix element piece
    EnsemReal matelem(int m, const ArrayInt& p_f, int n, const ArrayInt& p_i, 
		      int mu, int nu) const;

    //! Compute propagation factors
    EnsemReal propfact(int m, const ArrayInt& p_f, int n, const ArrayInt& p_i, int t) const;

  protected:
    //! Hide default onstructor
    MesonFFSys(); 

  private:
    MesonFFSysParam  param;

    //! Hold operator and wavefunctions
    struct OperWaveFunc_t
    {
      Handle< Operator<Complex> >                   op;
      std::list< Handle< WaveFunction<Complex> > >  wvfs;
      mutable Array<ManageEnergy>                   E;
      mutable Array<ManageZAmp>                     Z;
    };

    Array<OperWaveFunc_t>  source_operators;
    Array<OperWaveFunc_t>  sink_operators;
  
    struct ThreePtFunc_t
    {
      string           source_id;
      string           sink_id;
      Array<ArrayInt>  sink_mom;
      mutable Handle<Manage3PtFuncCache>  corr;
    };

    ThreePtFunc_t      threept;

    EnsemReal    op_renorm;

    std::list<ArrayInt>  row_list;
    std::list<ArrayInt>  col_list;

    int              nUnknown;
    int              nbin;
    int              t_i;
    int              t_f;
    Array<PiPf>      pi_pf;
    ArrayInt         sz;
  };



  // Constructor
  MesonFFSys::MesonFFSys(const MesonFFSysParam& param_) : 
    param(param_), 
    threept(new Meson3PtFunc(param.threept.pattern), param.threept.cache_file, param.threept.cfg_file)
  {
    t_i = param.fit_range[0][0];   //! start of fit range
    t_f = param.fit_range[0][1];   //! end of fit range

    // Read insertion-operator renorm factor
    read(param.insertion_operator.renorm);
    nbin = op_renorm.size();
    cout << __func__ << ": found nbins=" << nbin << endl;

    // Read energies
    E_rho.resize(param.E.pattern_rho.size());
    for(int i=0; i < param.E.pattern_rho.size(); ++i)
      E_rho[i].create(new LocalEnergyFunc(param.E.pattern_rho[i], "rho"));

    // Read amplitudes
    Z_rho.resize(param.twopt.pattern_Z_rho.size());
    for(int i=0; i < param.twopt.pattern_Z_rho.size(); ++i)
      Z_rho[i].create(new LocalZFunc(param.twopt.pattern_Z_rho[i], "rho"), param.twopt.avg_2pt_func);

    // Check ground state energies  
    ArrayInt zero(Nd-1);
    zero = 0;
    if (E_rho[0][zero].size() != nbins())
    {
      cerr << __func__ << ": energies inconsistent size with Z_V" << endl;
      exit(1);
    }

    nUnknown = maxNUnknown;   // NEED TO SET THIS
  
    cout << __func__ << ": maximum number of unknowns = " << maxNUnknown << endl;
  }


  //! Convenience - make a threept arg
  ThreePtArg MesonFFSys::threePtArg(const PiPf& pi_pf, int g, int nu) const
  {
    ThreePtArg  arg;
    arg.pi_pf  = pi_pf;
    arg.g      = g;
    arg.src.resize(1);
    arg.src[0] = nu;

//  std::cout << "3pt[]:  g= " << arg.g 
//	    << "   p_f= " << arg.pi_pf.p_f 
//	    << "   p_i= " << arg.pi_pf.p_i
//	    << "   nu= "  << nu
//	    << std::endl;

    return arg;
  }

  //! Convenience - make a twopt arg
  TwoPtArg MesonFFSys::twoPtArg(const ArrayInt& p) const
  {
    TwoPtArg  arg;
    arg.p = p;
    return arg;
  }



  //! Set number of FF
  void MesonFFSys::setNUnknowns(int n) 
  {
    if (n > maxNUnknowns())
    {
      cerr << __func__ << ": too many unknowns\n";
      exit(1);
    }
    nUnknown = n;
  }


  // Reset Qsq to be used
  void MesonFFSys::setQsq(const QsqVal& qsq)
  {
    pi_pf.resize(qsq.pi_pf.size());

    int n = 0;
    for(list<PiPf>::const_iterator nqq=qsq.pi_pf.begin(); 
	nqq != qsq.pi_pf.end(); 
	++nqq,++n)
    {
      pi_pf[n] = *nqq;
    }

    // Index sizes
    sz.resize(4);
    sz[0] = t_f - t_i + 1;   // fitting time range
    sz[1] = Nd;              // directions
    sz[2] = Nd-1;            // number of spatial directions
    sz[3] = pi_pf.size();    // number of momenta pairs
  }


  //! Reset time fit interval
  void MesonFFSys::setTimeFit(int ti, int tf) 
  {
    t_i = ti; t_f = tf;

    sz[0] = t_f - t_i + 1;   // reset time range
  }


  //! Get time fit interval
  ArrayInt MesonFFSys::getTimeFit() const
  {
    ArrayInt tt(2);
    tt[0] = t_i;
    tt[1] = t_f;

    return tt;
  }


  //! Return the matrix element insertion id
  std::string MesonFFSys::matrixElementInsertionId() const
  {
    return param.insertion_operator;
  }


  //! Return the matrix element key for a given source and sink key
  std::string MesonFFSys::matrixElementId(const std::string& source, const std::string& sink) const
  {
    return source + "|" + matrixElementInsertionId() + "|" + sink;
  }


  //! Should we use this matrix element or ignore it
  bool MesonFFSys::useMatrixElementP(const std::string& source, const std::string& sink) const
  {
    return (TheMatrixElementFactory::Instance().exist(matrixElementId(source, sink))) ? true : false;
  }



  //! Get time fit interval
  const Array<ArrayInt>& MesonFFSys::getFitRange() const
  {
    return param.fit_range;
  }


  //! Compute or return energy
  EnsemReal& MesonFFSys::energy_i(int level, const ArrayInt& p0) const
  {
//  cout << __func__ << endl;

    ArrayInt p = canonicalOrder(p0);
    ArrayInt zero(Nd-1);
    zero = 0;

#if 0
    if (p != zero && ! E_rho[level].exist(p))
      E_rho[level].insert(p, energyDisp(E_rho[level][zero], p, param.lattice));
#endif

    return E_rho[level][p];
  }


  //! Compute or return energy
  EnsemReal& MesonFFSys::energy_f(int level, const ArrayInt& p0) const
  {
//  cout << __func__ << endl;

    ArrayInt p = canonicalOrder(p0);
    ArrayInt zero(Nd-1);
    zero = 0;

#if 0
    if (p != zero && ! E_pion[level].exist(p))
      E_pion[level].insert(p, energyDisp(E_pion[level][zero], p, param.lattice));
#endif

    return E_pion[level][p];
  }


  //! Energy of initial state at momentum p
  EnsemReal& MesonFFSys::getMassI(int n) const
  {
    ArrayInt zero(Nd-1);
    zero = 0;
    return E_rho[n][zero];
  }


  //! Energy of final state at momentum p
  EnsemReal& MesonFFSys::getMassF(int m) const
  {
    ArrayInt zero(Nd-1);
    zero = 0;
    return E_pion[m][zero];
  }

 
  //! Compute matrix element piece
  EnsemReal MesonFFSys::matelem(int m, const ArrayInt& p_f, int n, const ArrayInt& p_i, 
				int mu, int nu) const
  {
    if (nu >= Nd-1)
    {
      cerr << __func__ << ": currently only support mu=0..3, nu<3\n";
      exit(1);
    }
  
    const double xi  = param.lattice.aniso.xi;
    ArrayDouble pp_i = contMom(p_i,param.lattice.latt_size[0]) / xi;
    ArrayDouble pp_f = contMom(p_f,param.lattice.latt_size[0]) / xi;

    ArrayInt zero(Nd-1);
    zero = 0;

    const EnsemReal& m_f = getMassF(m);
    const EnsemReal& m_i = getMassI(n);
    const EnsemReal& E_f = energy_f(m,p_f);
    const EnsemReal& E_i = energy_i(n,p_i);
    EnsemReal        rat = Real(2) / (m_i + m_f);

    // Transition matrix element
    EnsemReal f;
    f.resize(nbins());
    f = 0;

    // Really dumb method
    TensorReal anti = antiSymTensor<4>();
    Array<int> n1(4); 

    for(int alpha=0; alpha < Nd; ++alpha)
      for(int delta=0; delta < Nd; ++delta)
      {
	n1[0] = mu; n1[1] = alpha; n1[2] = delta; n1[3] = nu;

	Real a = peekObs(anti,n1);

	// printf("%s: mu=%d alpha=%d delta=%d nu=%d:     a= %d\n",
	//   __func__,mu,alpha,delta,nu,int(toFloat(a)));

	if (toFloat(a) != 0)
	{
	  if (alpha == Nd-1)
	  {
	    f += a * rat*E_f*Real(pp_i[delta]);
	  }
	  else if (delta == Nd-1)
	  {
	    f += a * rat*Real(pp_f[alpha])*E_i;
	  }
	  else 
	  {
	    f += a * rat*Real(pp_f[alpha])*Real(pp_i[delta]);
	  }

	  //printf("%s: mu=%d alpha=%d delta=%d nu=%d: pf=%d,%d,%d, pi=%d,%d,%d: a=%d f=%f\n",
	  // __func__,mu,alpha,delta,nu,
	  //  p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2],
	  // int(toFloat(a)),toFloat(mean(f)));

	}
      }

    return f;
  }



  //! Compute propagation factors
  EnsemReal MesonFFSys::propfact(int m, const ArrayInt& p_f, int n, const ArrayInt& p_i, int t) const
  {
    int T_i = param.t_source;
    int T_f = param.t_sink;
    int T   = param.lattice.latt_size[param.lattice.decay_dir];

    EnsemReal g, g_i, g_f;
    // now we've folded over the second half of the lattice we simplify the images (in fact they should be negligable)
    g_i = exp(-Real(t - T_i)*energy_i(n,p_i));       
    g_f = exp(-Real(abs(T_f - t))*energy_f(m,p_f));  
    g  = Z_pion[m][twoPtArg(p_f)] * Z_rho[n][twoPtArg(p_i)] * g_i * g_f;
    g /= Real(4) * energy_f(m,p_f) * energy_i(n,p_i);

    return g;
  }


  // Given source coordinates, return the corresponding matrix element
  EnsemReal MesonFFSys::mat(std::list<ArrayInt>::const_iterator ind, 
			    std::list<ArrayInt>::const_iterator nFF) const
  {
//  cout << __func__ << endl;

    // NOTE: encoding here:   nFF has indices like  int[numTemporalFF][2]
    // So, if m and n encode the energy level, then one indexes like
    //
    // (nFF & 1) == 0 spatial V
    // (nFF & 3) == 1 temporal V
    // nFF = <spatial/temporal V> + 2*( n*n + m)
    //
    int t  = ind[0] + t_i;
    int mu = ind[1];
    int nu = ind[2];
    int p  = ind[3];

    ArrayInt p_i = pi_pf[p].p_i;
    ArrayInt p_f = pi_pf[p].p_f;

    int m, n;
    decompose(m, n, nFF); // Determine sink,source energy levels

    EnsemReal f;
    f.resize(nbins());
    f = 0;

    // spatial V
    if ((nFF == 0)  && (mu < Nd-1))
    {
      f = matelem(m, p_f, n, p_i, mu, nu) * propfact(m, p_f, n, p_i, t);
    }
    if ((nFF == 2)  && (mu < Nd-1))
    {
      f = matelem(m, p_f, n, p_i, mu, nu) * propfact(m, p_f, n, p_i, t);
    }
    if ((nFF == 4)  && (mu < Nd-1))
    {
      f = matelem(m, p_f, n, p_i, mu, nu) * propfact(m, p_f, n, p_i, t);
    }
    if ((nFF == 6) && (mu < Nd-1))
    {
      f = matelem(m, p_f, n, p_i, mu, nu) * propfact(m, p_f, n, p_i, t);
    }

    // temporal V
    if ((nFF == 1) && (mu == Nd-1))
    {
      f = matelem(m, p_f, n, p_i, mu, nu) * propfact(m, p_f, n, p_i, t);
    }
    if ((nFF == 3)  && (mu == Nd-1))
    {
      f = matelem(m, p_f, n, p_i, mu, nu) * propfact(m, p_f, n, p_i, t);
    }
    if ((nFF == 5)  && (mu == Nd-1))
    {
      f = matelem(m, p_f, n, p_i, mu, nu) * propfact(m, p_f, n, p_i, t);
    }
    if ((nFF == 7)  && (mu == Nd-1))
    {
      f = matelem(m, p_f, n, p_i, mu, nu) * propfact(m, p_f, n, p_i, t);
    }


    // Final kinetic factor
    f /= Z_V;

//  cout << "exiting mat" << endl;

    return f;
  }


  // Given source coordinates, return the corresponding rhs element
  EnsemReal MesonFFSys::rhs(std::list<ArrayInt>::const_iterator ind) const
  {
//  cout << __func__ << endl;

    int t  = ind[0] + t_i;
    int mu = ind[1];
    int nu = ind[2];
    int p  = ind[3];

    if (nu < 0 || nu >= Nd-1)
    {
      cerr << __func__ << ": nu out of bounds\n";
      exit(1);
    }

    EnsemVectorComplexF& three = threept[threePtArg(pi_pf[p], 1 << mu, nu)];
    EnsemReal thr;

    // Compensate for minus sign in rho_y from sequential source construction
    // This sign is independent of mu or momenta
    int sg = (nu == 1)? -1 : 1;

    // Extract appropriate part of 3pt. Fold in rho_y sign above.
    if (mu == Nd-1)
      thr = - Real(0.5)*Real(sg)*real(   peekObs(three,t) + peekObs(three, 48-t)    );
    else
      thr = - Real(0.5)*Real(sg)*imag(   peekObs(three,t) - peekObs(three, 48-t)    );
    //the minus sign accounts for the use of euclidean antisym tensor rather than mink.

//  cout << "exiting rhs" << endl;

    return thr;
  }


  // Setup all the operators and wavefunctions
  MesonFFSys::setup()
  {
    //
    // Create the source operator and wave function overlaps
    //
    source_operators.resize(param.source_operators.size());
    for(int src_op=0; src_op < param.source_operators.size(); ++src_op)
    {
      string source_id = param.source_operators[src_op].operator_id;
      cout << "Source operator = " << source_id << endl;

//    const string operator_path = "/FileTest";
      source_operators[src_op].op = TheMesonOperatorFactory::Instance().createObject(source_id,
										     xml_in,
										     operator_path);

      // Construct the wave functions
      source_operators[src_op].wvfs = source_operators[src_op].op->overlaps();
    }

    //
    // Create the sink operator and wave function overlaps
    //
    sink_operators.resize(param.sink_operators.size());
    for(int snk_op=0; snk_op < param.sink_operators.size(); ++snk_op)
    {
      string sink_id = param.sink_operators[snk_op].operator_id;
      cout << "Sink operator = " << sink_id << endl;

//    const string operator_path = "/FileTest";
      sink_operators[snk_op].op = TheMesonOperatorFactory::Instance().createObject(sink_id,
										   xml_in,
										   operator_path);

      // Construct the wave functions
      sink_operators[src_op].wvfs = sink_operators[src_op].op->overlaps();
    }


    //
    // Loop over source and sink wavefunctions, find the name (quantum number)
    // of the state and find allowed transition matrix elements
    //
    for(int src_op=0; src_op < param.source_operators.size(); ++src_op)
    {
      string source_op_name = source_operators[src_op].op->operatorName();
      cout << "source_op_name = " << source_op_name << endl;

      for(int snk_op=0; snk_op < param.sink_operators.size(); ++snk_op)
      {
	string sink_op_name = sink_operators[src_op].op->operatorName();
	cout << "sink_op_name = " << sink_op_name << endl;

	typedef list< Handle< WaveFunction<Complex> > > ListType_t;

	for(ListType_t::const_iterator src_wvf=param.source_operators[src_op].wvfs.begin();
	    src_wvf != param.source_operators[src_op].wvfs.end(); 
	    ++src_wvf)
	{
	  const WaveFunction<Complex>& source_wvf = *(*src_wvf);
	  string source_particle = source_wvf.particleName();
	  cout << "source_particle = " << source_particle << endl;

	  for(ListType_t::const_iterator src_wvf=param.sink_operators[snk_op].wvfs.begin();
	      snk_wvf != param.sink_operators[snk_op].wvfs.end(); 
	      ++snk_wvf)
	  {
	    const WaveFunction<Complex>& sink_wvf = *(*snk_wvf);
	    string sink_particle = sink_wvf.particleName();
	    cout << "sink_particle = " << sink_particle << endl;
    
	    // Use if the quantum numbers match matrix element
	    string mat_elem_id = matrixElementId(source_particle, sink_particle);
	    if (! useMatrixElementP(source_particle, sink_particle))
	    {
	      cout << __func__ << ": no matrix element for key = " << mat_elem_id << endl;
	      continue;
	    }

	    // We know the matrix element exist and intend to use it
	    Handle<MesonMatrixElement> mes_mat_elem(
	      TheMesonOperatorFactory::Instance().createObject(mat_elem_id,
							       xml_in,
							       operator_path));

	  }
	}
      }
    }
  }

} // namespace FF


//using namespace ADATXML;
using namespace ENSEM;
using namespace FF;
using namespace std;


int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    cerr << "Usage:"
	 << argv[0] << ":  <input xml file>   <output xml file>" << endl;
    exit(1);
  }

  XMLReader xml_in(argv[1]);
  MesonFFSysParam linParam(xml_in, "/MesonFFFit");

  XMLFileWriter xml_out(argv[2]);
  push(xml_out, "MesonFFFit");

  try
  {
    cout << "construct linear system" << endl;
    MesonFFSys sys(linParam);

    cout << "sys.nbins=" << sys.nbins() << endl;

    // The momenta list
    std::list<QsqVal> qsq_val;

    // Use some method to construct Qsq list
    if (linParam.mom_range.mom_rangeP)
    {
      cout << "Construct qsq list over momenta range" << endl;
      QsqParam qsqParam;
      qsqParam.nq_norm_sq_max = linParam.mom_range.mom2_max;
      qsqParam.lattice        = linParam.lattice;
      qsqParam.mass_i         = toDouble(mean(sys.getMassI(0)));
      qsqParam.mass_f         = toDouble(mean(sys.getMassF(0)));
  
      for(int i=0; i < linParam.mom_range.sink_mom.size(); ++i)
	qsqParam.nf_list.push_back(linParam.mom_range.sink_mom[i]);

      // Construct list
      qsq_val = constructQsq(qsqParam);

#if 1
      // Cannot use 0 momentum
      //   {
      //     ArrayInt zero(Nd-1);
      //    zero = 0;
      //    if (qsq_val.front().pi_pf.front().p_f == zero &&
      //	  qsq_val.front().pi_pf.front().p_i == zero)
      //   {
      //	cout << "First momenta has p_f=p_i=zero - drop since it will cause problems\n" << endl;
      //	qsq_val.pop_front();
      //    }
      //  }
#endif

      // Print out list
      {
	std::ofstream f(linParam.qsq_file_name.c_str());
	f << qsq_val;
	f.close();
      }
    }
    else if (linParam.fixed_mom.fixed_momP)
    {
      cout << "Use fixed momentum list" << endl;

      ArrayInt p_f = linParam.fixed_mom.pi_pf[0].p_f;
      ArrayInt p_i = linParam.fixed_mom.pi_pf[0].p_i;
      const double mass_f = toDouble(mean(sys.getMassF(0)));
      const double mass_i = toDouble(mean(sys.getMassI(0)));

      QsqVal  qsq;
      qsq.Qsq_c = -qsqContDispPhys(mass_f, p_f, mass_i, p_i, linParam.lattice);
      qsq.Qsq_l = -qsqLattDispPhys(mass_f, p_f, mass_i, p_i, linParam.lattice);
      qsq.qsq   = norm2(p_f - p_i);

      for(int i=0; i < linParam.fixed_mom.pi_pf.size(); ++i)
	qsq.pi_pf.push_back(linParam.fixed_mom.pi_pf[i]);

      qsq_val.push_back(qsq);
    }
    else
    {
      cerr << argv[0] << ": invalid momentum specification\n";
      exit(1);
    }

    push(xml_out, "VaryUnknowns");
    int i = sys.nUnknowns();
    {
      push(xml_out, "elem");
      sys.setNUnknowns(i);
      cout << "\nNumber of unknowns set to  " << sys.nUnknowns() << endl;
      cerr << "\nNumber of unknowns set to  " << sys.nUnknowns() << endl;
      ff_driver(xml_out, sys, qsq_val);
      pop(xml_out);
    }
    pop(xml_out);

  }
  catch(std::bad_cast) 
  {
    cerr << argv[0] << ": caught cast error" << endl;
    exit(1);
  }
  catch(const std::string& e) 
  {
    cerr << argv[0] << ": Caught Exception: " << e << endl;
    exit(1);
  }
  catch(std::exception& e) 
  {
    cerr << argv[0] << ": Caught standard library exception: " << e.what() << endl;
    exit(1);
  }
  catch(...)
  {
    cerr << argv[0] << ": caught generic exception during measurement" << endl;
    exit(1);
  }

  pop(xml_out);  // MesonFF
  xml_out.close();
}
