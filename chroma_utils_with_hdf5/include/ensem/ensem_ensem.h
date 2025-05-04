// -*- C++ -*-
// $Id: ensem_ensem.h,v 2.0 2008/12/05 04:43:33 edwards Exp $

/*! \file
 * \brief Ensemble layer
 */

#define ENSEM_USE_JACKKNIFE
#undef ENSEM_USE_BOOTSTRAP

namespace ENSEM {

//-------------------------------------------------------------------------------------
//! The types of ensembles
enum EnsemType_t
{
  ENSEM_JACKKNIFE = 10,
  ENSEM_BOOTSTRAP = 11
};



//-------------------------------------------------------------------------------------
/*! \addtogroup escalar Ensemble scalar
 * \ingroup fiberbundle
 *
 * Ensemble scalars
 *
 * @{
 */

//! Scalar ensem
/*! All objects are of EScalar or Ensem type */
template<class T> class EScalar
{
public:
  EScalar() {}
  ~EScalar() {}

  //---------------------------------------------------------
  //! construct dest = const
  EScalar(const typename WordType<T>::Type_t& rhs)
    {
      typedef typename InternalScalar<T>::Type_t  Scalar_t;
      elem() = Scalar_t(rhs);
    }

  //! construct dest = 0
  EScalar(const Zero& rhs)
    {
      assign(rhs);
    }

  //! conversion by constructor  EScalar<T> = EScalar<T1>
  template<class T1>
  EScalar(const EScalar<T1>& rhs)
    {
      assign(rhs);
    }

  //! construct dest = T1
  template<class T1>
  EScalar(const T1& rhs)
    {
      elem() = rhs;
    }


protected:
  //---------------------------------------------------------
  //! construct dest = const
  EScalar& assign(const typename WordType<T>::Type_t& rhs)
    {
      typedef typename InternalScalar<T>::Type_t  Scalar_t;
      elem() = Scalar_t(rhs);
      return *this;
    }


  //! construct dest = 0
  EScalar& assign(const Zero& rhs)
    {
      zero_rep(*this);
      return *this;
    }


  //! conversion by assignment  EScalar<T> = EScalar<T1>
  template<class T1>
  EScalar& assign(const EScalar<T1>& rhs)
    {
      elem() = rhs.elem();
      return *this;
    }


public:
  //---------------------------------------------------------
  //! dest = const
  /*! Fill with an integer constant. Will be promoted to underlying word type */
  inline
  EScalar& operator=(const typename WordType<T>::Type_t& rhs)
    {
      return assign(rhs);
    }

  //! dest = zero
  inline
  EScalar& operator=(const Zero& rhs)
    {
      return assign(rhs);
    }

  //! EScalar = EScalar
  /*! Set equal to another EScalar */
  inline
  EScalar& operator=(const EScalar& rhs) 
    {
      return assign(rhs);
    }

  //! EScalar = EScalar
  /*! Set equal to another EScalar */
  template<class T1>
  inline
  EScalar& operator=(const EScalar<T1>& rhs) 
    {
      return assign(rhs);
    }

  //! EScalar += EScalar
  template<class T1>
  inline
  EScalar& operator+=(const EScalar<T1>& rhs) 
    {
      elem() += rhs.elem();
      return *this;
    }

  //! EScalar -= EScalar
  template<class T1>
  inline
  EScalar& operator-=(const EScalar<T1>& rhs) 
    {
      elem() -= rhs.elem();
      return *this;
    }

  //! EScalar *= EScalar
  template<class T1>
  inline
  EScalar& operator*=(const EScalar<T1>& rhs) 
    {
      elem() *= rhs.elem();
      return *this;
    }

  //! EScalar /= EScalar
  template<class T1>
  inline
  EScalar& operator/=(const EScalar<T1>& rhs) 
    {
      elem() /= rhs.elem();
      return *this;
    }

  //! EScalar %= EScalar
  template<class T1>
  inline
  EScalar& operator%=(const EScalar<T1>& rhs) 
    {
      elem() %= rhs.elem();
      return *this;
    }

  //! EScalar |= EScalar
  template<class T1>
  inline
  EScalar& operator|=(const EScalar<T1>& rhs) 
    {
      elem() |= rhs.elem();
      return *this;
    }

  //! EScalar &= EScalar
  template<class T1>
  inline
  EScalar& operator&=(const EScalar<T1>& rhs) 
    {
      elem() &= rhs.elem();
      return *this;
    }

  //! EScalar ^= EScalar
  template<class T1>
  inline
  EScalar& operator^=(const EScalar<T1>& rhs) 
    {
      elem() ^= rhs.elem();
      return *this;
    }

  //! EScalar <<= EScalar
  template<class T1>
  inline
  EScalar& operator<<=(const EScalar<T1>& rhs) 
    {
      elem() <<= rhs.elem();
      return *this;
    }

  //! EScalar >>= EScalar
  template<class T1>
  inline
  EScalar& operator>>=(const EScalar<T1>& rhs) 
    {
      elem() >>= rhs.elem();
      return *this;
    }

  //! Copy constructor
  EScalar(const EScalar& rhs) : F(rhs.F) {}


public:
  //! Return ref to an element via int
  const EScalar operator[](int off) const {return F[off];}

public:
  //! Has this object been initialized (resized and such)
  inline bool initP() const {return elem().initP();}

  //! Number of bins
  inline int size() const {return 1;}

  //! Number of elements of observable
  inline int numElem() const {return elem().numElem();}

  inline void resize(const EScalar& a) 
    {
      elem().resize(a.elem());
    }
  inline void resizeObs(int n1) 
    {
      elem().resize(n1);
    }
  inline void resizeObs(int n1, int n2)
    {
      elem().resize(n1,n2);
    }
  inline void resizeObs(int n1, int n2, int n3)
    {
      elem().resize(n1,n2,n3);
    }
  inline void resizeObs(const Array<int>& nz)
    {
      elem().resize(nz);
    }

public:
  T& elem() {return F;}
  const T& elem() const {return F;}

private:
  T F;
};


 
// Input
//! Ascii input
template<class T>
inline
std::istream& operator>>(std::istream& s, EScalar<T>& d)
{
  s >> d.elem();
  return s;
}

//! Ascii output
template<class T> 
inline  
std::ostream& operator<<(std::ostream& s, const EScalar<T>& d)
{
  return s << d.elem();
}


//! Text input
template<class T>
inline
ENSEM::TextReader& operator>>(ENSEM::TextReader& s, EScalar<T>& d)
{
  return s >> d.elem();
}

//! Text output
template<class T> 
inline  
ENSEM::TextWriter& operator<<(ENSEM::TextWriter& s, const EScalar<T>& d)
{
  return s << d.elem();
}


//! XML output
template<class T>
inline
ADATXML::XMLWriter& operator<<(ADATXML::XMLWriter& xml, const EScalar<T>& d)
{
  return xml << d.elem();
}

//! XML input
template<class T>
inline
void read(ADATXML::XMLReader& xml, const std::string& path, EScalar<T>& d)
{
  read(xml, path, d.elem());
}

//! XML output
template<class T>
inline
void write(ADATXML::XMLWriter& xml, const std::string& s, const EScalar<T>& d)
{
  xml.openTag(s);
  xml << d;
  xml.closeTag();
}

/*! @} */  // end of group escalar



//-------------------------------------------------------------------------------------
// Forward decl
double rescaleFactor(int N, EnsemType_t ensem_type);


//! Ensemble class
/*!
 * \addtogroup eensem Ensemble class
 * \ingroup fiberbundle 
 * @{
 */
template<class T> class Ensem
{
public:
  Ensem() {ensem_type=ENSEM_JACKKNIFE; N=0; F=0;}
  ~Ensem() {delete[] F;}

  //---------------------------------------------------------
  //! Force the ensemble to have certain statistics
  //! Return the proper rescaling factor
  double rescaleFactor() const {return ENSEM::rescaleFactor(N, ensem_type);}

  //---------------------------------------------------------
  //! Force the ensemble to have certain statistics
  void setEnsemType(EnsemType_t p) {ensem_type = p;}

  //! Return the ensemble type
  EnsemType_t getEnsemType() const {return ensem_type;}

  //---------------------------------------------------------
  //! conversion by constructor  Ensem<T> = Ensem<T1>
  template<class T1>
  Ensem(const Ensem<T1>& rhs) : N(0), F(0)
    {
      checkResize("convert Ensem", rhs.size(), rhs.getEnsemType());
      assign(rhs);
    }


protected:
  //---------------------------------------------------------
  //! conversion by constructor  Ensem<T> = EScalar<T1>
  template<class T1>
  Ensem& assign(const EScalar<T1>& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem();
      return *this;
    }

  //! conversion by constructor  Ensem<T> = Ensem<T1>
  template<class T1>
  Ensem& assign(const Ensem<T1>& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem(i);
      return *this;
    }

  //! construct Ensem = const
  Ensem& assign(const typename WordType<T>::Type_t& rhs)
    {
      typedef typename InternalScalar<T>::Type_t  Scalar_t;
      for(int i=0; i < N; ++i)
	elem(i) = Scalar_t(rhs);
      return *this;
    }

  //! construct Ensem = 0
  Ensem& assign(const Zero& rhs)
    {
      zero_rep(*this);
      return *this;
    }

public:
  //---------------------------------------------------------
  //! dest = const
  /*! Fill with an integer constant. Will be promoted to underlying word type */
  inline
  Ensem& operator=(const typename WordType<T>::Type_t& rhs)
    {
      checkSize("Ensem = const");
      return assign(rhs);
    }

  //! dest = zero
  inline
  Ensem& operator=(const Zero& rhs)
    {
      checkSize("Ensem = zero");
      return assign(rhs);
    }

  //---------------------------------------------------------
  //! Ensem = EScalar
  /*! Set equal to an EScalar */
  template<class T1>
  inline
  Ensem& operator=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem = EScalar");
      return assign(rhs);
    }

  //! Ensem += EScalar
  template<class T1>
  inline
  Ensem& operator+=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem += EScalar");
      for(int i=0; i < N; ++i)
	elem(i) += rhs.elem();

      return *this;
    }

  //! Ensem -= EScalar
  template<class T1>
  inline
  Ensem& operator-=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem -= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) -= rhs.elem();

      return *this;
    }

  //! Ensem *= EScalar
  template<class T1>
  inline
  Ensem& operator*=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem *= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) *= rhs.elem();

      return *this;
    }

  //! Ensem /= EScalar
  template<class T1>
  inline
  Ensem& operator/=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem /= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) /= rhs.elem();

      return *this;
    }

  //! Ensem %= EScalar
  template<class T1>
  inline
  Ensem& operator%=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem *= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) %= rhs.elem();

      return *this;
    }

  //! Ensem |= EScalar
  template<class T1>
  inline
  Ensem& operator|=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem |= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) |= rhs.elem();

      return *this;
    }

  //! Ensem &= EScalar
  template<class T1>
  inline
  Ensem& operator&=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem &= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) &= rhs.elem();

      return *this;
    }

  //! Ensem ^= EScalar
  template<class T1>
  inline
  Ensem& operator^=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem ^= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) ^= rhs.elem();

      return *this;
    }

  //! Ensem <<= EScalar
  template<class T1>
  inline
  Ensem& operator<<=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem <<= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) <<= rhs.elem();

      return *this;
    }

  //! Ensem >>= EScalar
  template<class T1>
  inline
  Ensem& operator>>=(const EScalar<T1>& rhs) 
    {
      checkSize("Ensem >>= EScalar");
      for(int i=0; i < N; ++i)
	elem(i) >>= rhs.elem();

      return *this;
    }

  //---------------------------------------------------------
  //! Ensem = Ensem
  /*! Set equal to another Ensem */
  inline
  Ensem& operator=(const Ensem& rhs) 
    {
      checkResize("Ensem = Ensem", rhs.size(), rhs.getEnsemType());
      ensem_type = rhs.getEnsemType();
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem(i);

      return *this;
    }

  //! Ensem = Ensem
  /*! Set equal to another Ensem */
  template<class T1>
  inline
  Ensem& operator=(const Ensem<T1>& rhs) 
    {
      checkResize("Ensem = Ensem<T1>", rhs.size(), rhs.getEnsemType());
      ensem_type = rhs.getEnsemType();
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem(i);

      return *this;
    }

  //! Ensem += Ensem
  template<class T1>
  inline
  Ensem& operator+=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem += Ensem", rhs.size(), rhs.getEnsemType());
      for(int i=0; i < N; ++i)
	elem(i) += rhs.elem(i);

      return *this;
    }

  //! Ensem -= Ensem
  template<class T1>
  inline
  Ensem& operator-=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem -= Ensem", rhs.size(), rhs.getEnsemType());
      for(int i=0; i < N; ++i)
	elem(i) -= rhs.elem(i);

      return *this;
    }

  //! Ensem *= Ensem
  template<class T1>
  inline
  Ensem& operator*=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem *= Ensem", rhs.size(), rhs.getEnsemType());

      Ensem<T>  d;
      d.checkResize(__func__, rhs.size(), rhs.getEnsemType());

      Ensem<T>  ll(rescaleEnsemDown(*this));
      Ensem<T1> rr(rescaleEnsemDown(rhs));

      for(int i=0; i < d.size(); ++i)
	d.elem(i) = ll.elem(i) * rr.elem(i);

      return (*this = rescaleEnsemUp(d));
    }

  //! Ensem /= Ensem
  template<class T1>
  inline
  Ensem& operator/=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem /= Ensem", rhs.size(), rhs.getEnsemType());

      Ensem<T>  d;
      d.checkResize(__func__, rhs.size(), rhs.getEnsemType());

      Ensem<T>  ll(rescaleEnsemDown(*this));
      Ensem<T1> rr(rescaleEnsemDown(rhs));

      for(int i=0; i < d.size(); ++i)
	d.elem(i) = ll.elem(i) / rr.elem(i);

      return (*this = rescaleEnsemUp(d));
    }

#if 0
  // THESE FUNCTIONS DO NOT MAKE SENSE IN AN ENSEMBLE
  // COMMENT OUT FOR NOW UNTIL THEY MIGHT BE NEEDED
  // POSSIBLY STILL NEED RESCALING

  //! Ensem %= Ensem
  template<class T1>
  inline
  Ensem& operator%=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem %= Ensem", rhs.size(), rhs.getEnsemType());
      for(int i=0; i < N; ++i)
	elem(i) %= rhs.elem(i);

      return *this;
    }

  //! Ensem |= Ensem
  template<class T1>
  inline
  Ensem& operator|=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem |= Ensem", rhs.size(), rhs.getEnsemType());
      for(int i=0; i < N; ++i)
	elem(i) |= rhs.elem(i);

      return *this;
    }

  //! Ensem &= Ensem
  template<class T1>
  inline
  Ensem& operator&=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem &= Ensem", rhs.size(), rhs.getEnsemType());
      for(int i=0; i < N; ++i)
	elem(i) &= rhs.elem(i);

      return *this;
    }

  //! Ensem ^= Ensem
  template<class T1>
  inline
  Ensem& operator^=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem ^= Ensem", rhs.size(), rhs.getEnsemType());
      for(int i=0; i < N; ++i)
	elem(i) ^= rhs.elem(i);

      return *this;
    }

  //! Ensem <<= Ensem
  template<class T1>
  inline
  Ensem& operator<<=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem <<= Ensem", rhs.size(), rhs.getEnsemType());
      for(int i=0; i < N; ++i)
	elem(i) <<= rhs.elem(i);

      return *this;
    }

  //! Ensem >>= Ensem
  template<class T1>
  inline
  Ensem& operator>>=(const Ensem<T1>& rhs) 
    {
      checkSize("Ensem >>= Ensem", rhs.size(), rhs.getEnsemType());
      for(int i=0; i < N; ++i)
	elem(i) >>= rhs.elem(i);

      return *this;
    }
#endif

  //---------------------------------------------------------
  //! Copy constructor
  Ensem(const Ensem& rhs) : N(0), F(0)
    {
      checkResize("Ensem convert ensem", rhs.size(), rhs.getEnsemType());
      assign(rhs);
    }


  //---------------------------------------------------------
public:
public:
  inline void checkSize(const char *s) const
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkSize: %s\n", s);
#endif

      if (size() == 0)
      {
	std::cerr << s << ": Invalid Ensem size" << std::endl;
	exit(1);
      }
    }

  inline void checkSize(const char *s, int n1, EnsemType_t e1) const
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkSize: %s, Ensem[%d]\n", s, n1);
#endif

      if (size() == 0 && size() == n1)
      {
	std::cerr << s << ": Invalid Ensem dest and/or source size" << std::endl;
	exit(1);
      }
    }

  inline void checkResize(const char *s, int n1, EnsemType_t e1)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, Ensem[%d]\n", s, n1);
#endif

      if (n1 == 0)
      {
	std::cerr << "checkResize: " << s << ": invalid Ensem source sizes: n1=" 
		  << n1 << std::endl;
	exit(1);
      }
      resize(n1);
      setEnsemType(e1);
    }

  inline void checkResize(const char *s, int n1, EnsemType_t e1, int n2, EnsemType_t e2)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, Ensem[%d,%d]\n", s, n1, n2);
#endif

      if (n1 == 0 || n1 != n2)
      {
	std::cerr << "checkResize: " << s << ": invalid Ensem source sizes: n1,n2=" 
		  << n1 << " " << n2 << std::endl;
	exit(1);
      }
      resize(n1);

      if (e1 != e2)
      {
	std::cerr << "checkResize: " << s << ": invalid Ensemble types: e1,e2=" 
		  << e1 << " " << e2 << std::endl;
	exit(1);
      }
      setEnsemType(e1);
    }

  inline void checkResize(const char *s, 
			  int n1, EnsemType_t e1, 
			  int n2, EnsemType_t e2, 
			  int n3, EnsemType_t e3)
    {
#if ENSEM_DEBUG >= 3
      fprintf(stdout,"checkResize: %s, Ensem[%d,%d,%d]\n", s, n1, n2, n3);
#endif

      if (n1 == 0 || n1 != n2 || n1 != n3)
      {
	std::cerr << "checkResize: " << s << ": invalid Ensem source sizes: n1,n2,n3=" 
		  << n1 << " " << n2 << " " << n3 << std::endl;
	exit(1);
      }
      resize(n1);

      if (e1 != e2 || e1 != e3)
      {
	std::cerr << "checkResize: " << s << ": invalid Ensemble types: e1,e2,n3=" 
		  << e1 << " " << e2 << " " << e3 << std::endl;
	exit(1);
      }
      setEnsemType(e1);
    }

  //! Has this object been initialized (resized and such)
  bool initP() const 
    {
      return ((N != 0) ? true : false) & elem().initP();
    }

  //! Number of configurations
  inline int size() const {return N;}

  //! Number of elements of observable
  inline int numElem() const {return elem(0).numElem();}

  //! Resize the number of configs
  inline void resize(int n)
    {
      if (N > 0)
      {
	delete[] F;
	N = 0;
	F = 0;
      }
      N = n;      
      F = new(std::nothrow) T[n];
      if ( F == 0x0 ) 
      { 
	std::cerr << "Unable to allocate memory in Ensem::resize(): n = " << n << std::endl;
	exit(1);
      }
#if ENSEM_DEBUG >= 1
      fprintf(stdout,"Ensem[%d]=0x%x, this=0x%x, bytes/elem=%d\n",
	      N,F,this,sizeof(T));
#endif
    }

public:
  inline void resize(const Ensem& a) 
    {
      resize(a.size());
      for(int i=0; i < N; ++i)
	elem(i).resize(a.elem(i));
    }
  inline void resizeObs(int n1) 
    {
      checkSize(__func__);
      for(int i=0; i < N; ++i)
	elem(i).resize(n1);
    }
  inline void resizeObs(int n1, int n2)
    {
      checkSize(__func__);
      for(int i=0; i < N; ++i)
	elem(i).resize(n1,n2);
    }
  inline void resizeObs(int n1, int n2, int n3)
    {
      checkSize(__func__);
      for(int i=0; i < N; ++i)
	elem(i).resize(n1,n2,n3);
    }
  inline void resizeObs(const Array<int>& nz)
    {
      checkSize(__func__);
      for(int i=0; i < N; ++i)
	elem(i).resize(nz);
    }

public:
  inline T& elem(int i) {return F[i];}
  inline const T& elem(int i) const {return F[i];}

public:
  //! Debugging info
  void print_info(char *name)
    {
      fprintf(stdout,"Info: %s = Ensem[%d]=0x%x, this=0x%x\n",
	      name,N,(void *)F,this);
    }


  //---------------------------------------------------------
public:
  //! The backdoor
  /*! 
   * Used by optimization routines (e.g., SSE) that need the memory address of data.
   * BTW: to make this a friend would be a real pain since functions are templatized.
   */
  inline const T* getF() const {return F;}
  inline T* getF() {return F;}

private:
  int N;
  T *F;
  EnsemType_t ensem_type;
};



#if 0
//! Stream input
template<class T>
inline
std::istream& operator>>(std::istream& s, Ensem<T>& d)
{
  for(int i=0; i < d.size(); ++i)
    s >> d.elem(i);

  return s;
}

//! Stream output
template<class T>
inline
std::ostream& operator<<(std::ostream& s, const Ensem<T>& d)
{
  for(int i=0; i < d.size(); ++i)
    s << d.elem(i) << "\n";
  return s;
}
#endif


//! Text input
template<class T>
inline
ENSEM::TextReader& operator>>(ENSEM::TextReader& s, Ensem<T>& d)
{
  for(int i=0; i < d.size(); ++i)
    s >> d.elem(i);
  return s;
}

//! Text output
template<class T> 
inline  
ENSEM::TextWriter& operator<<(ENSEM::TextWriter& s, const Ensem<T>& d)
{
  for(int i=0; i < d.size(); ++i)
    s << d.elem(i) << "\n";
  return s;
}


/*! @} */   // end of group eensem


//-----------------------------------------------------------------------------
// Traits classes to support operations of simple scalars (floating constants, 
// etc.) on EnsemTypes
//-----------------------------------------------------------------------------

template<class T>
struct WordType<EScalar<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};

template<class T>
struct WordType<Ensem<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};


// Internally used scalars
template<class T>
struct InternalScalar<EScalar<T> > {
  typedef EScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

template<class T>
struct InternalScalar<Ensem<T> > {
  typedef EScalar<typename InternalScalar<T>::Type_t>  Type_t;
};


// Trait to make a primitive scalar leaving grid along
template<class T>
struct PrimitiveScalar<EScalar<T> > {
  typedef EScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

template<class T>
struct PrimitiveScalar<Ensem<T> > {
  typedef Ensem<typename PrimitiveScalar<T>::Type_t>  Type_t;
};


// Trait to make a ensemble scalar leaving other indices alone
template<class T>
struct EnsemScalar<EScalar<T> > {
  typedef EScalar<typename EnsemScalar<T>::Type_t>  Type_t;
};

template<class T>
struct EnsemScalar<Ensem<T> > {
  typedef EScalar<typename EnsemScalar<T>::Type_t>  Type_t;
};


// Internally used real scalars
template<class T>
struct RealScalar<EScalar<T> > {
  typedef EScalar<typename RealScalar<T>::Type_t>  Type_t;
};

template<class T>
struct RealScalar<Ensem<T> > {
  typedef Ensem<typename RealScalar<T>::Type_t>  Type_t;
};


// Traits class to label IO types
template<class T> 
struct EnsbcIO<EScalar<T> > {
  enum {type = EnsbcIO<T>::type};
};

template<class T> 
struct EnsbcIO<Ensem<T> > {
  enum {type = EnsbcIO<T>::type};
};



//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(EScalar) -> EScalar
template<class T1, class Op>
struct UnaryReturn<EScalar<T1>, Op> {
  typedef EScalar<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};

// Default unary(Ensem) -> Ensem
template<class T1, class Op>
struct UnaryReturn<Ensem<T1>, Op> {
  typedef Ensem<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};

// Default binary(EScalar,EScalar) -> EScalar
template<class T1, class T2, class Op>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, Op> {
  typedef EScalar<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Currently, the only trinary operator is ``where'', so return 
// based on T2 and T3
// Default trinary(EScalar,EScalar,EScalar) -> EScalar
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<EScalar<T1>, EScalar<T2>, EScalar<T3>, Op> {
  typedef EScalar<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};

// Default binary(Ensem,Ensem) -> Ensem
template<class T1, class T2, class Op>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, Op> {
  typedef Ensem<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(EScalar,Ensem) -> Ensem
template<class T1, class T2, class Op>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, Op> {
  typedef Ensem<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(Ensem,EScalar) -> Ensem
template<class T1, class T2, class Op>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, Op> {
  typedef Ensem<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};


// Currently, the only trinary operator is ``where'', so return 
// based on T2 and T3

// Default trinary(Ensem,Ensem,Ensem) -> Ensem
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<Ensem<T1>, Ensem<T2>, Ensem<T3>, Op> {
  typedef Ensem<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};


// Default trinary(Ensem,EScalar,Ensem) -> Ensem
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<Ensem<T1>, EScalar<T2>, Ensem<T3>, Op> {
  typedef Ensem<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};

// Default trinary(Ensem,Ensem,EScalar) -> Ensem
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<Ensem<T1>, Ensem<T2>, EScalar<T3>, Op> {
  typedef Ensem<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};

// Default trinary(EScalar,Ensem,Ensem) -> Ensem
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<EScalar<T1>, Ensem<T2>, Ensem<T3>, Op> {
  typedef Ensem<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};


// Default trinary(EScalar,EScalar,Ensem) -> Ensem
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<EScalar<T1>, EScalar<T2>, Ensem<T3>, Op> {
  typedef Ensem<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};

// Default trinary(OSscalar,Ensem,EScalar) -> Ensem
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<EScalar<T1>, Ensem<T2>, EScalar<T3>, Op> {
  typedef Ensem<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};

// Default trinary(Ensem,EScalar,EScalar) -> Ensem
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<Ensem<T1>, EScalar<T2>, EScalar<T3>, Op> {
  typedef Ensem<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};



// Specific EScalar cases
// Local operations

template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpAddAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpAddAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpSubtractAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpSubtractAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpMultiplyAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpMultiplyAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpDivideAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpDivideAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpModAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpModAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpBitwiseOrAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpBitwiseOrAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpBitwiseAndAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpBitwiseAndAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpBitwiseXorAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpBitwiseXorAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpLeftShiftAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpLeftShiftAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpRightShiftAssign> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpRightShiftAssign>::Type_t>  &Type_t;
};
 
// Specific Ensem cases
// Local operations
template<class T>
struct UnaryReturn<Ensem<T>, OpNot> {
  typedef Ensem<typename UnaryReturn<T, OpNot>::Type_t>  Type_t;
};


template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpAddAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpAddAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpSubtractAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpSubtractAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpMultiplyAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpMultiplyAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpDivideAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpDivideAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpModAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpModAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseOrAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseOrAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseAndAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseAndAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseXorAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseXorAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpLeftShiftAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLeftShiftAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpRightShiftAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpRightShiftAssign>::Type_t>  &Type_t;
};
 

// Mixed Ensem & EScalar cases
// Global operations
// Local operations
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpAddAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpAddAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpSubtractAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpSubtractAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpMultiplyAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpMultiplyAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpDivideAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpDivideAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpModAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpModAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseOrAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseOrAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseAndAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseAndAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseXorAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseXorAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpLeftShiftAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLeftShiftAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpRightShiftAssign> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpRightShiftAssign>::Type_t>  &Type_t;
};




//-----------------------------------------------------------------------------
// Ensem scalar operations
//-----------------------------------------------------------------------------

/*! \addtogroup escalar */
/*! @{ */

// Ensem scalar

// ! EScalar
template<class T>
struct UnaryReturn<EScalar<T>, OpNot> {
  typedef EScalar<typename UnaryReturn<T, OpNot>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<EScalar<T1>, OpNot>::Type_t
operator!(const EScalar<T1>& l)
{
  return ! l.elem();
}


// +EScalar
template<class T1>
inline typename UnaryReturn<EScalar<T1>, OpUnaryPlus>::Type_t
operator+(const EScalar<T1>& l)
{
  return +l.elem();
}


// -EScalar
template<class T1>
inline typename UnaryReturn<EScalar<T1>, OpUnaryMinus>::Type_t
operator-(const EScalar<T1>& l)
{
  return -l.elem();
}


// EScalar + EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpAdd>::Type_t
operator+(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem()+r.elem();
}


// EScalar - EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpSubtract>::Type_t
operator-(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() - r.elem();
}


// EScalar * EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpMultiply>::Type_t
operator*(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() * r.elem();
}

// Optimized  adj(EScalar)*EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpAdjMultiply>::Type_t
adjMultiply(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return adjMultiply(l.elem(), r.elem());
}

// Optimized  EScalar*adj(EScalar)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return multiplyAdj(l.elem(), r.elem());
}

// Optimized  adj(EScalar)*adj(EScalar)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return adjMultiplyAdj(l.elem(), r.elem());
}


// EScalar / EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpDivide>::Type_t
operator/(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() / r.elem();
}


// EScalar << EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpLeftShift> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpLeftShift>::Type_t
operator<<(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() << r.elem();
}


// EScalar >> EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpRightShift> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpRightShift>::Type_t>  Type_t;
};
 

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpRightShift>::Type_t
operator>>(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() >> r.elem();
}


// EScalar % EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpMod>::Type_t
operator%(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() % r.elem();
}


// EScalar ^ EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpBitwiseXor>::Type_t
operator^(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() ^ r.elem();
}


// EScalar & EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpBitwiseAnd>::Type_t
operator&(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() & r.elem();
}


// EScalar | EScalar
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpBitwiseOr>::Type_t
operator|(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() | r.elem();
}



// Comparisons

// EScalar < EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpLT> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpLT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpLT>::Type_t
operator<(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() < r.elem();
}


// EScalar <= EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpLE> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpLE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpLE>::Type_t
operator<=(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() <= r.elem();
}


// EScalar > EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpGT> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpGT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpGT>::Type_t
operator>(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() > r.elem();
}


// EScalar >= EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpGE> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpGE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpGE>::Type_t
operator>=(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() >= r.elem();
}


// EScalar == EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpEQ> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpEQ>::Type_t
operator==(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() == r.elem();
}


// EScalar != EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpNE> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpNE>::Type_t
operator!=(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() != r.elem();
}


// EScalar && EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpAnd> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpAnd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpAnd>::Type_t
operator&&(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() && r.elem();
}


// EScalar || EScalar
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, OpOr> {
  typedef EScalar<typename BinaryReturn<T1, T2, OpOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, OpOr>::Type_t
operator||(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return l.elem() || r.elem();
}



//-----------------------------------------------------------------------------
// Functions

// EScalar = adj(EScalar)
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnAdjoint>::Type_t
adj(const EScalar<T1>& s1)
{
  return adj(s1.elem());
}


// EScalar = conj(EScalar)
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnConjugate>::Type_t
conj(const EScalar<T1>& s1)
{
  return conj(s1.elem());
}


// EScalar = transpose(EScalar)
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnTranspose>::Type_t
transpose(const EScalar<T1>& s1)
{
  return transpose(s1.elem());
}


// EScalar = Trace(EScalar)
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnTrace>::Type_t
trace(const EScalar<T1>& s1)
{
  return trace(s1.elem());
}


// EScalar = Re(Trace(EScalar))
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnRealTrace>::Type_t
trace_real(const EScalar<T1>& s1)
{
  return trace_real(s1.elem());
}


// EScalar = Im(Trace(EScalar))
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnImagTrace>::Type_t
trace_imag(const EScalar<T1>& s1)
{
  return trace_imag(s1.elem());
}


// EScalar = Re(EScalar)
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnReal>::Type_t
real(const EScalar<T1>& s1)
{
  return real(s1.elem());
}


// EScalar = Im(EScalar)
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnImag>::Type_t
imag(const EScalar<T1>& s1)
{
  return imag(s1.elem());
}


// ArcCos
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnArcCos>::Type_t
acos(const EScalar<T1>& s1)
{
  return acos(s1.elem());
}


// ArcSin
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnArcSin>::Type_t
asin(const EScalar<T1>& s1)
{
  return asin(s1.elem());
}


// ArcTan
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnArcTan>::Type_t
atan(const EScalar<T1>& s1)
{
  return atan(s1.elem());
}


// Cos
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnCos>::Type_t
cos(const EScalar<T1>& s1)
{
  return cos(s1.elem());
}


// Exp
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnExp>::Type_t
exp(const EScalar<T1>& s1)
{
  return exp(s1.elem());
}


// Fabs
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnFabs>::Type_t
fabs(const EScalar<T1>& s1)
{
  return fabs(s1.elem());
}


// Log
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnLog>::Type_t
log(const EScalar<T1>& s1)
{
  return log(s1.elem());
}


// Sin
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnSin>::Type_t
sin(const EScalar<T1>& s1)
{
  return sin(s1.elem());
}


// Sqrt
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnSqrt>::Type_t
sqrt(const EScalar<T1>& s1)
{
  return sqrt(s1.elem());
}


// Tan
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnTan>::Type_t
tan(const EScalar<T1>& s1)
{
  return tan(s1.elem());
}


//! EScalar = pow(EScalar, EScalar)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, FnPow>::Type_t
pow(const EScalar<T1>& s1, const EScalar<T2>& s2)
{
  return pow(s1.elem(), s2.elem());
}

//! EScalar = pow(EScalar, int)
template<class T1>
inline EScalar<T1>
pow(const EScalar<T1>& s1, int s2)
{
  EScalar<T1>  d = s1;

  if (s2 > 0)
  {
    for(int j=1; j < s2; ++j)
    {
//	printf("pow: j=%d, s2=%d\n", j,s2);
      d.elem() *= s1.elem();
    }
  }
  else
  {
    std::cerr << __func__ << ": pow - illegal size\n";
    exit(1);
  }

  return d;
}


//! EScalar = atan2(EScalar, EScalar)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, FnArcTan2>::Type_t
atan2(const EScalar<T1>& s1, const EScalar<T2>& s2)
{
  return atan2(s1.elem(), s2.elem());
}


//! EScalar = cmplx(EScalar , EScalar)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, FnCmplx>::Type_t
cmplx(const EScalar<T1>& s1, const EScalar<T2>& s2)
{
  typename BinaryReturn<EScalar<T1>, EScalar<T2>, FnCmplx>::Type_t  d;

  return cmplx(s1.elem(), s2.elem());
}


//! EScalar = i * EScalar
template<class T>
inline typename UnaryReturn<EScalar<T>, FnTimesI>::Type_t
timesI(const EScalar<T>& s1)
{
  return timesI(s1.elem());
}

//! EScalar = -i * EScalar
template<class T>
inline typename UnaryReturn<EScalar<T>, FnTimesMinusI>::Type_t
timesMinusI(const EScalar<T>& s1)
{
  return timesMinusI(s1.elem());
}

//! bool = isZero(EScalar)
template<class T>
bool
isZero(const EScalar<T>& s1)
{
  return isZero(s1.elem());
}

//! bool = isNaN(EScalar)
template<class T>
bool
isNaN(const EScalar<T>& s1)
{
  return isNaN(s1.elem());
}

//! bool = isInf(EScalar)
template<class T>
bool
isInf(const EScalar<T>& s1)
{
  return isInf(s1.elem());
}

//! bool = isFinite(EScalar)
template<class T>
bool
isFinite(const EScalar<T>& s1)
{
  return isFinite(s1.elem());
}

//! dest [float type] = source [seed type]
template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnSeedToFloat>::Type_t
seedToFloat(const EScalar<T1>& s1)
{
  return seedToFloat(s1.elem());
}


//! EScalar = outerProduct(EScalar, EScalar)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const EScalar<T1>& l, const EScalar<T2>& r)
{
  return outerProduct(l.elem(),r.elem());
}


//! Contract over a specific list of indices
template<class T>
struct UnaryReturn<EScalar<T>, FnContract> {
  typedef EScalar<typename UnaryReturn<T, FnContract>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<EScalar<T>, FnContract>::Type_t
contract(const EScalar<T>& s1, const Array<int>& nn)
{
  return contract(s1.elem(),nn);
}


// EScalar = shift(EScalar, int)
template<class T1>
struct UnaryReturn<EScalar<T1>, FnShift> {
  typedef EScalar<typename UnaryReturn<T1, FnShift>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnShift>::Type_t
shift(const EScalar<T1>& s1, int sh)
{
  return shift(s1.elem(), sh);
}


// EScalar = cshift(EScalar, int)
template<class T1>
struct UnaryReturn<EScalar<T1>, FnCshift> {
  typedef EScalar<typename UnaryReturn<T1, FnCshift>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnCshift>::Type_t
cshift(const EScalar<T1>& s1, int sh)
{
  return cshift(s1.elem(), sh);
}


//-----------------------------------------------------------------------------
//! ENSEM Int to int primitive in conversion routine
template<class T> 
inline int 
toInt(const EScalar<T>& s) 
{
  return toInt(s.elem());
}

//! ENSEM Real to float primitive in conversion routine
template<class T> 
inline float
toFloat(const EScalar<T>& s) 
{
  return toFloat(s.elem());
}

//! ENSEM Double to double primitive in conversion routine
template<class T> 
inline double
toDouble(const EScalar<T>& s) 
{
  return toDouble(s.elem());
}

//! ENSEM Boolean to bool primitive in conversion routine
template<class T> 
inline bool
toBool(const EScalar<T>& s) 
{
  return toBool(s.elem());
}


//------------------------------------------
//! dest [float type] = source [int type]
template<class T, class T1>
inline
void cast_rep(T& d, const EScalar<T1>& s1)
{
  cast_rep(d, s1.elem());
}


//! dest [float type] = source [int type]
template<class T, class T1>
inline
void recast_rep(EScalar<T>& d, const EScalar<T1>& s1)
{
  cast_rep(d.elem(), s1.elem());
}


//------------------------------------------
// Global sum over site indices only
template<class T>
struct UnaryReturn<EScalar<T>, FnSum > {
  typedef EScalar<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<EScalar<T>, FnSum>::Type_t
sum(const EScalar<T>& s1)
{
//  return sum(s1.elem());
  return s1.elem();
}


// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T>
struct UnaryReturn<EScalar<T>, FnNorm2> {
  typedef EScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<EScalar<T>, FnLocalNorm2> {
  typedef EScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<EScalar<T>, FnLocalNorm2>::Type_t
localNorm2(const EScalar<T>& s1)
{
  return localNorm2(s1.elem());
}



//! EScalar = InnerProduct(adj(EScalar)*EScalar)
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, FnInnerProduct> {
  typedef EScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, EScalar<T2>, FnLocalInnerProduct> {
  typedef EScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, EScalar<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const EScalar<T1>& s1, const EScalar<T2>& s2)
{
  return localInnerProduct(s1.elem(), s2.elem());
}


//! EScalar = where(EScalar, EScalar, EScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<EScalar<T1>, EScalar<T2>, EScalar<T3>, FnWhere> {
  typedef EScalar<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<EScalar<T1>, EScalar<T2>, EScalar<T3>, FnWhere>::Type_t
where(const EScalar<T1>& a, const EScalar<T2>& b, const EScalar<T3>& c)
{
  return where(a.elem(), b.elem(), c.elem());
}



//-----------------------------------------------------------------------------
//! Extract components
template<class T1>
struct UnaryReturn<EScalar<T1>, FnPeekEnsem> {
  typedef EScalar<typename UnaryReturn<T1, FnPeekEnsem>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnPeekEnsem>::Type_t
peekEnsem(const EScalar<T1>& l, int nbin)
{
  return l.elem();
}


//! Extract components
template<class T1>
struct UnaryReturn<EScalar<T1>, FnPeekObsVector> {
  typedef EScalar<typename UnaryReturn<T1, FnPeekObsVector>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnPeekObsVector>::Type_t
peekObs(const EScalar<T1>& l, int row)
{
  return peekObs(l.elem(), row);
}


template<class T1>
struct UnaryReturn<EScalar<T1>, FnPeekObsTensor> {
  typedef EScalar<typename UnaryReturn<T1, FnPeekObsTensor>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnPeekObsTensor>::Type_t
peekObs(const EScalar<T1>& l, const Array<int>& nn)
{
  return peekObs(l.elem(), nn);
}


//! Extract components
template<class T1>
struct UnaryReturn<EScalar<T1>, FnPeekSpinVector> {
  typedef EScalar<typename UnaryReturn<T1, FnPeekSpinVector>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnPeekSpinVector>::Type_t
peekSpin(const EScalar<T1>& l, int row)
{
  return peekSpin(l.elem(), row);
}


template<class T1>
struct UnaryReturn<EScalar<T1>, FnPeekSpinMatrix> {
  typedef EScalar<typename UnaryReturn<T1, FnPeekSpinMatrix>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<EScalar<T1>, FnPeekSpinMatrix>::Type_t
peekSpin(const EScalar<T1>& l, int row, int col)
{
  return peekSpin(l.elem(), row, col);
}


//-----------------------------------------------------------------------------
//! Insert components
template<class T1, class T2>
inline EScalar<T1>&
pokeEnsem(EScalar<T1>& l, const EScalar<T2>& r, int nbin)
{
  l.elem() = r.elem();
  return l;
}


//! Insert components
template<class T1, class T2>
inline EScalar<T1>&
pokeObs(EScalar<T1>& l, const EScalar<T2>& r, int row)
{
  pokeObs(l.elem(), r.elem(), row);
  return l;
}


template<class T1, class T2>
inline EScalar<T1>&
pokeObs(EScalar<T1>& l, const EScalar<T2>& r, const Array<int>& nn)
{
  pokeObs(l.elem(), r.elem(), nn);
  return l;
}

//! Insert spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline EScalar<T1>&
pokeSpin(EScalar<T1>& l, const EScalar<T2>& r, int row)
{
  pokeSpin(l.elem(),r.elem(),row);
  return l;
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline EScalar<T1>&
pokeSpin(EScalar<T1>& l, const EScalar<T2>& r, int row, int col)
{
  pokeSpin(l.elem(),r.elem(),row,col);
  return l;
}



//-----------------------------------------------------------------------------
// Broadcast operations
//! dest = 0
template<class T> 
inline
void zero_rep(EScalar<T>& d) 
{
  zero_rep(d.elem());
}


//-----------------------------------------------------------------------------
// Random numbers
//! dest  = random  
template<class T>
inline void
fill_random(EScalar<T>& d)
{
  fill_random(d.elem());
}


//! dest  = gaussian  
template<class T>
inline void
fill_gaussian(EScalar<T>& d, EScalar<T>& r1, EScalar<T>& r2)
{
  fill_gaussian(d.elem(), r1.elem(), r2.elem());
}


//-----------------------------------------------------------------------------
// Gamma operations
//-----------------------------------------------------------------------------
template<int N, int m, class T, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, EScalar<T>, OpGammaConstMultiply> {
  typedef EScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! GammaConst * EScalar
template<class T2,int N,int m>
inline typename BinaryReturn<GammaConst<N,m>, EScalar<T2>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<N,m> & l,const EScalar<T2> & r)
{
  return l * r.elem();
}

template<class T, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<EScalar<T>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef EScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! EScalar * GammaConst
template<class T1,int N,int m>
inline typename BinaryReturn<EScalar<T1>, GammaConst<N,m>, OpMultiplyGammaConst>::Type_t
operator*(const EScalar<T1> & l,const GammaConst<N,m> & r)
{
  return l.elem() * r;
}


template<class T, int N, class OpGammaTypeMultiply>
struct BinaryReturn<GammaType<N>, EScalar<T>, OpGammaTypeMultiply> {
  typedef EScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! GammaType * EScalar
template<class T,int N>
inline typename BinaryReturn<GammaType<N>, EScalar<T>, OpGammaTypeMultiply>::Type_t
operator*(const GammaType<N> & l,const EScalar<T> & r)
{
  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_gamma0_arg<T,N>, 
		 &Multiply_gamma1_arg<T,N>,
		 &Multiply_gamma2_arg<T,N>,
		 &Multiply_gamma3_arg<T,N>,
		 &Multiply_gamma4_arg<T,N>,
		 &Multiply_gamma5_arg<T,N>,
		 &Multiply_gamma6_arg<T,N>,
		 &Multiply_gamma7_arg<T,N>,
		 &Multiply_gamma8_arg<T,N>,
		 &Multiply_gamma9_arg<T,N>,
		 &Multiply_gamma10_arg<T,N>,
		 &Multiply_gamma11_arg<T,N>,
		 &Multiply_gamma12_arg<T,N>,
		 &Multiply_gamma13_arg<T,N>,
		 &Multiply_gamma14_arg<T,N>,
		 &Multiply_gamma15_arg<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).
  return s[l.elem()](r.elem());
}

template<class T, int N, class OpMultiplyGammaType>
struct BinaryReturn<EScalar<T>, GammaType<N>, OpMultiplyGammaType> {
  typedef EScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! EScalar * GammaType
template<class T,int N>
inline typename BinaryReturn<EScalar<T>, GammaType<N>, OpMultiplyGammaType>::Type_t
operator*(const EScalar<T> & l,const GammaType<N> & r)
{
  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_arg_gamma0<T,N>, 
		 &Multiply_arg_gamma1<T,N>,
		 &Multiply_arg_gamma2<T,N>,
		 &Multiply_arg_gamma3<T,N>,
		 &Multiply_arg_gamma4<T,N>,
		 &Multiply_arg_gamma5<T,N>,
		 &Multiply_arg_gamma6<T,N>,
		 &Multiply_arg_gamma7<T,N>,
		 &Multiply_arg_gamma8<T,N>,
		 &Multiply_arg_gamma9<T,N>,
		 &Multiply_arg_gamma10<T,N>,
		 &Multiply_arg_gamma11<T,N>,
		 &Multiply_arg_gamma12<T,N>,
		 &Multiply_arg_gamma13<T,N>,
		 &Multiply_arg_gamma14<T,N>,
		 &Multiply_arg_gamma15<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).
  return s[r.elem()](l.elem());
}


//-----------------------------------------------------------------------------
// GammaDP operations
//-----------------------------------------------------------------------------
template<int N, int m, class T, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, EScalar<T>, OpGammaConstDPMultiply> {
  typedef EScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! GammaConstDP * EScalar
template<class T2,int N,int m>
inline typename BinaryReturn<GammaConstDP<N,m>, EScalar<T2>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<N,m> & l,const EScalar<T2> & r)
{
  return l * r.elem();
}

template<class T, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<EScalar<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef EScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! EScalar * GammaConstDP
template<class T1,int N,int m>
inline typename BinaryReturn<EScalar<T1>, GammaConstDP<N,m>, OpMultiplyGammaConstDP>::Type_t
operator*(const EScalar<T1> & l,const GammaConstDP<N,m> & r)
{
  return l.elem() * r;
}


template<class T, int N, class OpGammaTypeDPMultiply>
struct BinaryReturn<GammaTypeDP<N>, EScalar<T>, OpGammaTypeDPMultiply> {
  typedef EScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! GammaTypeDP * EScalar
template<class T,int N>
inline typename BinaryReturn<GammaTypeDP<N>, EScalar<T>, OpGammaTypeDPMultiply>::Type_t
operator*(const GammaTypeDP<N> & l,const EScalar<T> & r)
{
  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_gamma0_DP_arg<T,N>, 
		 &Multiply_gamma1_DP_arg<T,N>,
		 &Multiply_gamma2_DP_arg<T,N>,
		 &Multiply_gamma3_DP_arg<T,N>,
		 &Multiply_gamma4_DP_arg<T,N>,
		 &Multiply_gamma5_DP_arg<T,N>,
		 &Multiply_gamma6_DP_arg<T,N>,
		 &Multiply_gamma7_DP_arg<T,N>,
		 &Multiply_gamma8_DP_arg<T,N>,
		 &Multiply_gamma9_DP_arg<T,N>,
		 &Multiply_gamma10_DP_arg<T,N>,
		 &Multiply_gamma11_DP_arg<T,N>,
		 &Multiply_gamma12_DP_arg<T,N>,
		 &Multiply_gamma13_DP_arg<T,N>,
		 &Multiply_gamma14_DP_arg<T,N>,
		 &Multiply_gamma15_DP_arg<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).
  return s[l.elem()](r.elem());
}

template<class T, int N, class OpMultiplyGammaTypeDP>
struct BinaryReturn<EScalar<T>, GammaTypeDP<N>, OpMultiplyGammaTypeDP> {
  typedef EScalar<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! EScalar * GammaTypeDP
template<class T,int N>
inline typename BinaryReturn<EScalar<T>, GammaTypeDP<N>, OpMultiplyGammaTypeDP>::Type_t
operator*(const EScalar<T> & l,const GammaTypeDP<N> & r)
{
  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_arg_gamma0_DP<T,N>, 
		 &Multiply_arg_gamma1_DP<T,N>,
		 &Multiply_arg_gamma2_DP<T,N>,
		 &Multiply_arg_gamma3_DP<T,N>,
		 &Multiply_arg_gamma4_DP<T,N>,
		 &Multiply_arg_gamma5_DP<T,N>,
		 &Multiply_arg_gamma6_DP<T,N>,
		 &Multiply_arg_gamma7_DP<T,N>,
		 &Multiply_arg_gamma8_DP<T,N>,
		 &Multiply_arg_gamma9_DP<T,N>,
		 &Multiply_arg_gamma10_DP<T,N>,
		 &Multiply_arg_gamma11_DP<T,N>,
		 &Multiply_arg_gamma12_DP<T,N>,
		 &Multiply_arg_gamma13_DP<T,N>,
		 &Multiply_arg_gamma14_DP<T,N>,
		 &Multiply_arg_gamma15_DP<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).
  return s[r.elem()](l.elem());
}


/*! @} */  // end of group escalar

//-----------------------------------------------------------------------------
// Ensemble operations
//-----------------------------------------------------------------------------

/*! \addtogroup eensem */
/*! @{ */

// Ensem operations


//---------------------------------------------------------
//! Force the ensemble to have certain statistics
//! Return the proper rescaling factor
inline
double rescaleFactor(int N, EnsemType_t ensem_type)
{
  /* Rescaling factor */
  double factor;

  switch (ensem_type)
  {
  case ENSEM_JACKKNIFE:
    factor = -(N - 1);
    break;

  case ENSEM_BOOTSTRAP:
    factor = sqrt(double(N - 1));
    break;

  default:
    std::cerr << __func__ << ": unknown ensemble type = " << ensem_type;
	exit(1);
  }

  return factor;
}


//! Return a new ensemble that is rescaled
template<class T>
inline
Ensem<T> rescaleEnsem(const Ensem<T>& src, double factor)
{
  src.checkSize(__func__);

  const int num = src.size();

  Ensem<T> dst;
  dst.resize(num);
  dst.setEnsemType(src.getEnsemType());

  T  avg;
  avg.resize(src.elem(0));   // resize according to src
  avg = zero;

  for(int i=0; i < num; ++i)
    avg += src.elem(i);

  typedef typename InternalScalar<T>::Type_t  Scalar_t;
  avg /= Scalar_t(num);

  Scalar_t  ff(factor);
  for(int i=0; i < num; ++i)
  {
//    dst.elem(i) = avg + (src.elem(i) - avg)*ff;
    T foo = src.elem(i) - avg;
    foo *= ff;
    dst.elem(i) = avg + foo;
  }

  return dst;
}


//! Return a new ensemble that is rescaled
template<class T>
inline
Ensem<T> rescaleEnsemUp(const Ensem<T>& src)
{
  return rescaleEnsem(src, src.rescaleFactor());
}


//! Return a new ensemble that is rescaled
template<class T>
inline
Ensem<T> rescaleEnsemDown(const Ensem<T>& src)
{
  return rescaleEnsem(src, 1.0 / src.rescaleFactor());
}


//! Drop data from the beginning of the ensemble
template<class T>
inline
Ensem<T> dropEnsem(const Ensem<T>& src, int nd)
{
  src.checkSize(__func__);
  typedef typename InternalScalar<T>::Type_t  Scalar_t;

  Ensem<T> d;

  if (nd <= 0)
  {
    d = src;
  }
  else
  {
    int num    = src.size();
    int num_nd = num - nd;
    int nobs   = src.numElem();
    if (num_nd <= 0) num_nd = 0;
  
    d.setEnsemType(d.getEnsemType());
    d.resize(num_nd); 
    d.resizeObs(nobs);

    for(int i=nd,k=0; i < num; ++i, ++k)
      d.elem(k) = src.elem(i);
  }

  return d;
}



//! Bin the ensemble
template<class T>
inline
Ensem<T> binEnsem(const Ensem<T>& src, int nb)
{
  src.checkSize(__func__);
  typedef typename InternalScalar<T>::Type_t  Scalar_t;

  Ensem<T> d;

  if (nb <= 1)
  {
    d = src;
  }
  else
  {
    int num    = src.size();
    int num_nb = num / nb;
    int nobs   = src.numElem();
    if (num % nb != 0) 
      ++num_nb;
    
    d.setEnsemType(d.getEnsemType());
    d.resize(num_nb); 
    d.resizeObs(nobs);

    int i=0, k=0;
    while(k < num_nb)
    {
      T  avg;
      avg.resize(nobs);
      zero_rep(avg);

      int ib = 0, ncnt = 0;
      while((ib < nb) && (i < num))
      {
	avg += src.elem(i);
	++ib;
	++i;
	++ncnt;
      }

      d.elem(k) = avg / Scalar_t(ncnt);
      k++;
    }
  }

  return d;
}



//! Ensem = ! Ensem
template<class T1>
inline typename UnaryReturn<Ensem<T1>, OpNot>::Type_t
operator!(const Ensem<T1>& l)
{
  typename UnaryReturn<Ensem<T1>, OpNot>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = ! l.elem(i);
  return d;
}


//! Ensem = +Ensem
template<class T1>
inline typename UnaryReturn<Ensem<T1>, OpUnaryPlus>::Type_t
operator+(const Ensem<T1>& l)
{
  typename UnaryReturn<Ensem<T1>, OpUnaryPlus>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = +l.elem(i);
  return d;
}


//! Ensem = -Ensem
template<class T1>
inline typename UnaryReturn<Ensem<T1>, OpUnaryMinus>::Type_t
operator-(const Ensem<T1>& l)
{
  typename UnaryReturn<Ensem<T1>, OpUnaryMinus>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = -l.elem(i);
  return d;
}


// Ensem + Ensem
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpAdd>::Type_t
operator+(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) + r.elem(i);
  return d;
}

// Ensem + EScalar
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpAdd>::Type_t
operator+(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) + r.elem();
  return d;
}

// EScalar + Ensem
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpAdd>::Type_t
operator+(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpAdd>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() + r.elem(i);
  return d;
}


// Ensem - Ensem
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpSubtract>::Type_t
operator-(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) - r.elem(i);
  return d;
}

// Ensem - EScalar
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpSubtract>::Type_t
operator-(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) - r.elem();
  return d;
}

// EScalar - Ensem
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpSubtract>::Type_t
operator-(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpSubtract>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() - r.elem(i);
  return d;
}


// Ensem * Ensem
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpMultiply>::Type_t
operator*(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  Ensem<T1> ll(rescaleEnsemDown(l));
  Ensem<T2> rr(rescaleEnsemDown(r));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = ll.elem(i) * rr.elem(i);

  return rescaleEnsemUp(d);
}

// Ensem * EScalar
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpMultiply>::Type_t
operator*(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) * r.elem();
  return d;
}

// EScalar * Ensem
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpMultiply>::Type_t
operator*(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpMultiply>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() * r.elem(i);
  return d;
}


// Optimized  adj(Ensem)*Ensem
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpAdjMultiply>::Type_t
adjMultiply(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpAdjMultiply>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  Ensem<T1> ll(rescaleEnsemDown(l));
  Ensem<T2> rr(rescaleEnsemDown(r));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = adjMultiply(ll.elem(i), rr.elem(i));

  return rescaleEnsemUp(d);
}

// Optimized  adj(Ensem)*EScalar
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpAdjMultiply>::Type_t
adjMultiply(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpAdjMultiply>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = adjMultiply(l.elem(i), r.elem());
  return d;
}

// Optimized  adj(EScalar)*Ensem
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpAdjMultiply>::Type_t
adjMultiply(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpAdjMultiply>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = adjMultiply(l.elem(), r.elem(i));
  return d;
}


// Optimized  Ensem*adj(Ensem)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpMultiplyAdj>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  Ensem<T1> ll(rescaleEnsemDown(l));
  Ensem<T2> rr(rescaleEnsemDown(r));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = multiplyAdj(ll.elem(i), rr.elem(i));

  return rescaleEnsemUp(d);
}

// Optimized  Ensem*adj(EScalar)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpMultiplyAdj>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = multiplyAdj(l.elem(i), r.elem());
  return d;
}

// Optimized  EScalar*adj(Ensem)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpMultiplyAdj>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = multiplyAdj(l.elem(), r.elem(i));
  return d;
}


// Optimized  adj(Ensem)*adj(Ensem)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpAdjMultiplyAdj>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  Ensem<T1> ll(rescaleEnsemDown(l));
  Ensem<T2> rr(rescaleEnsemDown(r));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = adjMultiplyAdj(ll.elem(i), rr.elem(i));

  return rescaleEnsemUp(d);
}

// Optimized  adj(Ensem)*adj(EScalar)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpAdjMultiplyAdj>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = adjMultiplyAdj(l.elem(i), r.elem(i));
  return d;
}

// Optimized  adj(EScalar)*adj(Ensem)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpAdjMultiplyAdj>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = adjMultiplyAdj(l.elem(i), r.elem(i));
  return d;
}


// Ensem / Ensem
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpDivide>::Type_t
operator/(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpDivide>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  Ensem<T1> ll(rescaleEnsemDown(l));
  Ensem<T2> rr(rescaleEnsemDown(r));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = ll.elem(i) / rr.elem(i);

  return rescaleEnsemUp(d);
}

// Ensem / EScalar
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpDivide>::Type_t
operator/(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpDivide>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) / r.elem();
  return d;
}


// EScalar / Ensem
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpDivide>::Type_t
operator/(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpDivide>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  Ensem<T2> rr(rescaleEnsemDown(r));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() / rr.elem(i);

  return rescaleEnsemUp(d);
}


//! Ensem << Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpLeftShift> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpLeftShift>::Type_t
operator<<(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpLeftShift>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) << r.elem(i);
  return d;
}

//! Ensem << EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpLeftShift> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpLeftShift>::Type_t
operator<<(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpLeftShift>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) << r.elem();
  return d;
}

//! EScalar << Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpLeftShift> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpLeftShift>::Type_t
operator<<(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpLeftShift>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() << r.elem(i);
  return d;
}


//! Ensem >> Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpRightShift> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpRightShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpRightShift>::Type_t
operator>>(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpRightShift>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) >> r.elem(i);
  return d;
}

// Ensem >> EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpRightShift> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpRightShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpRightShift>::Type_t
operator>>(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpRightShift>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) >> r.elem();
  return d;
}

// EScalar >> Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpRightShift> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpRightShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpRightShift>::Type_t
operator>>(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpRightShift>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() >> r.elem(i);
  return d;
}


// Ensem % Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpMod> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpMod>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpMod>::Type_t
operator%(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpMod>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) % r.elem(i);
  return d;
}

// Ensem % EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpMod> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpMod>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpMod>::Type_t
operator%(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpMod>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) % r.elem();
  return d;
}

// EScalar % Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpMod> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpMod>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpMod>::Type_t
operator%(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpMod>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() % r.elem(i);
  return d;
}


// Ensem ^ Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseXor> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseXor>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseXor>::Type_t
operator^(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseXor>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) ^ r.elem(i);
  return d;
}

// Ensem ^ EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseXor> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseXor>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseXor>::Type_t
operator^(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseXor>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) ^ r.elem();
  return d;
}

// EScalar ^ Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseXor> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseXor>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseXor>::Type_t
operator^(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseXor>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() ^ r.elem(i);
  return d;
}


// Ensem & EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseAnd> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseAnd>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseAnd>::Type_t
operator&(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseAnd>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) & r.elem(i);
  return d;
}

// Ensem & EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseAnd> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseAnd>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseAnd>::Type_t
operator&(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseAnd>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) & r.elem();
  return d;
}

// EScalar & Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseAnd> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseAnd>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseAnd>::Type_t
operator&(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseAnd>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() & r.elem(i);
  return d;
}


// Ensem | EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseOr> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseOr>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseOr>::Type_t
operator|(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpBitwiseOr>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) | r.elem(i);
  return d;
}

// Ensem | EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseOr> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseOr>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseOr>::Type_t
operator|(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpBitwiseOr>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) | r.elem();
  return d;
}

// EScalar | Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseOr> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpBitwiseOr>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseOr>::Type_t
operator|(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpBitwiseOr>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() | r.elem(i);
  return d;
}


// Comparisons

// Ensem < Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpLT> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpLT>::Type_t
operator<(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpLT>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) < r.elem(i);
  return d;
}

// Ensem < EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpLT> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpLT>::Type_t
operator<(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpLT>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) < r.elem();
  return d;
}

// EScalar < Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpLT> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpLT>::Type_t
operator<(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpLT>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() < r.elem(i);
  return d;
}


// Ensem <= Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpLE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpLE>::Type_t
operator<=(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpLE>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) <= r.elem(i);
  return d;
}

// Ensem <= EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpLE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpLE>::Type_t
operator<=(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpLE>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) <= r.elem();
  return d;
}

// EScalar <= Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpLE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpLE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpLE>::Type_t
operator<=(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpLE>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() <= r.elem(i);
  return d;
}


// Ensem > Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpGT> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpGT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpGT>::Type_t
operator>(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) > r.elem(i);
  return d;
}

// Ensem > EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpGT> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpGT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpGT>::Type_t
operator>(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) > r.elem();
  return d;
}

// Ensem > Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpGT> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpGT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpGT>::Type_t
operator>(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() > r.elem(i);
  return d;
}



// Ensem >= Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpGE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpGE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpGE>::Type_t
operator>=(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) >= r.elem(i);
  return d;
}

// Ensem >= EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpGE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpGE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpGE>::Type_t
operator>=(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) >= r.elem();
  return d;
}

// EScalar >= Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpGE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpGE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpGE>::Type_t
operator>=(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() >= r.elem(i);
  return d;
}


// Ensem == Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpEQ> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpEQ>::Type_t
operator==(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) == r.elem(i);
  return d;
}

// Ensem == EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpEQ> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpEQ>::Type_t
operator==(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) == r.elem();
  return d;
}

// EScalar == Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpEQ> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpEQ>::Type_t
operator==(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpGT>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() == r.elem(i);
  return d;
}


// Ensem != Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpNE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpNE>::Type_t
operator!=(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpNE>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) != r.elem(i);
  return d;
}

// Ensem != EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpNE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpNE>::Type_t
operator!=(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpNE>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) != r.elem();
  return d;
}

// Ensem != Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpNE> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpNE>::Type_t
operator!=(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpNE>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() != r.elem(i);
  return d;
}


// Ensem && Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpAnd> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpAnd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpAnd>::Type_t
operator&&(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpAnd>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) && r.elem(i);
  return d;
}

// Ensem && EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpAnd> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpAnd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpAnd>::Type_t
operator&&(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpAnd>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) && r.elem();
  return d;
}

// EScalar && Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpAnd> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpAnd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpAnd>::Type_t
operator&&(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpAnd>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() && r.elem(i);
  return d;
}


// Ensem || Ensem
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, OpOr> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpOr>::Type_t
operator||(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, OpOr>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) || r.elem(i);
  return d;
}

// Ensem || EScalar
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, OpOr> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpOr>::Type_t
operator||(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, OpOr>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) || r.elem();
  return d;
}

// EScalar || Ensem
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, OpOr> {
  typedef Ensem<typename BinaryReturn<T1, T2, OpOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpOr>::Type_t
operator||(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, OpOr>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem() || r.elem(i);
  return d;
}


//-----------------------------------------------------------------------------
// Functions

// Ensem = adj(Ensem)
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnAdjoint>::Type_t
adj(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnAdjoint>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = adj(s1.elem(i));
  return d;
}


// Ensem = conj(Ensem)
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnConjugate>::Type_t
conj(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnConjugate>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = conj(s1.elem(i));
  return d;
}


// Ensem = transpose(Ensem)
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnTranspose>::Type_t
transpose(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnTranspose>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = transpose(s1.elem(i));
  return d;
}


// Ensem = Trace(Ensem)
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnTrace>::Type_t
trace(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnTrace>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = trace(s1.elem(i));
  return d;
}


// Ensem = Re(Trace(Ensem))
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnRealTrace>::Type_t
trace_real(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnRealTrace>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = trace_real(s1.elem(i));
  return d;
}


// Ensem = Im(Trace(Ensem))
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnImagTrace>::Type_t
trace_imag(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnImagTrace>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = trace_imag(s1.elem(i));
  return d;
}


// Ensem = Re(Ensem)
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnReal>::Type_t
real(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnReal>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = real(s1.elem(i));
  return d;
}


// Ensem = Im(Ensem)
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnImag>::Type_t
imag(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnImag>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = imag(s1.elem(i));
  return d;
}


// ArcCos
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnArcCos>::Type_t
acos(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnArcCos>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = acos(ss1.elem(i));

  return rescaleEnsemUp(d);
}

// ArcSin
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnArcSin>::Type_t
asin(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnArcSin>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = asin(ss1.elem(i));

  return rescaleEnsemUp(d);
}

// ArcTan
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnArcTan>::Type_t
atan(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnArcTan>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = atan(ss1.elem(i));

  return rescaleEnsemUp(d);
}

// Cos
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnCos>::Type_t
cos(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnCos>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = cos(s1.elem(i));
  return d;
}

// Exp
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnExp>::Type_t
exp(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnExp>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = exp(ss1.elem(i));

  return rescaleEnsemUp(d);
}

// Fabs
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnFabs>::Type_t
fabs(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnFabs>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = fabs(s1.elem(i));
  return d;
}

// Log
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnLog>::Type_t
log(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnLog>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = log(ss1.elem(i));

  return rescaleEnsemUp(d);
}

// Sin
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnSin>::Type_t
sin(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnSin>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = sin(ss1.elem(i));

  return rescaleEnsemUp(d);
}

// Sqrt
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnSqrt>::Type_t
sqrt(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnSqrt>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = sqrt(ss1.elem(i));

  return rescaleEnsemUp(d);
}

// Tan
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnTan>::Type_t
tan(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnTan>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = tan(ss1.elem(i));

  return rescaleEnsemUp(d);
}


//! Ensem = pow(Ensem, int)
template<class T1>
inline Ensem<T1>
pow(const Ensem<T1>& s1, int s2)
{
  Ensem<T1>  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  if (s2 > 0)
  {
    for(int i=0; i < d.size(); ++i)
    {
      d.elem(i) = ss1.elem(i);
      for(int j=1; j < s2; ++j)
      {
//	printf("pow: i=%d j=%d, s2=%d\n", i,j,s2);
	d.elem(i) *= ss1.elem(i);
      }
    }
  }
  else
  {
    std::cerr << __func__ << ": pow - illegal size\n";
    exit(1);
  }

  return rescaleEnsemUp(d);
}

//! Ensem = pow(Ensem, Ensem)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnPow>::Type_t
pow(const Ensem<T1>& s1, const Ensem<T2>& s2)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType(), s2.size(), s2.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));
  Ensem<T2> ss2(rescaleEnsemDown(s2));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = pow(ss1.elem(i), ss2.elem(i));

  return rescaleEnsemUp(d);
}

//! Ensem = pow(Ensem, EScalar)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnPow>::Type_t
pow(const Ensem<T1>& s1, const EScalar<T2>& s2)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = pow(ss1.elem(i), s2.elem());

  return rescaleEnsemUp(d);
}

//! Ensem = pow(EScalar, Ensem)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnPow>::Type_t
pow(const EScalar<T1>& s1, const Ensem<T2>& s2)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnPow>::Type_t  d;
  d.checkResize(__func__, s2.size(), s2.getEnsemType());

  Ensem<T2> ss2(rescaleEnsemDown(s2));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = pow(s1.elem(i), ss2.elem(i));

  return rescaleEnsemUp(d);
}


//! Ensem = atan2(Ensem, Ensem)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnArcTan2>::Type_t
atan2(const Ensem<T1>& s1, const Ensem<T2>& s2)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType(), s2.size(), s2.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));
  Ensem<T2> ss2(rescaleEnsemDown(s2));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = atan2(ss1.elem(i), ss2.elem(i));

  return rescaleEnsemUp(d);
}

//! Ensem = atan2(Ensem, EScalar)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnArcTan2>::Type_t
atan2(const Ensem<T1>& s1, const EScalar<T2>& s2)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  Ensem<T1> ss1(rescaleEnsemDown(s1));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = atan2(ss1.elem(i), s2.elem());

  return rescaleEnsemUp(d);
}

//! Ensem = atan2(EScalar, Ensem)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnArcTan2>::Type_t
atan2(const EScalar<T1>& s1, const Ensem<T2>& s2)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnArcTan2>::Type_t  d;
  d.checkResize(__func__, s2.size(), s2.getEnsemType());

  Ensem<T2> ss2(rescaleEnsemDown(s2));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = atan2(s1.elem(), ss2.elem(i));

  return rescaleEnsemUp(d);
}


//! Ensem = (Ensem , Ensem)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnCmplx>::Type_t
cmplx(const Ensem<T1>& s1, const Ensem<T2>& s2)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnCmplx>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType(), s2.size(), s2.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = cmplx(s1.elem(i), s2.elem(i));

  return d;
}

//! Ensem = (Ensem , EScalar)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnCmplx>::Type_t
cmplx(const Ensem<T1>& s1, const EScalar<T2>& s2)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnCmplx>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = cmplx(s1.elem(i), s2.elem());

  return d;
}

//! Ensem = (EScalar , Ensem)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnCmplx>::Type_t
cmplx(const EScalar<T1>& s1, const Ensem<T2>& s2)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnCmplx>::Type_t  d;
  d.checkResize(__func__, s2.size(), s2.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = cmplx(s1.elem(), s2.elem(i));

  return d;
}


//! Ensem = i * Ensem
template<class T>
inline typename UnaryReturn<Ensem<T>, FnTimesI>::Type_t
timesI(const Ensem<T>& s1)
{
  typename UnaryReturn<Ensem<T>, FnTimesI>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = timesI(s1.elem(i));

  return d;
}

//! Ensem = -i * Ensem
template<class T>
inline typename UnaryReturn<Ensem<T>, FnTimesMinusI>::Type_t
timesMinusI(const Ensem<T>& s1)
{
  typename UnaryReturn<Ensem<T>, FnTimesMinusI>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = timesMinusI(s1.elem(i));

  return d;
}


//! bool = isNaN(Ensem)
template<class T>
bool
isNaN(const Ensem<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.size(); ++i)
    d |= isNaN(s1.elem(i));

  return d;
}


//! bool = isZero(Ensem)
template<class T>
bool
isZero(const Ensem<T>& s1)
{
  bool d = true;

  for(int i=0; i < s1.size(); ++i)
    d &= isZero(s1.elem(i));

  return d;
}


//! bool = isInf(Ensem)
template<class T>
bool
isInf(const Ensem<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.size(); ++i)
    d |= isInf(s1.elem(i));

  return d;
}


//! bool = isFinite(Ensem)
template<class T>
bool
isFinite(const Ensem<T>& s1)
{
  bool d = false;

  for(int i=0; i < s1.size(); ++i)
    d |= isFinite(s1.elem(i));

  return d;
}


//! dest [float type] = source [seed type]
template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnSeedToFloat>::Type_t
seedToFloat(const Ensem<T1>& s1)
{
  typename UnaryReturn<Ensem<T1>, FnSeedToFloat>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = seedToFloat(s1.elem(i));
  return d;
}


//! Ensem = outerProduct(Ensem, Ensem)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnOuterProduct>::Type_t
outerProduct(const Ensem<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnOuterProduct>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType(), r.size(), r.getEnsemType());

  Ensem<T1> ll(rescaleEnsemDown(l));
  Ensem<T2> rr(rescaleEnsemDown(r));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = outerProduct(ll.elem(i), rr.elem(i));

  return rescaleEnsemUp(d);
}

//! Ensem = outerProduct(Ensem, EScalar)
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const Ensem<T1>& l, const EScalar<T2>& r)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnOuterProduct>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  Ensem<T1> ll(rescaleEnsemDown(l));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = outerProduct(ll.elem(i), r.elem());

  return rescaleEnsemUp(d);
}

//! Ensem = outerProduct(EScalar, Ensem)
template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnOuterProduct>::Type_t
outerProduct(const EScalar<T1>& l, const Ensem<T2>& r)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnOuterProduct>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  Ensem<T1> rr(rescaleEnsemDown(r));

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = outerProduct(l.elem(), rr.elem(i));

  return rescaleEnsemUp(d);
}

//! Contract over a specific list of indices
template<class T>
struct UnaryReturn<Ensem<T>, FnContract> {
  typedef Ensem<typename UnaryReturn<T, FnContract>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<Ensem<T>, FnContract>::Type_t
contract(const Ensem<T>& s1, const Array<int>& nn)
{
  typename UnaryReturn<Ensem<T>, FnContract>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = contract(s1.elem(i), nn);

  return d;
}


// Ensem = shift(Ensem, int)
template<class T1>
struct UnaryReturn<Ensem<T1>, FnShift> {
  typedef Ensem<typename UnaryReturn<T1, FnShift>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnShift>::Type_t
shift(const Ensem<T1>& s1, int sh)
{
  typename UnaryReturn<Ensem<T1>, FnShift>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = shift(s1.elem(i), sh);

  return d;
}


// Ensem = cshift(Ensem, int)
template<class T1>
struct UnaryReturn<Ensem<T1>, FnCshift> {
  typedef Ensem<typename UnaryReturn<T1, FnCshift>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnCshift>::Type_t
cshift(const Ensem<T1>& s1, int sh)
{
  typename UnaryReturn<Ensem<T1>, FnCshift>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = cshift(s1.elem(i), sh);

  return d;
}


//-----------------------------------------------------------------------------
//! Extract components
template<class T1>
struct UnaryReturn<Ensem<T1>, FnPeekEnsem> {
  typedef EScalar<typename UnaryReturn<T1, FnPeekEnsem>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnPeekEnsem>::Type_t
peekEnsem(const Ensem<T1>& l, int nbin)
{
  l.checkSize(__func__);
  return l.elem(nbin);
}


//! Extract components
template<class T1>
struct UnaryReturn<Ensem<T1>, FnPeekObsVector> {
  typedef Ensem<typename UnaryReturn<T1, FnPeekObsVector>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnPeekObsVector>::Type_t
peekObs(const Ensem<T1>& s1, int row)
{
  typename UnaryReturn<Ensem<T1>, FnPeekObsVector>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = peekObs(s1.elem(i), row);

  return d;
}

//! Extract components
template<class T1>
struct UnaryReturn<Ensem<T1>, FnPeekObsTensor> {
  typedef Ensem<typename UnaryReturn<T1, FnPeekObsTensor>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<Ensem<T1>, FnPeekObsTensor>::Type_t
peekObs(const Ensem<T1>& s1, const Array<int>& nn)
{
  typename UnaryReturn<Ensem<T1>, FnPeekObsTensor>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = peekObs(s1.elem(i), nn);

  return d;
}


//-----------------------------------------------------------------------------
//! Insert components
template<class T1>
inline Ensem<T1>&
pokeEnsem(Ensem<T1>& l, const EScalar<T1>& r, int nbin)
{
  l.checkSize(__func__);
  l.elem(nbin) = r.elem();
  return l;
}


//! Insert components
template<class T1, class T2>
inline Ensem<T1>&
pokeObs(Ensem<T1>& l, const Ensem<T2>& r, int row)
{
  l.checkSize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < l.size(); ++i)
    pokeObs(l.elem(i), r.elem(i), row);

  return l;
}

//! Insert components
template<class T1, class T2>
inline Ensem<T1>&
pokeObs(Ensem<T1>& l, const Ensem<T2>& r, const Array<int>& nn)
{
  l.checkSize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < l.size(); ++i)
    pokeObs(l.elem(i), r.elem(i), nn);

  return l;
}


//! Insert spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline Ensem<T1>&
pokeSpin(Ensem<T1>& l, const Ensem<T2>& r, int row)
{
  l.checkSize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < l.size(); ++i)
    pokeSpin(l.elem(i), r.elem(i), row);

  return l;
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline EScalar<T1>&
pokeSpin(Ensem<T1>& l, const Ensem<T2>& r, int row, int col)
{
  l.checkSize(__func__, r.size(), r.getEnsemType());

  for(int i=0; i < l.size(); ++i)
    pokeSpin(l.elem(i), r.elem(i), row, col);

  return l;
}

//-----------------------------------------------------------------------------
//! dest = 0
template<class T> 
inline void 
zero_rep(Ensem<T>& d) 
{
  d.checkSize(__func__);
  for(int i=0; i < d.size(); ++i)
    zero_rep(d.elem(i));
}


//! dest  = random  
template<class T>
inline void
fill_random(Ensem<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.size(); ++i)
    fill_random(d.elem(i));
}


// Global sum over site indices only
template<class T>
struct UnaryReturn<Ensem<T>, FnSum> {
  typedef EScalar<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<Ensem<T>, FnSum>::Type_t
sum(const Ensem<T>& s1)
{
  typename UnaryReturn<Ensem<T>, FnSum>::Type_t  d;
  s1.checkSize(__func__);

  d.elem() = s1.elem(0);
  for(int i=1; i < s1.size(); ++i)
    d.elem() += s1.elem(i);

  return d;
}


// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T>
struct UnaryReturn<Ensem<T>, FnNorm2> {
  typedef EScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<Ensem<T>, FnLocalNorm2> {
  typedef Ensem<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<Ensem<T>, FnLocalNorm2>::Type_t
localNorm2(const Ensem<T>& s1)
{
  typename UnaryReturn<Ensem<T>, FnLocalNorm2>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < s1.size(); ++i)
    d.elem(i) = localNorm2(s1.elem(i));

  return d;
}


//! EScalar = InnerProduct(adj(Ensem)*Ensem)
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, FnInnerProduct> {
  typedef EScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, Ensem<T2>, FnLocalInnerProduct> {
  typedef Ensem<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const Ensem<T1>& s1, const Ensem<T2>& s2)
{
  typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnLocalInnerProduct>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType(), s2.size(), s2.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = localInnerProduct(s1.elem(i), s2.elem(i));

  return d;
}

//! EScalar = InnerProduct(adj(Ensem)*EScalar)
template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, FnInnerProduct> {
  typedef EScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<Ensem<T1>, EScalar<T2>, FnLocalInnerProduct> {
  typedef Ensem<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const Ensem<T1>& s1, const EScalar<T2>& s2)
{
  typename BinaryReturn<Ensem<T1>, EScalar<T2>, FnLocalInnerProduct>::Type_t  d;
  d.checkResize(__func__, s1.size(), s1.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = localInnerProduct(s1.elem(i), s2.elem());

  return d;
}

//! EScalar = InnerProduct(adj(EScalar)*Ensem)
template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, FnInnerProduct> {
  typedef EScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<EScalar<T1>, Ensem<T2>, FnLocalInnerProduct> {
  typedef Ensem<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const EScalar<T1>& s1, const Ensem<T2>& s2)
{
  typename BinaryReturn<EScalar<T1>, Ensem<T2>, FnLocalInnerProduct>::Type_t  d;
  d.checkResize(__func__, s2.size(), s2.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = localInnerProduct(s1.elem(), s2.elem(i));

  return d;
}


//! Ensem = where(Ensem, Ensem, Ensem)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<Ensem<T1>, Ensem<T2>, Ensem<T3>, FnWhere> {
  typedef Ensem<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<Ensem<T1>, Ensem<T2>, Ensem<T3>, FnWhere>::Type_t
where(const Ensem<T1>& a, const Ensem<T2>& b, const Ensem<T3>& c)
{
  typename TrinaryReturn<Ensem<T1>, Ensem<T2>, Ensem<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, a.size(), a.getEnsemType(), b.size(), b.getEnsemType(), c.size(), c.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = where(a.elem(i), b.elem(i), c.elem(i));

  return d;
}

//! Ensem = where(Ensem, Ensem, EScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<Ensem<T1>, Ensem<T2>, EScalar<T3>, FnWhere> {
  typedef Ensem<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<Ensem<T1>, Ensem<T2>, EScalar<T3>, FnWhere>::Type_t
where(const Ensem<T1>& a, const Ensem<T2>& b, const EScalar<T3>& c)
{
  typename TrinaryReturn<Ensem<T1>, Ensem<T2>, EScalar<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, a.size(), a.getEnsemType(), b.size(), b.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = where(a.elem(i), b.elem(i), c.elem());

  return d;
}

//! Ensem = where(Ensem, EScalar, Ensem)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<Ensem<T1>, EScalar<T2>, Ensem<T3>, FnWhere> {
  typedef Ensem<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<Ensem<T1>, EScalar<T2>, Ensem<T3>, FnWhere>::Type_t
where(const Ensem<T1>& a, const EScalar<T2>& b, const Ensem<T3>& c)
{
  typename TrinaryReturn<Ensem<T1>, EScalar<T2>, Ensem<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, a.size(), a.getEnsemType(), c.size(), c.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = where(a.elem(i), b.elem(), c.elem(i));

  return d;
}

//! Ensem = where(Ensem, EScalar, EScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<Ensem<T1>, EScalar<T2>, EScalar<T3>, FnWhere> {
  typedef Ensem<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<Ensem<T1>, EScalar<T2>, EScalar<T3>, FnWhere>::Type_t
where(const Ensem<T1>& a, const EScalar<T2>& b, const EScalar<T3>& c)
{
  typename TrinaryReturn<Ensem<T1>, EScalar<T2>, EScalar<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, a.size(), a.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = where(a.elem(i), b.elem(), c.elem());

  return d;
}

//! Ensem = where(EScalar, Ensem, Ensem)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<EScalar<T1>, Ensem<T2>, Ensem<T3>, FnWhere> {
  typedef Ensem<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<EScalar<T1>, Ensem<T2>, Ensem<T3>, FnWhere>::Type_t
where(const EScalar<T1>& a, const Ensem<T2>& b, const Ensem<T3>& c)
{
  typename TrinaryReturn<EScalar<T1>, Ensem<T2>, Ensem<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, b.size(), b.getEnsemType(), c.size(), c.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = where(a.elem(), b.elem(i), c.elem(i));

  return d;
}

//! Ensem = where(EScalar, Ensem, EScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<EScalar<T1>, Ensem<T2>, EScalar<T3>, FnWhere> {
  typedef Ensem<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<EScalar<T1>, Ensem<T2>, EScalar<T3>, FnWhere>::Type_t
where(const EScalar<T1>& a, const Ensem<T2>& b, const EScalar<T3>& c)
{
  typename TrinaryReturn<EScalar<T1>, Ensem<T2>, EScalar<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, b.size(), b.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = where(a.elem(), b.elem(i), c.elem());

  return d;
}

//! Ensem = where(EScalar, EScalar, Ensem)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<EScalar<T1>, EScalar<T2>, Ensem<T3>, FnWhere> {
  typedef Ensem<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<EScalar<T1>, EScalar<T2>, Ensem<T3>, FnWhere>::Type_t
where(const EScalar<T1>& a, const EScalar<T2>& b, const Ensem<T3>& c)
{
  typename TrinaryReturn<EScalar<T1>, EScalar<T2>, Ensem<T3>, FnWhere>::Type_t  d;
  d.checkResize(__func__, c.size(), c.getEnsemType());

  for(int i=0; i < d.size(); ++i)
    d.elem(i) = where(a.elem(), b.elem(), c.elem(i));

  return d;
}


//-----------------------------------------------------------------------------
// Gamma operations
//-----------------------------------------------------------------------------
template<int N, int m, class T, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, Ensem<T>, OpGammaConstMultiply> {
  typedef Ensem<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! GammaConst * Ensem
template<class T,int N,int m>
inline typename BinaryReturn<GammaConst<N,m>, Ensem<T>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<N,m> & l,const Ensem<T> & r)
{
  typename BinaryReturn<GammaConst<N,m>, Ensem<T>, OpGammaConstMultiply>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  // gamma mults only move around re/im - no need to rescale
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l * r.elem(i);

  return d;
}

template<class T, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<Ensem<T>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef Ensem<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! Ensem * GammaConst
template<class T,int N,int m>
inline typename BinaryReturn<Ensem<T>, GammaConst<N,m>, OpMultiplyGammaConst>::Type_t
operator*(const Ensem<T> & l,const GammaConst<N,m> & r)
{
  typename BinaryReturn<Ensem<T>, GammaConst<N,m>, OpMultiplyGammaConst>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  // gamma mults only move around re/im - no need to rescale
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) * r;

  return d;
}


template<class T, int N, class OpGammaTypeMultiply>
struct BinaryReturn<GammaType<N>, Ensem<T>, OpGammaTypeMultiply> {
  typedef Ensem<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! GammaType * Ensem
template<class T,int N>
inline typename BinaryReturn<GammaType<N>, Ensem<T>, OpGammaTypeMultiply>::Type_t
operator*(const GammaType<N> & l,const Ensem<T> & r)
{
  typename BinaryReturn<GammaType<N>, Ensem<T>, OpGammaTypeMultiply>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_gamma0_arg<T,N>, 
		 &Multiply_gamma1_arg<T,N>,
		 &Multiply_gamma2_arg<T,N>,
		 &Multiply_gamma3_arg<T,N>,
		 &Multiply_gamma4_arg<T,N>,
		 &Multiply_gamma5_arg<T,N>,
		 &Multiply_gamma6_arg<T,N>,
		 &Multiply_gamma7_arg<T,N>,
		 &Multiply_gamma8_arg<T,N>,
		 &Multiply_gamma9_arg<T,N>,
		 &Multiply_gamma10_arg<T,N>,
		 &Multiply_gamma11_arg<T,N>,
		 &Multiply_gamma12_arg<T,N>,
		 &Multiply_gamma13_arg<T,N>,
		 &Multiply_gamma14_arg<T,N>,
		 &Multiply_gamma15_arg<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).

  // gamma mults only move around re/im - no need to rescale
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = s[l.elem()](r.elem(i));

  return d;
}

template<class T, int N, class OpMultiplyGammaType>
struct BinaryReturn<Ensem<T>, GammaType<N>, OpMultiplyGammaType> {
  typedef Ensem<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! Ensem * GammaType
template<class T,int N>
inline typename BinaryReturn<Ensem<T>, GammaType<N>, OpMultiplyGammaType>::Type_t
operator*(const Ensem<T> & l,const GammaType<N> & r)
{
  typename BinaryReturn<Ensem<T>, GammaType<N>, OpMultiplyGammaType>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_arg_gamma0<T,N>, 
		 &Multiply_arg_gamma1<T,N>,
		 &Multiply_arg_gamma2<T,N>,
		 &Multiply_arg_gamma3<T,N>,
		 &Multiply_arg_gamma4<T,N>,
		 &Multiply_arg_gamma5<T,N>,
		 &Multiply_arg_gamma6<T,N>,
		 &Multiply_arg_gamma7<T,N>,
		 &Multiply_arg_gamma8<T,N>,
		 &Multiply_arg_gamma9<T,N>,
		 &Multiply_arg_gamma10<T,N>,
		 &Multiply_arg_gamma11<T,N>,
		 &Multiply_arg_gamma12<T,N>,
		 &Multiply_arg_gamma13<T,N>,
		 &Multiply_arg_gamma14<T,N>,
		 &Multiply_arg_gamma15<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).

  // gamma mults only move around re/im - no need to rescale
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = s[r.elem()](l.elem(i));

  return d;
}


//-----------------------------------------------------------------------------
// GammaDP operations
//-----------------------------------------------------------------------------
template<int N, int m, class T, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, Ensem<T>, OpGammaConstDPMultiply> {
  typedef Ensem<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! GammaConstDP * Ensem
template<class T,int N,int m>
inline typename BinaryReturn<GammaConstDP<N,m>, Ensem<T>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<N,m> & l,const Ensem<T> & r)
{
  typename BinaryReturn<GammaConstDP<N,m>, Ensem<T>, OpGammaConstDPMultiply>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  // gamma mults only move around re/im - no need to rescale
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l * r.elem(i);

  return d;
}

template<class T, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<Ensem<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef Ensem<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! Ensem * GammaConstDP
template<class T,int N,int m>
inline typename BinaryReturn<Ensem<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP>::Type_t
operator*(const Ensem<T> & l,const GammaConstDP<N,m> & r)
{
  typename BinaryReturn<Ensem<T>, GammaConstDP<N,m>, OpMultiplyGammaConstDP>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  // gamma mults only move around re/im - no need to rescale
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = l.elem(i) * r;

  return d;
}


template<class T, int N, class OpGammaTypeDPMultiply>
struct BinaryReturn<GammaTypeDP<N>, Ensem<T>, OpGammaTypeDPMultiply> {
  typedef Ensem<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! GammaTypeDP * Ensem
template<class T,int N>
inline typename BinaryReturn<GammaTypeDP<N>, Ensem<T>, OpGammaTypeDPMultiply>::Type_t
operator*(const GammaTypeDP<N> & l,const Ensem<T> & r)
{
  typename BinaryReturn<GammaTypeDP<N>, Ensem<T>, OpGammaTypeDPMultiply>::Type_t  d;
  d.checkResize(__func__, r.size(), r.getEnsemType());

  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_gamma0_DP_arg<T,N>, 
		 &Multiply_gamma1_DP_arg<T,N>,
		 &Multiply_gamma2_DP_arg<T,N>,
		 &Multiply_gamma3_DP_arg<T,N>,
		 &Multiply_gamma4_DP_arg<T,N>,
		 &Multiply_gamma5_DP_arg<T,N>,
		 &Multiply_gamma6_DP_arg<T,N>,
		 &Multiply_gamma7_DP_arg<T,N>,
		 &Multiply_gamma8_DP_arg<T,N>,
		 &Multiply_gamma9_DP_arg<T,N>,
		 &Multiply_gamma10_DP_arg<T,N>,
		 &Multiply_gamma11_DP_arg<T,N>,
		 &Multiply_gamma12_DP_arg<T,N>,
		 &Multiply_gamma13_DP_arg<T,N>,
		 &Multiply_gamma14_DP_arg<T,N>,
		 &Multiply_gamma15_DP_arg<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).

  // gamma mults only move around re/im - no need to rescale
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = s[l.elem()](r.elem(i));

  return d;
}

template<class T, int N, class OpMultiplyGammaTypeDP>
struct BinaryReturn<Ensem<T>, GammaTypeDP<N>, OpMultiplyGammaTypeDP> {
  typedef Ensem<typename UnaryReturn<T, OpUnaryPlus>::Type_t>  Type_t;
};

//! Ensem * GammaTypeDP
template<class T,int N>
inline typename BinaryReturn<Ensem<T>, GammaTypeDP<N>, OpMultiplyGammaTypeDP>::Type_t
operator*(const Ensem<T> & l,const GammaTypeDP<N> & r)
{
  typename BinaryReturn<Ensem<T>, GammaTypeDP<N>, OpMultiplyGammaTypeDP>::Type_t  d;
  d.checkResize(__func__, l.size(), l.getEnsemType());

  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_arg_gamma0_DP<T,N>, 
		 &Multiply_arg_gamma1_DP<T,N>,
		 &Multiply_arg_gamma2_DP<T,N>,
		 &Multiply_arg_gamma3_DP<T,N>,
		 &Multiply_arg_gamma4_DP<T,N>,
		 &Multiply_arg_gamma5_DP<T,N>,
		 &Multiply_arg_gamma6_DP<T,N>,
		 &Multiply_arg_gamma7_DP<T,N>,
		 &Multiply_arg_gamma8_DP<T,N>,
		 &Multiply_arg_gamma9_DP<T,N>,
		 &Multiply_arg_gamma10_DP<T,N>,
		 &Multiply_arg_gamma11_DP<T,N>,
		 &Multiply_arg_gamma12_DP<T,N>,
		 &Multiply_arg_gamma13_DP<T,N>,
		 &Multiply_arg_gamma14_DP<T,N>,
		 &Multiply_arg_gamma15_DP<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).

  // gamma mults only move around re/im - no need to rescale
  for(int i=0; i < d.size(); ++i)
    d.elem(i) = s[r.elem()](l.elem(i));

  return d;
}


/*! @} */  // end of group eensem

} // namespace ENSEM
