// -*- C++ -*-
// $Id: ensem_scalar_specific.h,v 2.0 2008/12/05 04:43:34 edwards Exp $
//
// ENSEM interface
//
// Routines specific to a scalar platform 

namespace ENSEM {

// Use separate defs here. This will cause subroutine calls under g++

//-------------------------------------------------------------------------------------
//! Antisymmetric tensor
/*! 
 * For now, instead of making a struct/tag for special ops like
 * Gamma matrices, make a function that constructs a real Tensor
 * object and simply manually do the contractions using contract
 *
 * This is not optimal and certainly wasteful, but its quick and
 * easy to get going.
 */
template<int N>
inline TensorReal
antiSymTensor()
{
  std::cerr << __func__ << "Not implemented" << std::endl;

  TensorReal d;
  return d;   // junk to make the compiler happy
}


//! 3rd-rank anti-symmetric tensor
template<>
inline TensorReal
antiSymTensor<3>()
{
  TensorReal d;

  const int N = 3;
  Array<int> nn(N);
  nn[0] = nn[1] = nn[2] = N;

  d.resizeObs(nn);
  d = zero;
  
  // d = \epsilon^{i,j,k}
  // Permutations: +(0,1,2)+(1,2,0)+(2,0,1)-(1,0,2)-(0,2,1)-(2,1,0)

  d.elem().elem(((0)*N+1)*N+2) = 1.0;
  d.elem().elem(((1)*N+2)*N+0) = 1.0;
  d.elem().elem(((2)*N+0)*N+1) = 1.0;

  d.elem().elem(((1)*N+0)*N+2) = -1.0;
  d.elem().elem(((0)*N+2)*N+1) = -1.0;
  d.elem().elem(((2)*N+1)*N+0) = -1.0;

  return d;
}


//! 4th-rank anti-symmetric tensor
template<>
inline TensorReal
antiSymTensor<4>()
{
  TensorReal d;

  const int N = 4;
  Array<int> nn(N);
  nn = N;

  d.resizeObs(nn);
  d = zero;
  
  // d = \epsilon^{i,j,k,l}
  // Permutations: +(0,1,2,3)-(1,2,3,0)+(2,3,0,1)-(3,0,1,2)   // 1
  //               -(1,0,2,3)+(0,2,3,1)-(2,3,1,0)+(3,1,0,2)   // 2
  //               -(2,1,0,3)+(1,0,3,2)-(0,3,2,1)+(3,2,1,0)   // 3
  //               -(3,1,2,0)+(1,2,0,3)-(2,0,3,1)+(0,3,1,2)   // 4
  //               -(0,2,1,3)+(2,1,3,0)-(1,3,0,2)+(3,0,2,1)   // 5
  //               -(0,1,3,2)+(1,3,2,0)-(3,2,0,1)+(2,0,1,3)   // 6

  // 1
  d.elem().elem((((0)*N+1)*N+2)*N+3) =  1.0;
  d.elem().elem((((1)*N+2)*N+3)*N+0) = -1.0;
  d.elem().elem((((2)*N+3)*N+0)*N+1) =  1.0;
  d.elem().elem((((3)*N+0)*N+1)*N+2) = -1.0;

  // 2
  d.elem().elem((((1)*N+0)*N+2)*N+3) = -1.0;
  d.elem().elem((((0)*N+2)*N+3)*N+1) =  1.0;
  d.elem().elem((((2)*N+3)*N+1)*N+0) = -1.0;
  d.elem().elem((((3)*N+1)*N+0)*N+2) =  1.0;

  // 3
  d.elem().elem((((2)*N+1)*N+0)*N+3) = -1.0;
  d.elem().elem((((1)*N+0)*N+3)*N+2) =  1.0;
  d.elem().elem((((0)*N+3)*N+2)*N+1) = -1.0;
  d.elem().elem((((3)*N+2)*N+1)*N+0) =  1.0;

  // 4
  d.elem().elem((((3)*N+1)*N+2)*N+0) = -1.0;
  d.elem().elem((((1)*N+2)*N+0)*N+3) =  1.0;
  d.elem().elem((((2)*N+0)*N+3)*N+1) = -1.0;
  d.elem().elem((((0)*N+3)*N+1)*N+2) =  1.0;

  // 5
  d.elem().elem((((0)*N+2)*N+1)*N+3) = -1.0;
  d.elem().elem((((2)*N+1)*N+3)*N+0) =  1.0;
  d.elem().elem((((1)*N+3)*N+0)*N+2) = -1.0;
  d.elem().elem((((3)*N+0)*N+2)*N+1) =  1.0;

  // 6
  d.elem().elem((((0)*N+1)*N+3)*N+2) = -1.0;
  d.elem().elem((((1)*N+3)*N+2)*N+0) =  1.0;
  d.elem().elem((((3)*N+2)*N+0)*N+1) = -1.0;
  d.elem().elem((((2)*N+0)*N+1)*N+3) =  1.0;

  return d;
}


//! Symmetric Kronecker-delta
/*! 
 *  A two index beast of size N
 */
template<int N>
inline TensorReal
symTensor()
{
  TensorReal d;
  Array<int> nn(2); 
  nn = N;
  d.resizeObs(nn);
  d = zero;

  for(int i=0; i < N; ++i)
    d.elem().elem((i)*N+i) = 1.0;

  return d;
}


//-----------------------------------------------------------------------------
//! Contract over 1 index
template<class T>
inline typename UnaryReturn<EScalar<T>, FnContract>::Type_t
contract(const EScalar<T>& s1, int n1)
{
  Array<int> nn(1);  nn[0] = n1;
  return contract(s1,nn);
}

//! Contract over 2 indices
template<class T>
inline typename UnaryReturn<EScalar<T>, FnContract>::Type_t
contract(const EScalar<T>& s1, int n1, int n2)
{
  Array<int> nn(2);  nn[0] = n1; nn[1] = n2;
  return contract(s1,nn);
}


//! Contract over 1 index
template<class T>
inline typename UnaryReturn<Ensem<T>, FnContract>::Type_t
contract(const Ensem<T>& s1, int n1)
{
  Array<int> nn(1);  nn[0] = n1;
  return contract(s1,nn);
}

//! Contract over 2 indices
template<class T>
inline typename UnaryReturn<Ensem<T>, FnContract>::Type_t
contract(const Ensem<T>& s1, int n1, int n2)
{
  Array<int> nn(2);  nn[0] = n1; nn[1] = n2;
  return contract(s1,nn);
}


//-----------------------------------------------------------------------------
//! dest  = random  
template<class T>
void 
random(EScalar<T>& d)
{
  fill_random(d.elem());
}


//! dest  = random
template<class T>
void 
random(Ensem<T>& d)
{
  d.checkSize(__func__);
  for(int i=0; i < d.size(); ++i) 
    fill_random(d.elem(i));
}


//! dest  = gaussian   under a subset
template<class T>
void gaussian(Ensem<T>& d)
{
  d.checkSize(__func__);
  Ensem<T>  r1, r2;
  r1.resize(d);
  r2.resize(d);

  random(r1);
  random(r2);

  for(int i=0; i < d.size(); ++i) 
    fill_gaussian(d.elem(i), r1.elem(i), r2.elem(i));
}



//-----------------------------------------------------------------------------
// Global sums
//! EScalar = norm2(trace(adj(source)*source))
/*!
 * return  num(trace(adj(source)*source))
 *
 * Sum over the ensemble
 * Allow a global sum that sums over all indices
 */
template<class T>
inline typename UnaryReturn<Ensem<T>, FnNorm2>::Type_t
norm2(const Ensem<T>& s1)
{
  return sum(localNorm2(s1));
}


//! EScalar = innerProduct(adj(source1)*source2)
/*!
 * return  sum(trace(adj(source1)*source2))
 *
 * Sum over the ensemble
 */
template<class T1, class T2>
inline typename BinaryReturn<Ensem<T1>, Ensem<T2>, FnInnerProduct>::Type_t
innerProduct(const Ensem<T1>& s1, const Ensem<T2>& s2)
{
  return sum(localInnerProduct(s1,s2));
}


//! EScalar = norm2(trace(adj(Array<source>)*Array<source>))
/*!
 * return  \sum_{Array} \sum_x(trace(adj(Array<source>)*Array<source>))
 *
 * Sum over the ensemble
 * Allow a global sum that sums over all indices
 */
template<class T>
inline typename UnaryReturn<EScalar<T>, FnNorm2>::Type_t
norm2(const Array< EScalar<T> >& s1)
{
  typename UnaryReturn<EScalar<T>, FnNorm2>::Type_t  d;

  // Possibly loop entered
  zero_rep(d.elem());

  for(int n=0; n < s1.size(); ++n)
    d.elem() += localNorm2(s1[n].elem());

  return d;
}

//! EScalar = norm2(Array<Ensem>) under an explicit subset
/*!
 * return  \sum_{Array} \sum_x(trace(adj(Array<source>)*Array<source>))
 *
 * Sum over the ensemble
 * Allow a global sum that sums over all indices
 */
template<class T>
inline typename UnaryReturn<Ensem<T>, FnNorm2>::Type_t
norm2(const Array< Ensem<T> >& s1)
{
  typename UnaryReturn<Ensem<T>, FnNorm2>::Type_t  d;

  // Possibly loop entered
  zero_rep(d.elem());

  for(int n=0; n < s1.size(); ++n)
  {
    const Ensem<T>& ss1 = s1[n];
    for(int i=0; i < ss1.size(); ++i) 
      d.elem() += localNorm2(ss1.elem(i));
  }

  return d;
}


//-----------------------------------------------------------------------------
// Statistics
//! Calculate mean, err and err/mean for data
template<class T>
struct UnaryReturn<Ensem<T>, FnMean> {
  typedef EScalar<typename UnaryReturn<T, FnMean>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<Ensem<T>, FnMean>::Type_t
mean(const Ensem<T>& src)
{
  src.checkSize(__func__);
  const int num = src.size();

  typename UnaryReturn<T, FnMean>::Type_t  avg;
  avg.resize(src.elem(0));  // resize according to srcb
  avg = zero;

  for(int i=0; i < num; ++i)
    avg += src.elem(i);

  typedef typename InternalScalar<T>::Type_t  Scalar_t;
  avg /= Scalar_t(num);

  return avg;
}


//! Calculate mean, err and err/mean for data
template<class T>
struct UnaryReturn<Ensem<T>, FnVariance> {
  typedef EScalar<typename UnaryReturn<T, FnVariance>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<Ensem<T>, FnVariance>::Type_t
variance(const Ensem<T>& src)
{
  src.checkSize(__func__);
  const int num = src.size();

  typename UnaryReturn<T, FnMean>::Type_t  avg;
  avg.resize(src.elem(0));   // resize according to src
  avg = zero;

  for(int i=0; i < num; ++i)
    avg += src.elem(i);

  typedef typename InternalScalar<T>::Type_t  Scalar_t;
  avg /= Scalar_t(num);

  T  diff;
  diff.resize(avg);

  typename UnaryReturn<T, FnVariance>::Type_t  var;
  var.resize(real(avg));
  var = zero;

  for(int i=0; i < num; ++i)
  {
    diff = (src.elem(i) - avg);
    var += real(conj(diff) * diff);
  }

  var /= Scalar_t((num - 1) * num);

  return var;
}


//! Do a whole calc
template<class T>
inline std::ostream&
calc(std::ostream& s, const Ensem<T>& src)
{
  src.checkSize(__func__);
  const int num = src.size();

  typename UnaryReturn<Ensem<T>, FnMean>::Type_t  avg = mean(src);
  typename UnaryReturn<Ensem<T>, FnVariance>::Type_t  var = sqrt(variance(src));
  int len = avg.numElem();

  for(int k = 0; k < len; ++k)
  {
    s << k << "\t" << peekObs(avg,k) << "  " << peekObs(var,k) << "\n";
  }

  return s;
}


//-----------------------------------------------------------------------------
// Input and output of various flavors that are architecture specific

// XML output
template<class T>  
ADATXML::XMLWriter& operator<<(ADATXML::XMLWriter& xml, const Ensem<T>& d)
{
  xml.openTag("Ensem");

  XMLWriterAPI::AttributeList alist;

  for(int i=0; i < d.size(); ++i) 
  {
    alist.clear();
    alist.push_back(XMLWriterAPI::Attribute("cfg", i));

    xml.openTag("elem", alist);
    xml << d.elem(i);
    xml.closeTag();
  }

  xml.closeTag(); // Ensem

  return xml;
}


//! XML output 
template<class T>
inline
void write(ADATXML::XMLWriter& xml, const std::string& s, const Ensem<T>& d)
{
  xml.openTag(s);
  xml << d;
  xml.closeTag();
}



//-----------------------------------------------------------------------------
//! Binary input
template<class T>
inline
void read(ADATIO::BinaryReader& bin, EScalar< OScalar<T> >& d)
{
#if 0
  int num, len, type;

  read(bin, num);
  read(bin, len);
  read(bin, type);

  if (num != 1)
  {
    fprintf(stderr, "BinaryReader is an invalid EScalar type\n");
    exit(1);
  }

  if (type != EnsbcIO< EScalar< OScalar<T> > >::type)
  {
    fprintf(stderr, "BinaryReader is an invalid OScalar type\n");
    exit(1);
  }
#endif

  read(bin, d.elem());
}

//! Binary input
template<class T>
inline
void read(ADATIO::BinaryReader& bin, EScalar< OVector<T> >& d)
{
#if 0
  int num, len, type;

  read(bin, num);
  read(bin, len);
  read(bin, type);

  if (num != 1)
  {
    fprintf(stderr, "BinaryReader is an invalid EScalar type\n");
    exit(1);
  }

  if (type != EnsbcIO< EScalar< OVector<T> > >::type)
  {
    fprintf(stderr, "BinaryReader is an invalid OScalar type\n");
    exit(1);
  }
#else
  int len;
  read(bin, len);
#endif

  d.resizeObs(len);
  read(bin, d.elem());
}

//! Binary input
template<class T>
inline
void read(ADATIO::BinaryReader& bin, EScalar< OTensor<T> >& d)
{
#if 0
  int num, len, type;

  read(bin, num);
  read(bin, len);
  read(bin, type);

  if (num != 1)
  {
    fprintf(stderr, "BinaryReader is an invalid ETensor type\n");
    exit(1);
  }

  if (type != EnsbcIO< EScalar< OTensor<T> > >::type)
  {
    fprintf(stderr, "BinaryReader is an invalid OScalar type\n");
    exit(1);
  }
#endif

  Array<int> nz;
  read(bin, nz);

  d.resizeObs(nz);
  read(bin, d.elem());
}


//! Binary input
template<class T>  
inline
void read(ADATIO::BinaryReader& bin, Ensem< OScalar<T> >& d)
{
  int num, len, type;

  read(bin, num);
  read(bin, len);
  read(bin, type);

  if (type != EnsbcIO< Ensem< OScalar<T> > >::type)
  {
    fprintf(stderr, "BinaryReader is an invalid EScalar type\n");
    exit(1);
  }

  d.checkResize(__func__, num, ENSEM_JACKKNIFE);
  for(int n = 0; n < num; ++n)
    read(bin, d.elem(n));
}

//! Binary input
template<class T>  
inline
void read(ADATIO::BinaryReader& bin, Ensem< OVector<T> >& d)
{
  int num, len, type;

  read(bin, num);
  read(bin, len);
  read(bin, type);

  if (type != EnsbcIO< Ensem< OVector<T> > >::type)
  {
    fprintf(stderr, "BinaryReader is an invalid OVector type\n");
    exit(1);
  }

  d.checkResize(__func__, num, ENSEM_JACKKNIFE);   // at the moment only jackknife
  d.resizeObs(len);
  for(int n = 0; n < num; ++n)
    read(bin, d.elem(n));
}


//! Binary output
template<class T>
inline
void write(ADATIO::BinaryWriter& bin, const EScalar< OScalar<T> >& d)
{
#if 0
  int num  = 1;
  int len  = 1;
  int type = EnsbcIO< EScalar< OScalar<T> > >::type;

  write(bin, num);
  write(bin, len);
  write(bin, type);
#endif

  write(bin, d.elem());
}

//! Binary output
template<class T>
inline
void write(ADATIO::BinaryWriter& bin, const EScalar< OVector<T> >& d)
{
#if 0
  int num  = 1;
  int len  = d.elem().size();
  int type = EnsbcIO< EScalar< OVector<T> > >::type;

  write(bin, num);
  write(bin, len);
  write(bin, type);
#else
  int len  = d.elem().size();
  write(bin, len);
#endif

  write(bin, d.elem());
}

//! Binary output
template<class T>
inline
void write(ADATIO::BinaryWriter& bin, const EScalar< OTensor<T> >& d)
{
#if 0
  int num  = 1;
  int len  = d.elem().numElem();
  int type = EnsbcIO< EScalar< OTensor<T> > >::type;

  write(bin, num);
  write(bin, len);
  write(bin, type);
#endif

  write(bin, d.elem().size());   // writes an array of dimensions
  write(bin, d.elem());
}


//! Binary output
template<class T>  
inline
void write(ADATIO::BinaryWriter& bin, const Ensem< OScalar<T> >& d)
{
  d.checkSize(__func__);

  int num  = d.size();
  int len  = d.numElem();
  int type = EnsbcIO< EScalar< OScalar<T> > >::type;

  write(bin, num);
  write(bin, len);
  write(bin, type);

  for(int n=0; n < num; ++n) 
    write(bin, d.elem(n));
}

//! Binary output
template<class T>  
inline
void write(ADATIO::BinaryWriter& bin, const Ensem< OVector<T> >& d)
{
  d.checkSize(__func__);

  int num  = d.size();
  int len  = d.elem(0).size();
  int type = EnsbcIO< EScalar< OVector<T> > >::type;

  write(bin, num);
  write(bin, len);
  write(bin, type);

  for(int n = 0; n < num; ++n)
    write(bin, d.elem(n));
}

//! Binary output
template<class T>  
inline
void write(ADATIO::BinaryWriter& bin, const Ensem< OTensor<T> >& d)
{
  d.checkSize(__func__);
  int num  = d.size();
  int len  = d.elem(0).numElem();
  int type = EnsbcIO< EScalar< OTensor<T> > >::type;

  write(bin, num);
  write(bin, len);
  write(bin, type);
  write(bin, d.elem(0).size()); // write an array of dimensions

  for(int n = 0; n < num; ++n)
    write(bin, d.elem(n));
}


//-----------------------------------------------------------------------------
// Serialization functions

//! Read from text
template<class T>
inline ENSEM::TextReader&
operator>>(ENSEM::TextReader& s, EScalar< OScalar<T> >& d)
{
  int num, len, type, ncol;
  int junk;

  char cline[256];
  s.getIstream().getline(cline, 256);
  std::istringstream line(cline);

  line >> num >> len >> type >> junk >> ncol;

  if (num != 1)
  {
    fprintf(stderr, "stream is an invalid EScalar type\n");
    exit(1);
  }

  if (type != EnsbcIO< EScalar< OScalar<T> > >::type)
  {
    fprintf(stderr, "stream is an invalid OScalar type\n");
    exit(1);
  }

  s >> d.elem();
  return s;
}


//! Read from text
template<class T>
inline ENSEM::TextReader&
operator>>(ENSEM::TextReader& s, EScalar< OVector<T> >& d)
{
  int num, len, type, ncol;
  int junk;

  char cline[256];
  s.getIstream().getline(cline, 256);
  std::istringstream line(cline);

  line >> num >> len >> type >> junk >> ncol;

  if (num != 1)
  {
    fprintf(stderr, "stream is an invalid EScalar type\n");
    exit(1);
  }

  if (type != EnsbcIO< EScalar< OVector<T> > >::type)
  {
    fprintf(stderr, "stream is an invalid OVector type\n");
    exit(1);
  }

  d.resizeObs(len);
  s >> d.elem();
  return s;
}


//! Read from text
template<class T>
inline ENSEM::TextReader&
operator>>(ENSEM::TextReader& s, EScalar< OTensor<T> >& d)
{
  int num, len, type, ncol;
  int junk;

  char cline[256];
  s.getIstream().getline(cline, 256);
  std::istringstream line(cline);

  line >> num >> len >> type >> junk >> ncol;

  if (num != 1)
  {
    fprintf(stderr, "stream is an invalid EScalar type\n");
    exit(1);
  }

  if (type != EnsbcIO< EScalar< OTensor<T> > >::type)
  {
    fprintf(stderr, "stream is an invalid OTensor type\n");
    exit(1);
  }

  int ndim;
  line >> ndim;

  Array<int> nz(ndim);
  for(int i=0; i < ndim; ++i)
    line >> nz[i];

  std::cout << "Tensor dimensions: ndim=" << ndim << ": dim=";
  for(int i=0; i < ndim; ++i)
    std::cout << " " << nz[i];
  std::cout << std::endl;

  d.resizeObs(nz);
  s >> d.elem();

  return s;
}


//! Read from text
template<class T>
inline ENSEM::TextReader&
operator>>(ENSEM::TextReader& s, Ensem< OScalar<T> >& d)
{
  int num, len, type, ncol;
  int junk;

  char cline[256];
  s.getIstream().getline(cline, 256);
  std::istringstream line(cline);

  line >> num >> len >> type >> junk >> ncol;

  if (type != EnsbcIO< Ensem< OScalar<T> > >::type)
  {
    fprintf(stderr, "stream is an invalid OScalar type\n");
    exit(1);
  }

  d.checkResize(__func__, num, ENSEM_JACKKNIFE);
  for(int n = 0; n < num; ++n)
    s >> d.elem(n);

  return s;
}


//! Read from text
template<class T>
inline ENSEM::TextReader&
operator>>(ENSEM::TextReader& s, Ensem< OVector<T> >& d)
{
  int num, len, type, ncol;
  int junk;

  char cline[256];
  s.getIstream().getline(cline, 256);
  std::istringstream line(cline);

  line >> num >> len >> type >> junk >> ncol;

  if (type != EnsbcIO< Ensem< OVector<T> > >::type)
  {
    fprintf(stderr, "stream is an invalid OVector type\n");
    exit(1);
  }

  d.checkResize(__func__, num, ENSEM_JACKKNIFE);   // at the moment only jackknife
  d.resizeObs(len);
  for(int n = 0; n < num; ++n)
    s >> d.elem(n);

  return s;
}


//! Read from text
template<class T>
inline ENSEM::TextReader&
operator>>(ENSEM::TextReader& s, Ensem< OTensor<T> >& d)
{
  int num, len, type, ncol;
  int junk;

  char cline[256];
  s.getIstream().getline(cline, 256);
  std::istringstream line(cline);

  line >> num >> len >> type >> junk >> ncol;

  if (type != EnsbcIO< Ensem< OTensor<T> > >::type)
  {
    fprintf(stderr, "stream is an invalid OTensor type\n");
    exit(1);
  }

  int ndim;
  line >> ndim;

  Array<int> nz(ndim);
  for(int i=0; i < ndim; ++i)
    line >> nz[i];

  std::cout << "Tensor dimensions: ndim=" << ndim << ": dim=";
  for(int i=0; i < ndim; ++i)
    std::cout << " " << nz[i];
  std::cout << std::endl;

  d.checkResize(__func__, num);
  d.resizeObs(nz);
  for(int n = 0; n < num; ++n)
    s >> d.elem(n);

  return s;
}



// Read from an ensbc file
template<class T>
inline
void read(const std::string& filename, EScalar<T>& d)
{
  ENSEM::TextReader f(filename);
  if (! f.is_open())
  {
    std::cerr << "read: error opening file" << filename << std::endl;
    exit(1);
  }

  f >> d;

  f.close();
}

// Read from an ensbc file
template<class T>
inline
void read(const std::string& filename, Ensem<T>& d)
{
  ENSEM::TextReader f(filename);
  if (! f.is_open())
  {
    std::cerr << "read: error opening file" << filename << std::endl;
    exit(1);
  }

  f >> d;

  f.close();
}


//! Print out to text
template<class T>
inline ENSEM::TextWriter&
operator<<(ENSEM::TextWriter& s, const EScalar< OScalar<T> >& d)
{
  s << " 1 1 " 
    << EnsbcIO< EScalar< OScalar<T> > >::type
    << " 0 1\n";

  s << d.elem();
  return s;
}


//! Print out to text
template<class T>
inline ENSEM::TextWriter&
operator<<(ENSEM::TextWriter& s, const EScalar< OVector<T> >& d)
{
  const int len = d.elem().size();

  s << " 1 " << len << " "
    << EnsbcIO< EScalar< OVector<T> > >::type 
    << " 0 1\n";

  s << d.elem();
  return s;
}


//! Print out to text
template<class T>
inline ENSEM::TextWriter&
operator<<(ENSEM::TextWriter& s, const EScalar< OTensor<T> >& d)
{
  const int len = d.elem().numElem();

  s << " 1 " << len << " " 
    << EnsbcIO< EScalar< OTensor<T> > >::type 
    << " 0 1 ";
  s << d.elem().size().size();   // number of dimensions
  for(int i = 0; i < d.elem().size().size(); ++i)
    s << " " << d.elem().size()[i];  // dimensions
  s << "\n";

  s << d.elem();
  return s;
}


//! Print out to text
template<class T>
inline ENSEM::TextWriter&
operator<<(ENSEM::TextWriter& s, const Ensem< OScalar<T> >& d)
{
  d.checkSize(__func__);
  const int num = d.size();

  s << num << " 1 "
    << EnsbcIO< Ensem< OScalar<T> > >::type 
    << " 0 1" << "\n";

  for(int n = 0; n < num; ++n)
  {
    s << d.elem(n);
    if (n < num-1)
      s << "\n";
  }
  return s;
}


//! Print out to text
template<class T>
inline ENSEM::TextWriter&
operator<<(ENSEM::TextWriter& s, const Ensem< OVector<T> >& d)
{
  d.checkSize(__func__);
  const int num = d.size();
  const int len = d.elem(0).size();

  s << num << " " << len << " " 
    << EnsbcIO< Ensem< OVector<T> > >::type 
    << " 0 1" << "\n";

  for(int n = 0; n < num; ++n)
  {
    s << d.elem(n);
    if (n < num-1)
      s << "\n";
  }
  return s;
}


//! Print out to text
template<class T>
inline ENSEM::TextWriter&
operator<<(ENSEM::TextWriter& s, const Ensem< OTensor<T> >& d)
{
  d.checkSize(__func__);
  const int num = d.size();
  const int len = d.elem(0).numElem();

  s << num << " " << len << " "
    << EnsbcIO< Ensem< OTensor<T> > >::type 
    << " 0 1 ";
  s << d.elem(0).size().size();   // number of dimensions
  for(int i = 0; i < d.elem(0).size().size(); ++i)
    s << " " << d.elem(0).size()[i];  // dimensions
  s << "\n";

  for(int n = 0; n < num; ++n)
  {
    s << d.elem(n);
    if (n < num-1)
      s << "\n";
  }
  return s;
}


// Write to an ensbc file
template<class T>
inline
void write(const std::string& filename, const EScalar<T>& d)
{
  ENSEM::TextWriter f(filename);

  if (! f.is_open())
  {
    std::cerr << "write: error opening file: " << filename << std::endl;
    exit(1);
  }

  f << d << "\n";

  f.close();
}

// Write to an ensbc file
template<class T>
inline
void write(const std::string& filename, const Ensem<T>& d)
{
  ENSEM::TextWriter f(filename);

  if (! f.is_open())
  {
    std::cerr << "write: error opening file: " << filename << std::endl;
    exit(1);
  }

  f << d << "\n";

  f.close();
}

} // namespace ENSEM
