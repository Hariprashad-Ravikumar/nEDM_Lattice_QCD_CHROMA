// -*- C++ -*-
//
// ENSEM data parallel interface
//

namespace ENSEM {


//
// Conversion routines. These cannot be implicit conversion functions
// since they foul up the PETE defs in ENSEMOperators.h using primitive
// types
//

//! Make an int from an Integer
inline int 
toInt(const Integer& s) 
{
  return toInt(s.elem());
}

//! Make a float from a Real32
inline float
toFloat(const Real32& s) 
{
  return toFloat(s.elem());
}

//! Make a double from a Real64
inline double
toDouble(const Real64& s) 
{
  return toDouble(s.elem());
}

//! Make a bool from a Boolean
inline bool
toBool(const Boolean& s) 
{
  return toBool(s.elem());
}


//
// Return an equivalent ENSEM type given some simple machine type
//
template<>
struct SimpleScalar<float>
{
  typedef Real32   Type_t;
};

// Construct simple float word
template<>
struct SimpleScalar<int>
{
  typedef Integer   Type_t;
};

// Construct simple double word
template<>
struct SimpleScalar<double>
{
  typedef Real64   Type_t;
};

// Construct simple boolean word
template<>
struct SimpleScalar<bool>
{
  typedef Boolean   Type_t;
};


//
// Type constructors for ENSEM types within the type system. Namely,
// at some level like a primitive, sometimes scalar temporaries are needed
// These are the bottom most constructors given a machine type
//
// Construct simple word
template<>
struct InternalScalar<float>
{
  typedef float  Type_t;
};

template<>
struct InternalScalar<int>
{
  typedef int   Type_t;
};

template<>
struct InternalScalar<double>
{
  typedef double   Type_t;
};

template<>
struct InternalScalar<bool>
{
  typedef bool  Type_t;
};


// Makes a primitive scalar leaving grid alone
template<>
struct PrimitiveScalar<float>
{
  typedef float  Type_t;
};

template<>
struct PrimitiveScalar<int>
{
  typedef int   Type_t;
};

template<>
struct PrimitiveScalar<double>
{
  typedef double   Type_t;
};

template<>
struct PrimitiveScalar<bool>
{
  typedef bool  Type_t;
};



// Makes a ensem scalar leaving primitive indices alone
template<>
struct EnsemScalar<float>
{
  typedef float  Type_t;
};

template<>
struct EnsemScalar<int>
{
  typedef int   Type_t;
};

template<>
struct EnsemScalar<double>
{
  typedef double   Type_t;
};

template<>
struct EnsemScalar<bool>
{
  typedef bool  Type_t;
};



// Internally used real scalars
template<>
struct RealScalar<int> {
  typedef REAL32  Type_t;
};

template<>
struct RealScalar<float> {
  typedef REAL32  Type_t;
};

template<>
struct RealScalar<double> {
  typedef REAL64  Type_t;
};

} // namespace ENSEM


namespace ENSEM {

// XML readers
void read(ADATXML::XMLReader& xml, const std::string& s, Array<Integer>& d);
void read(ADATXML::XMLReader& xml, const std::string& s, Array<Real32>& d);
void read(ADATXML::XMLReader& xml, const std::string& s, Array<Real64>& d);
void read(ADATXML::XMLReader& xml, const std::string& s, Array<Boolean>& d);

// XML writers
void write(ADATXML::XMLWriter& xml, const std::string& s, const Array<Integer>& d);
void write(ADATXML::XMLWriter& xml, const std::string& s, const Array<Real32>& d);
void write(ADATXML::XMLWriter& xml, const std::string& s, const Array<Real64>& d);
void write(ADATXML::XMLWriter& xml, const std::string& s, const Array<Boolean>& d);

} // namespace Ensem

