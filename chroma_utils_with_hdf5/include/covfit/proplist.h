#include <iostream>
#include <cstdio>
#include <complex>

#include "xml_array.h"
#include "xml_array2d.h"

#ifndef PROPLIST_H
#define PROPLIST_H

namespace CovFit {

  // namespace composition
  using XMLArray::Array;
  using XMLArray::Array2d;

  // Types snarfed from ADAT. These should be merged somehow
  typedef float  Real;
  typedef double Double;
  typedef int    Integer;
  typedef bool   Boolean;

  typedef std::complex<Double> DComplex ;
  typedef std::complex<Real> Complex ;

  /**  
  struct Complex
  {
    Real  re;
    Real  im;
  };

  struct DComplex
  {
    Double  re;
    Double  im;
  };
  **/
  
  // CovFit specific
  typedef Array< Array<Double> > PropList ;
  typedef Array< Array<DComplex> > cPropList ;

  //! A splice routine usefull for Jackknife
  template <typename T> 
  void splice(Array<T>& out, const Array<T>& in, int start, int end){
    if((in.size()-(end-start+1))>0){
      //cout<<"splice size: "<<in.size()-(end-start+1)<<endl;
      out.resize(in.size()-(end-start+1));
      //cout<<"splice: "<<out.size()<<endl;
      for(int i(0);i<start;i++)
	out[i]=in[i] ;
      int i;
      for(int k(end+1),i=start;k<in.size();k++,i++)
	out[i] = in[k] ;
    }
    else{
      std::cerr<<"splice: Size missmatch\n" ;
      exit(312);
    }
  }


  //! A block routine                                                     
  template <typename T>
  void arrayblock(Array<T>& out, const Array<T>& in, int start, int end){
    if((in.size()>end)&&(start>-1)){
      out.resize(end-start+1);
      //cout<<"arrayblock: "<<out.size()<<endl;                        
      for(int i(0);i<out.size();i++)
        out[i]=in[start+i] ;
    }
    else{
      std::cerr<<"arrayblock: out of bounds\n" ;
      exit(3112);
    }
  }


  int ReadProplist(PropList& pr, const std::string& file);
  int ReadProplist(cPropList& pr, const std::string& file);

  void WriteProplist(PropList& pr, int Nx, const std::string& file);
  void WriteProplist(cPropList& pr, int Nx, const std::string& file);

  int ReadProplist(Array< Array< Array<Double> > >& pr, const std::string& file);
  int ReadProplist(Array< Array< Array<DComplex> > >& pr, const std::string& file);

  void WriteProplist(Array< Array< Array<Double> > >& pr, int Nx, const std::string& file);
  void WriteProplist(Array< Array< Array<DComplex> > >& pr, int Nx, const std::string& file);
  
  int ReadProplist(Array< Array< Array2d<Double> > >& pr, const std::string& file); 
 
  void WriteProplist(Array< Array< Array2d<Double> > >& pr, int Nx, const std::string& file);

  int ReadProplist(Array< Array< Array2d<DComplex> > >& pr, const std::string& file) ; 

  void WriteProplist(Array< Array< Array2d<DComplex> > >& pr, int Nx, const std::string& file); 
 
}  // namespace CovFit

#endif
