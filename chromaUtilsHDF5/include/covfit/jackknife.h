#include <iostream>
#include <xml_array.h>
#include <io/adat_io.h>
#include <io/adat_xmlio.h>
#include <covfit/statistics.h>
#include <covfit/proplist.h>
#include <covfit/Function.h>
#include <covfit/fitters/fitfunctions.h>
#include <covfit/fitter.h>


#ifndef JACKKNIFE_H
#define JACKKNIFE_H

namespace CovFit {
  using XMLArray::Array;
  using namespace ADATXML;
  using namespace std;

  class Value{
  public:
    Double v ;
    Double e ;

    Value(Double a, Double b):v(a),e(b){}
      /* This neat contructor makes possible to use the overloaded 
	 << operator to output convert Array<Double> to Value and then
	 write out the result!!  */
      Value(const Array<Double>& j):v(mean(j)),e(jackerr(j)){}
	
  } ;

  Value jackknife(const Array<Double>& d){
    return Value(mean(d),jackerr(d));
  }

  void read(XMLReader& xml, const string& path, Value& v)
  {
    XMLReader top(xml, path);
    
    read(top,"val" ,v.v);
    read(top,"err" ,v.e);
  }

  void write(XMLWriter& top, const string& path, const Value& c)
  {
    
    push(top,path);

    write(top,"val", c.v);
    write(top,"err", c.e);

    pop(top);
  }

  ostream& operator<<(ostream& s, const Value& p)
    {
      s << p.v << " +/- " << p.e;
    }



  template<typename T> 
    void jacklist(Array< Array<T> >& jc, const Array<T>& c, int block)
    {
      // care needs to be taken if the the division is not exact 
      int Ncnfs = c.size()/block ;
      jc.resize(c.size()) ;
      //Jackknife loop 
      for(int j(0);j<Ncnfs;j++)
	splice(jc[j],c,j,j+block-1) ; 
    
    }

  template<typename T> 
    void jackmean( Array<T>& c, const Array< Array<T> >& jc)
    {
      int Ncnfs = jc.size() ; 
      c.resize(Ncnfs);
      for(int j(0);j<Ncnfs;j++)
	c[j] = mean(jc[j]);
    }

  template<typename T> 
    void jackmean( Array<T>& mc, Array<T>& ec, const Array<T>& c, int block)
    {
      // care needs to be taken if the the division is not exact 
      int Ncnfs = c.size()/block ; 
      mc.resize(Ncnfs);
      ec.resize(Ncnfs);
      for(int j(0);j<Ncnfs;j++){
	Array<T> jc ;
	splice(jc,c,j,j+block-1) ; 
	mc[j] = mean(jc);
	ec[j] = jackerr(jc);
      }
    }

 template<typename T> 
   void jackmean( Array<T>& mc, const Array<T>& c, int block){
   Array<T> ec ;
   jackmean( mc, ec, c,  block) ;
 }

 template<typename T> 
   Array<T> effmass( const  Array<T>& y, const int& d ){
   Array<T> mass(y.size()) ;
   for(int t(0) ; t<y.size()-d;t++){
     mass[t] = log(y[t]/y[t+d])/double(d) ;
   }
   for(int k(1);k<d+1;k++)
     mass[y.size()-k] = mass[y.size()-d-1] ;

   return mass ;
 }

  template<typename T> 
   Array<T> effmass( const  Array<T>& y ){
    return effmass(y,1) ;
  }

 
 template<typename T> 
   Array< Array<T> >  effmass( const Array< Array<T> >& y, const int& d )
   {
     int Ncnfs = y.size() ;
     Array< Array<T> > mass(Ncnfs) ;
     for(int j(0);j<Ncnfs;j++)
       {
	 //cout<<j<<endl ;
	 mass[j] = effmass(y[j],d);
       }
     return mass ;
   }
 template<typename T> 
   Array< Array<T> >  effmass( const Array< Array<T> >& y ){
   return effmass(y,1) ;
 }

  template<typename T> 
    Array<T> effmass_cosh( const  Array<T>& y, const int& d ){
    Array<T> mass(y.size()) ;
    for(int t(d) ; t<y.size()-d;t++){
      mass[t] = acosh((y(t+d) + y(t-d))/(2.0*y(t)))/double(d) ;
    }
    // fix up boundary terms with something random...
    for(int k(1);k<d+1;k++){
      mass[y.size()-k] = mass[y.size()-d-1] ;
      mass[k-1] = mass[d] ;
    }
    
    return mass ;
  }

  template<typename T>
    Array<T> effmass_cosh( const  Array<T>& y ){
    return effmass_cosh(y,1);
  }

 template<typename T> 
   Array< Array<T> >  effmass_cosh( const Array< Array<T> >& y, const int& d )
   {
     int Ncnfs = y.size() ;
     Array< Array<T> > mass(Ncnfs) ;
     for(int j(0);j<Ncnfs;j++)
       {
	 //cout<<j<<endl ;
	 mass[j] = effmass_cosh(y[j],d);
       }
     return mass ;
   }

 template<typename T>
   Array< Array<T> >  effmass_cosh( const Array< Array<T> >& y ){
   effmass_cosh(y,1);
 }

 template<typename T> 
   Array<T> effmass_cosh_stag( const  Array<T>& y ){
   Array<T> mass(y.size()) ;
   for(int t(2) ; t<y.size()-2;t++){
     mass[t] = 0.5*acosh((y(t+2) + y(t-2))/(2.0*y(t))) ;
   }
   // fix up boundary terms with something random...
   mass[y.size()-1]  = mass[y.size()-2] = mass[y.size()-3] ;
   mass[0] = mass[1] = mass[2] ;
   return mass ;
 }

 template<typename T> 
   Array< Array<T> >  effmass_cosh_stag( const Array< Array<T> >& y )
   {
     int Ncnfs = y.size() ;
     Array< Array<T> > mass(Ncnfs) ;
     for(int j(0);j<Ncnfs;j++)
       {
	 //cout<<j<<endl ;
	 mass[j] = effmass_cosh_stag(y[j]);
       }
     return mass ;
   }


  template<typename T> 
    void jackfit( Array<Array<T> >& jfitpar, //jfitpar[0] carries the guess
		  Array<Double>& conf, //jfitpar[0] carries the guess
		  Array<Double>& chi2,
		  Fitter& fit,
		  Array<T>& x,
		  Array< Array<T> >& y,
		  Array< Array<T> >& e )
    {
      int Ncnfs = y.size() ;
      int npar = jfitpar[0].size();
      conf.resize(Ncnfs);
      chi2.resize(Ncnfs);
      Array<T> v(npar), ev(npar) ;
      v = jfitpar[0] ;
      for(int j(0);j<Ncnfs;j++){
	fit.setParams(v) ;
	fit.setData(x,y[j],e[j]) ;
	fit.Fit();
	fit.getParams(v, ev) ;
	jfitpar[j] = v ;
	conf[j] = fit.Confidence() ;
	chi2[j] = fit.ChiSq() ;
      }
    }




  template<typename T> 
    void jackfit( Array<Array<T> >& jfitpar, //jfitpar[0] carries the guess
		  Array<Double>& conf, //jfitpar[0] carries the guess
		  Array<Double>& chi2,
		  Fitter& fit,
		  Array<T>& x,
		  Array< Array<T> >& y,
		  const Array<T>& e )
    {
      Array< Array<T> > ee(y.size());
      for(int i(0);i<y.size();i++)
	ee[i] = e ;

      jackfit(jfitpar,conf,chi2,fit,x,y,ee);
    }

  void jackCovarMat(Array<CovarMat>& cm,
		    const Array<Double>& x,
		    const PropList& y,
		    int block)
  {
    
    int Ncnfs = y.size() ; 
    cm.resize(Ncnfs);
    for(int j(0);j<Ncnfs;j++){
      //cout<<"CODE LINE: "<<__LINE__<<"| jack sample: "<<j<<endl ;
      PropList jy ;
      splice(jy,y,j,j+block-1) ; 
      cm[j].Construct(x.size(),block) ;
      cm[j].SetDataList(x,jy);	
      cm[j].CalcCovarMat();
    }
  }
 
 template<typename T> 
   void jackcofit( Array<Array<T> >& jfitpar, //jfitpar[0] carries the guess
		   Array<Double>& conf, //jfitpar[0] carries the guess
		   Array<Double>& chi2,
		   coFitter& fit,
		   Array<CovarMat>& cm)
   {
      int Ncnfs = cm.size() ;
      int npar = jfitpar[0].size();
      conf.resize(Ncnfs);
      chi2.resize(Ncnfs);
      Array<T> v(npar), ev(npar) ;
      v = jfitpar[0] ;
      for(int j(0);j<Ncnfs;j++){
	fit.setParams(v) ;
	fit.SetCovarMat(cm[j]) ;
	fit.Fit();
	fit.getParams(v, ev) ;
	jfitpar[j] = v ;
	conf[j] = fit.Confidence() ;
	chi2[j] = fit.ChiSq() ;
      }
   }


  
}

#endif 
