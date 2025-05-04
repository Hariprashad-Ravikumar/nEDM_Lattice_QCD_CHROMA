//$Id: BaryonSfunc.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: BaryonSfunc.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2006/04/24 20:56:08  kostas
//Added the strange functions that apear in HBchiPT
//
//
/*! ChiPT basic function for 1-loop Baryon fits
*/


#ifndef CHIPT_BaryonSfunc
#define CHIPT_BaryonSfunc
#include "covfit/Function.h"

namespace CovFit {

class BaryonSfunc: public Function{
private:
    double Delta ;
    double mu2 ;
    double Delta2 ;
public:
  
  BaryonSfunc(double D, double sc):Function(0),Delta(D){
    mu2= (sc*sc) ; Delta2 = (D*D) ;
    npar=0 ; 
  }
  
  virtual ~BaryonSfunc(){}


  void f_init(){
    
  }
  void f_init(FILE *fp)
  {
    fprintf(stderr,"One loop chiPT for F fuction (apearing in HBchiPT\n");
    fprintf(stderr,"S(m,Delta,mu)  =") ;
    fprintf(stderr,"sqrt(D^2-m^2)*log((D - sqrt(D^2-m^2)/(D + sqrt(D^2-m^2))) - D*(log(m^2/mu^2)+1/3)\n");
  }
  
  // x is mpi^2 
  double f(double x, double *b)
    {

      double chiLog(log(x/mu2)); 
      double sqrtChiLog ;
      double root ;
      if(Delta2<x){
	root =  sqrt(x -  Delta2) ;
	sqrtChiLog = 2.0*root*atan(root/Delta) ;
	//return (x-Delta2)*(root*2.0*theta-Delta*chiLog)-0.5*x*Delta *chiLog ;
      }
      else{
	root =  sqrt(Delta2-x) ;
	sqrtChiLog = root*log((Delta-root)/(Delta+root)) ;
      }
      
      return sqrtChiLog - Delta*(chiLog+0.33333) ;

    }
  
  double df(double x, double *b, int i)
    {
      return 0.0 ;
    }
  
  double ddf(double x, double *b, int i, int j){
    return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
