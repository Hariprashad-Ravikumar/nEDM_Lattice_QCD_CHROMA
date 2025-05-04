//$Id: gAJfunc.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: gAJfunc.h,v $
//Revision 2.0  2008/12/05 04:43:33  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2006/07/12 00:20:26  kostas
//added some things about gA fits.
//
//Revision 1.1  2006/04/24 20:56:08  kostas
//Added the strange functions that apear in HBchiPT
//
//
/*!  ChiPT basic function for 1-loop Baryon fits
*/


#ifndef CHIPT_gAJfunc
#define CHIPT_gAJfunc
#include "covfit/Function.h"

namespace CovFit {

class gAJfunc: public Function{
private:
    double Delta ;
    double mu2 ;
    double Delta2 ;
public:
  
  gAJfunc(double D, double sc):Function(0),Delta(D/sc){
    mu2= (sc*sc) ; Delta2 = (Delta*Delta) ;
    npar=0 ; 
  }
  
  virtual ~gAJfunc(){}


  void f_init(){
    
  }
  void f_init(FILE *fp)
  {
    fprintf(stderr,"One loop chiPT for J fuction (apearing in HBchiPT for gA\n");
    fprintf(stderr,"J(m,Delta,mu)  =") ;
    fprintf(stderr,"(m^2-2D^2)log(m^2/mu^2) + 2Dsqrt(D^2-m^2)*log((D - sqrt(D^2-m^2)/(D + sqrt(D^2-m^2))))\n     x is m^2/mu^2\n  Delta, m and mu should be in the same units");
  }
  
  // x is mpi^2/mu^2 Delta is  in units of mu 
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
	//std::cout<<"\nroot   : "<<root<<std::endl;
	//std::cout<<"Log_arg: "<<((Delta-root)/(Delta+root))<<std::endl;
	//std::cout<<"Log    : "<<log((Delta-root)/(Delta+root))<<std::endl;
	//std::cout<<"rootLog: "<<sqrtChiLog<<std::endl;
	
      }
            
      return (x-2.0*Delta2)*chiLog + 2.0*Delta*sqrtChiLog ;

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
