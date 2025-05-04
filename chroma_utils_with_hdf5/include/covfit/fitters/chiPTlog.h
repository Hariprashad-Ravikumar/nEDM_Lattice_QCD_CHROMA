//$Id: chiPTlog.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: chiPTlog.h,v $
//Revision 2.0  2008/12/05 04:43:33  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2006/05/17 21:14:18  kostas
//added genneric chiPT log fit function
//
//
/*! genneric ChiPT logarithm
*/


#ifndef CHIPT_LOG
#define CHIPT_LOG
#include "covfit/Function.h"

namespace CovFit {

class ChiPTlog: public Function{
private:

public:
  
  ChiPTlog():Function(3){
    npar=4 ; 
    }
    

  virtual ~ChiPTlog(){}


  void f_init(){
    
  }
  void f_init(FILE *fp)
  {
    fprintf(stderr,"Generic One loop chiPT \n");
    fprintf(stderr,"Fpi = e0 + e1 x log x  + e2 x + e4 x^2\n");
    fprintf(stderr,"b[0] = e0   b[1] = e1   b[2] = e2 b[3] =  e4\n");
  }
  
  // x is mpi^2 b[0] = f b[1]= quadratic term 
  double f(double x, double *b)
    {
      return ( b[0] + b[1]*x*log(x) + b[2]*x + b[3]*x*x ) ;
    }
  
  double df(double x, double *b, int i)
    {
      if(i==0)
	return 1.0 ;
      else if( i==1 )
	return x*log(x) ;
      else if( i==2 )
	return x ;
      else
	return x*x ;
    }
  
  double ddf(double x, double *b, int i, int j)
  {
    return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
