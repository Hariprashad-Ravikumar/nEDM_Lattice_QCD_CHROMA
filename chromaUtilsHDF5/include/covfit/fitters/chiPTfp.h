//$Id: chiPTfp.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: chiPTfp.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2005/11/23 03:18:51  kostas
//moved the fitter header files into the main include tree
//
//Revision 1.3  2005/10/05 14:11:26  kostas
//added mpi^4 term in chiPTfpi and also added it to the makefiles
//
//Revision 1.2  2005/10/03 16:04:02  kostas
//fixed bugs
//
//Revision 1.1  2005/10/03 15:52:52  kostas
//fpi chipt fit function
//
/*! ChiPT formula for Fpi
*/


#ifndef CHIPT_FPI
#define CHIPT_FPI
#include "covfit/Function.h"

namespace CovFit {

class ChiPTfpi: public Function{
private:
  double inv_EIGHT_PI_SQ ;

public:
  
  ChiPTfpi():Function(3){
    npar=3 ; 
    inv_EIGHT_PI_SQ = 1.0/78.95683520871486895155 ;}
  
  virtual ~ChiPTfpi(){}


  void f_init(){
    
  }
  void f_init(FILE *fp)
  {
    fprintf(stderr,"One loop chiPT for Fpi\n");
    fprintf(stderr,"Fpi = f ( 1 - 1/(8*pi^2 f^2) x log x ) + e2 x + e4 x^2\n");
    fprintf(stderr,"b[0] = f   b[1] = e2   b[2] = e4\n");
  }
  
  // x is mpi^2 b[0] = f b[1]= quadratic term 
  double f(double x, double *b)
    {
      return ( b[0]*(1.0 - inv_EIGHT_PI_SQ*x/(b[0]*b[0])*log(x)) + 
	       b[1]*x + b[2]*x*x ) ;
    }
  
  double df(double x, double *b, int i)
    {
      if(i==0)
	return 1.0 + inv_EIGHT_PI_SQ*x/(b[0]*b[0])*log(x) ;
      else if( i==1 )
	return x ;
      else
	return x*x ;
    }
  
  double ddf(double x, double *b, int i, int j)
  {
    if((i==0) && (j == 0))
      return - 2.0 * inv_EIGHT_PI_SQ*x/(b[0]*b[0]*b[0])*log(x) ;
    else
      return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
