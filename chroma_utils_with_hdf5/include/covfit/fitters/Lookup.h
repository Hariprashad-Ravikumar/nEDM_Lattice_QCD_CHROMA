//$Id: Lookup.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: Lookup.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2006/03/14 18:39:38  kostas
//Lookup table function
//
//
/*! LookUp table function
*/


#ifndef CHIPT_LookUp
#define CHIPT_LookUp
#include "covfit/Function.h"

namespace CovFit {

class LookUp: public Function{
private:
  double Tol ;

  Array<Double> xx ;
  Array<Double> yy ;

public:
  
  LookUp(const Array<Double>& x_, const Array<Double>& y_):
    Function(0),xx(x_),yy(y_){
      //npar=3 ;
      Tol = 1.0e-6 ;
  }

  LookUp(const Array<Double>& x_, const Array<Double>& y_,const Double& f):
    Function(0),xx(x_),yy(y_){
      //npar=3 ;
      Tol =f ;
  }
  
    int Nval(){xx.size();} 

  virtual ~LookUp(){}


  void f_init(){
    
  }
  void f_init(FILE *fp)
  {
    fprintf(stderr,"Initialization of lookup table\n");
    fprintf(stderr," Using %i parameters and tolerance %e\n",npar/2,Tol) ;
  }
  
  double f(double x, double *b)
    {
      for( int i(0);i<xx.size();i++)
	if( fabs(xx[i] - x)<Tol )
	  return yy[i] ;

      return -100000.00 ;
    }
  
  double df(double x, double *b, int i)
    {
      return 1.0 ;
    }
  
  double ddf(double x, double *b, int i, int j)
  {
    return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
