//$Id: expon_altern.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: expon_altern.h,v $
//Revision 2.0  2008/12/05 04:43:33  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.3  2007/03/08 05:15:20  kostas
//fixed the problems alternating exps. had with spectrum
//
//Revision 1.2  2006/08/03 00:01:49  kostas
//fixed bugs
//
//Revision 1.1  2006/08/02 21:47:36  kostas
//added alternating exponential fits
//
/*! Definition of the sum of exponentials plus alternating exponentials function */


#ifndef ALT_EXPONENTIALS
#define ALT_EXPONENTIALS

#include "covfit/Function.h"

namespace CovFit {

class AltExponentials: public Function{
private:
  double center ; // the t=0 point
public:
  
  AltExponentials():Function(){}
  AltExponentials(int n):center(0.0),Function(4*n){f_init();}
  AltExponentials(double c):center(c),Function(){f_init();}
  AltExponentials(double c, int n):center(c),Function(4*n){f_init();}
  
  virtual ~AltExponentials(){}

  void setCenter(double c){
    center = c ;
  }

  void f_init(){
    if(npar%4!=0)
      std::cerr<<"AltExponentials::f_init : npar has to be multiple of four\n";
  }
  void f_init(FILE *fp)
  {
    int n ;
    fprintf(stderr,"sum of exponentials + alternating exponentials fitting\n");
    fprintf(stderr,"Enter the center:\n");
    fscanf(fp,"%lg",&center);
    fprintf(stderr,"Enter the number of exponentials (each exponential is acompanied by an alternating exponentials as well):\n");
    fscanf(fp,"%d",&n);
    npar = 4*n ;
    fprintf(stderr,"Using center: %g\n",center);
    fprintf(stderr,"Number of parameters: %i\n",npar);
  }
  
  double f(double x, double *b)
  {
    double z;
    int s(1), ix ;
    double xx(fabs(x-center));
    int i;
    z=0;
    ix = (int)(xx+0.5);
    s = 1 - 2*(ix%2) ;
    

    //printf("%i %i\n",ix,s);

    for(i=0;i<npar;i+=4) z += b[i]*exp(-b[i+1]*xx) + s*b[i+2]*exp(-b[i+3]*xx);
    return(z);
  }
  
  double df(double x, double *b, int i)
  {
    double xx(fabs(x-center));
    int s(1),ix ;
    
    if(i%4>1){
      ix = (int)(xx+0.5);
      s = 1 - 2*(ix%2) ;
    }
      

    if(i%2)
      return -s*xx*b[i-1]*exp(-b[i]*xx) ;
    else 
      return s*exp(-b[i+1]*xx) ;
  }
  
  double ddf(double x, double *b, int i, int j)
  {
    double xx(fabs(x-center));
    int s(1),ix ;
    
    if(i%4>1){
      ix = (int)(xx+0.5);
      s = 1 - 2*(ix%2) ;
    }
   

    if(i%2)
      if(j==i)
	return s*xx*xx*b[i-1]*exp(-b[i]*xx) ;
      else if(j==i-1)
	return -s*xx*exp(-b[i]*xx) ;
      else
	return 0.0 ;
    else
      if(j==i+1)
	return  -s*xx*exp(-b[i+1]*xx) ;
      else
	return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
