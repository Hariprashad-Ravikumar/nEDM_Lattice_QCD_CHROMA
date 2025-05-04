//$Id: exponper_altern_onemass.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: exponper_altern_onemass.h,v $
//Revision 2.0  2008/12/05 04:43:33  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2007/03/08 05:57:09  kostas
//added new fit function
//
//
//
/*! Definition of the sum of periodic exponentials plus alternating exponentials function with a common mass. Used for mixed staggered/Wilson-like mesons
*/


#ifndef ALT_EXPONENTIALS_PERIODIC_ONEMASS
#define ALT_EXPONENTIALS_PERIODIC_ONEMASS

#include "covfit/Function.h"

namespace CovFit {

class AltExponentialsPeriodicOneMass: public Function{
private:
  double center ; // the t=0 point
  double period ; // periodic extent of the lattice
public:
  
  

  AltExponentialsPeriodicOneMass():Function(){}
  AltExponentialsPeriodicOneMass(int n):center(0.0),period(0.0),Function(3*n){f_init();}
  AltExponentialsPeriodicOneMass(double c):center(c),period(0.0),Function(){f_init();}
  AltExponentialsPeriodicOneMass(double p, int n):center(0.0),period(p),Function(3*n){f_init();}
  AltExponentialsPeriodicOneMass(double c, double p):center(c),period(p),Function(){f_init();} 
  AltExponentialsPeriodicOneMass(double c, double p,int n):center(c),period(p),Function(3*n){f_init();}
  
  
  virtual ~AltExponentialsPeriodicOneMass(){}

  void setCenter(double c){
    center = c ;
  }

  void setPeriod(double p){
    period = p ;
  }

  void f_init(){
    if(npar%3!=0)
      std::cerr<<"AltExponentialsPeriodicOneMass::f_init : npar has to be multiple of three\n";
  }
  void f_init(FILE *fp)
  {
    int n ;
    fprintf(stderr,"sum of exponentials + alternating exponentials fitting\n");
    fprintf(stderr,"one mass per pair (A + (-1)^t B)(exp(-m t) + image\n");
    fprintf(stderr,"Enter the center:\n");
    fscanf(fp,"%lg",&center);
    fprintf(stderr,"Enter the period:\n");
    fscanf(fp,"%lg",&period);
    fprintf(stderr,"Enter the number of exponentials (each exponential is acompanied by an alternating exponential as well):\n");
    fscanf(fp,"%d",&n);
    npar = 3*n ;
    fprintf(stderr,"Using center: %g\n",center);
    fprintf(stderr,"Using period: %g\n",period);
    fprintf(stderr,"Number of parameters: %i\n",npar);
  }

#define NT period

  double f(double x, double *b)
  /*  b[0] = coefficient 1
   *  b[1] = coefficient 2 
   *  b[2] = mass 
   *   ........
   */
  {
    double z;
    int s(1), ix ;
    double xx(fabs(x-center));
    int i;
    z=0;
    ix = (int)(xx+0.5);
    s = 1 - 2*(ix%2) ;

    for(int i(0);i<npar;i+=3) 
      z += (b[i]+ s*b[i+1])*(exp(-b[i+2]*xx) + exp(-b[i+2]*(NT-xx)));

    return(z);

  }
  
  double df(double x, double *b, int i)
  {
    double xx(fabs(x-center));
    int s(1),ix ;
    ix = (int)(xx+0.5);
    s = 1 - 2*(ix%2) ;

    if(i%3 == 0)
      return  (exp(-b[i+2]*xx) + exp(-b[i+2]*(NT-xx))) ;
    else if (i%3 == 1)  
      return s*(exp(-b[i+1]*xx) + exp(-b[i+1]*(NT-xx))) ;
    else if(i%3 == 2)
      return (b[i-2]+ s*b[i-1])*
	(-xx*exp(-b[i]*xx) - (NT -xx)*exp(-b[i]*(NT-xx)));
    else
      return 0.0 ;
  }
  
  double ddf(double x, double *b, int i, int j)
  {
    double xx(fabs(x-center));
    int s(1),ix ;
    
    ix = (int)(xx+0.5);
    s = 1 - 2*(ix%2) ;
    
    if(i%3==0)
      if( j== i + 2 )
	return (-xx*exp(-b[j]*xx) -(NT-xx)*exp(-b[j]*(NT-xx))) ;
      else
	return 0.0 ;
    else if ( i%3 == 1)
      if( j == 1+i )
        return s*(-xx*exp(-b[j]*xx) -(NT-xx)*exp(-b[j]*(NT-xx))) ;
      else
        return 0.0 ;
    else if (i%3 == 2)
      if( j == i )
	return (b[i-2]+ s*b[i-1])*(xx*xx*exp(-b[i]*xx) + 
				   (NT-xx)*(NT-xx)*exp(-b[i]*(NT-xx)));
      else if(j == i - 1 )
	return s*(-xx*exp(-b[i]*xx) - (NT -xx)*exp(-b[i]*(NT-xx))) ;
      else if (j == i - 2 )
	return (-xx*exp(-b[i]*xx) - (NT -xx)*exp(-b[i]*(NT-xx))) ;
      else 
	return 0.0 ;
    else 
      return(0.0);
  }
  
#undef NT
  
};
  
} // namespace CovFit

#endif
