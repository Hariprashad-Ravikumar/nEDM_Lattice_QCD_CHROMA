//$Id: exponper_altern.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: exponper_altern.h,v $
//Revision 2.0  2008/12/05 04:43:33  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.2  2007/03/08 05:54:24  kostas
//changed to be able to do arbitrary number of exponentials
//
//Revision 1.1  2006/12/27 04:23:31  kostas
//Added periodic alternating exponential fit function
//
//
/*! Definition of the sum of periodic exponentials plus alternating exponentials function */


#ifndef ALT_EXPONENTIALS_PERIODIC
#define ALT_EXPONENTIALS_PERIODIC

#include "covfit/Function.h"

namespace CovFit {

class AltExponentialsPeriodic: public Function{
private:
  double center ; // the t=0 point
  double period ; // periodic extent of the lattice

  double pExp(double& x, double& b){
    return exp(-b*x) + exp(-b*(period-x)) ;
  }

  double DpExp_db(double& x, double& b){
    return -x*exp(-b*x) -(period-x)*exp(-b*(period-x)) ;
  }

  double D2pExp_db2(double& x, double& b){
    return x*x*exp(-b*x) + (period-x)*(period-x)*exp(-b*(period-x)) ;
  }


public:
  
  

  AltExponentialsPeriodic():Function(){}
  AltExponentialsPeriodic(int n):center(0.0),period(0.0),Function(4*n){f_init();}
  AltExponentialsPeriodic(double c):center(c),period(0.0),Function(){f_init();}
  AltExponentialsPeriodic(double p, int n):center(0.0),period(p),Function(4*n){f_init();}
  AltExponentialsPeriodic(double c, double p):center(c),period(p),Function(){f_init();} 
  AltExponentialsPeriodic(double c, double p,int n):center(c),period(p),Function(2*n){f_init();}
  
  
  virtual ~AltExponentialsPeriodic(){}

  void setCenter(double c){
    center = c ;
  }

  void setPeriod(double p){
    period = p ;
  }

  void f_init(){
    if(npar%4!=0)
      std::cerr<<"AltExponentialsPeriodic::f_init : npar has to be multiple of four\n";
  }
  void f_init(FILE *fp)
  {
    int n ;
    fprintf(stderr,"sum of exponentials + alternating exponentials fitting\n");
    fprintf(stderr,"Enter the center:\n");
    fscanf(fp,"%lg",&center);
    fprintf(stderr,"Enter the period:\n");
    fscanf(fp,"%lg",&period);
    fprintf(stderr,"Enter the number of exponentials (each exponential is acompanied by an alternating exponential as well):\n");
    fscanf(fp,"%d",&n);
    npar = 4*n ;
    fprintf(stderr,"Using center: %g\n",center);
    fprintf(stderr,"Using period: %g\n",period);
    fprintf(stderr,"Number of parameters: %i\n",npar);
  }

#define NT period

  
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

    for(i=0;i<npar;i+=4) 
      z += b[i]*pExp(xx,b[i+1]) + s*b[i+2]*pExp(xx,b[i+3]);

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
      return s*b[i-1]*DpExp_db(xx,b[i]) ;
    else
      return s*pExp(xx,b[i+1]) ;


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
        return s*b[i-1]*D2pExp_db2(xx,b[i]) ;
      else if(j==i-1)
        return s*DpExp_db(xx,b[i]) ;
      else
        return 0.0 ;
    else
      if(j==i+1)
        return  s*DpExp_db(xx,b[i+1]) ;
      else
        return 0.0 ;
  }
  
#undef NT
  
};
  
} // namespace CovFit

#endif
