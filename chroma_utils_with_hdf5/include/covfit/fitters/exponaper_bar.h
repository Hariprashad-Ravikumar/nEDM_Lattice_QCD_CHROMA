//$Id: exponaper_bar.h,v 2.1 2009/01/23 04:09:07 kostas Exp $
//$Log: exponaper_bar.h,v $
//Revision 2.1  2009/01/23 04:09:07  kostas
//new fit function
//
//Revision 2.0  2008/12/05 04:43:33  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.2  2008/09/26 15:45:33  edwards
//Fixed some more missing includes to make the intel compiler happy.
//
//Revision 1.1  2005/11/23 03:18:51  kostas
//moved the fitter header files into the main include tree
//
//Revision 1.5  2005/04/27 19:53:52  edwards
//Completely overhauled include and lib structure. Put all
//of Kostas' fit stuff into their own covfit subdir in include
//and lib. All in its own namespace. Similarly, all of
//ensem is in it's own dirs and namespaces. Unified use of
//XML Array as base array type. I've not touched TTT which
//uses someother classes.
//The new convention should be is people put there little mini
//packages in their own dirs and namespaces. One can compose
//a bigger namespace with namespace composition.
//
//Revision 1.4  2005/02/02 22:01:30  kostas
//yet an other fix
//
//Revision 1.3  2005/02/02 20:28:30  kostas
//fixed bag
//
//Revision 1.2  2005/02/02 20:27:06  kostas
//fixed baryon fit function
//
//Revision 1.1  2004/09/26 19:23:58  kostas
//added antiperiodic exponential fit for baryons (not sure if it works)
//

/*! Definition of the sum of exponentials with anti-periodic BC function for
    baryon spectrum
*/


#ifndef EXPONENTIALS_APERIODIC_BARYON
#define EXPONENTIALS_APERIODIC_BARYON
//#define FUNCTION_NAME Exponentials
#include "covfit/Function.h"

namespace CovFit {

class ExponentialsAPeriodicBar: public Function{
private:
  double center ; // the t=0 point
  double period ; // periodic extent of the lattice
public:
  
  ExponentialsAPeriodicBar():Function(){}
  ExponentialsAPeriodicBar(int n):center(0.0),period(0.0),Function(4*n){f_init();}
  ExponentialsAPeriodicBar(double c):center(c),period(0.0),Function(){f_init();}
  ExponentialsAPeriodicBar(double p, int n):center(0.0),period(p),Function(4*n){f_init();}
  ExponentialsAPeriodicBar(double c, double p):center(c),period(p),Function(){f_init();}
  ExponentialsAPeriodicBar(double c, double p,int n):center(c),period(p),Function(4*n){f_init();}
  
  virtual ~ExponentialsAPeriodicBar(){}

  void setCenter(double c){
    center = c ;
  }

  void f_init(){
    if(npar%4!=0)
      std::cerr<<"ExponentialsAPeriodicBar::f_init : npar has to be a multiple of 4\n";
  }
  void f_init(FILE *fp)
  {
    int n ;
    fprintf(stderr,"sum of exponentials fitting with anti-periodic b.c. for baryons\n");
    fprintf(stderr,"Enter the center:\n");
    fscanf(fp,"%lg",&center);
    fprintf(stderr,"Enter the periodic extent of the lattice:\n");
    fscanf(fp,"%lg",&period);
    fprintf(stderr,"Enter the number of exponentials:\n");
    fscanf(fp,"%d",&n);
    npar = 4*n ;
    fprintf(stderr,"Using center: %g\n",center);
    fprintf(stderr,"Number of parameters: %i\n",npar);
  }
  
  double f(double x, double *b)
  {
    double z;
    double xx(fabs(x-center));
    double yy(period-xx);
    int i;
    z=0;
    for(i=0;i<npar;i+=4) z += b[i]*exp(-b[i+1]*xx)+b[i+2]*exp(-b[i+3]*yy);
    return(z);
  }
  
  double df(double x, double *b, int i)
  {
    double xx(fabs(x-center));
    if(i%4>1) xx = period-xx ;


    switch (i%2){
    case  0:
      return exp(-b[i+1]*xx) ; 
    case  1:
      return -xx*b[i-1]*exp(-b[i]*xx) ;
    }
    return 0.0;
  }
  
  double ddf(double x, double *b, int i, int j)
  {
    double xx(fabs(x-center));
    if(i%4>1) xx = period-xx ;

    switch (i%2){
    case  0:
      {
	if(j==i+1)
	  return -xx*exp(-b[i+1]*xx);
	else
	  return 0.0 ;
      }
      break;
    case  1:
	{
	  if(j==i-1)
	    return -xx*exp(-b[i]*xx) ;
	  else if(j==i)
	    return xx*xx*b[i-1]*exp(-b[i]*xx) ;
	  else
	    return 0.0 ;
	}
      break ;
    }
    return 0.0;
  }
  
};

} // namespace CovFit

#endif
