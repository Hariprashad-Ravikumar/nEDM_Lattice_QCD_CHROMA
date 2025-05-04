//$Id: expon.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: expon.h,v $
//Revision 2.0  2008/12/05 04:43:33  edwards
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
//Revision 1.4  2005/04/27 19:53:52  edwards
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
//Revision 1.3  2004/03/18 14:55:14  kostas
//fixed bug
//
//Revision 1.2  2004/03/17 23:34:44  kostas
//small fix bugs
//
//Revision 1.1  2004/03/17 22:30:51  kostas
//modified how to compile fitfunctions
//
/*! Definition of the sum of exponentials function */


#ifndef EXPONENTIALS
#define EXPONENTIALS
//#define FUNCTION_NAME Exponentials
#include "covfit/Function.h"

namespace CovFit {

class Exponentials: public Function{
private:
  double center ; // the t=0 point
public:
  
  Exponentials():Function(){}
  Exponentials(int n):center(0.0),Function(2*n){f_init();}
  Exponentials(double c):center(c),Function(){f_init();}
  Exponentials(double c, int n):center(c),Function(2*n){f_init();}
  
  virtual ~Exponentials(){}

  void setCenter(double c){
    center = c ;
  }

  void f_init(){
    if(npar%2!=0)
      std::cerr<<"Exponentials::f_init : npar has to be even\n";
  }
  void f_init(FILE *fp)
  {
    int n ;
    fprintf(stderr,"sum of exponentials fitting\n");
    fprintf(stderr,"Enter the center:\n");
    fscanf(fp,"%lg",&center);
    fprintf(stderr,"Enter the number of exponentials:\n");
    fscanf(fp,"%d",&n);
    npar = 2*n ;
    fprintf(stderr,"Using center: %g\n",center);
    fprintf(stderr,"Number of parameters: %i\n",npar);
  }
  
  double f(double x, double *b)
  {
    double z;
    double xx(fabs(x-center));
    int i;
    z=0;
    for(i=0;i<npar;i+=2) z += b[i]*exp(-b[i+1]*xx);
    return(z);
  }
  
  double df(double x, double *b, int i)
  {
    double xx(fabs(x-center));
    if(i%2)
      return -xx*b[i-1]*exp(-b[i]*xx) ;
    else 
      return exp(-b[i+1]*xx) ;
  }
  
  double ddf(double x, double *b, int i, int j)
  {
    double xx(fabs(x-center));
    if(i%2)
      if(j==i)
	return xx*xx*b[i-1]*exp(-b[i]*xx) ;
      else if(j==i-1)
	return -xx*exp(-b[i]*xx) ;
      else
	return 0.0 ;
    else
      if(j==i+1)
	return  -xx*exp(-b[i+1]*xx) ;
      else
	return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
