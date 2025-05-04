// -*- C++ -*-
//$Id: polynomium.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: polynomium.h,v $
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
//Revision 1.4  2004/03/17 22:30:51  kostas
//modified how to compile fitfunctions
//
//Revision 1.3  2004/02/26 20:37:56  kostas
//fixed problems with automake
//
//Revision 1.2  2004/02/23 02:35:47  edwards
//Fix ups. Should now compile
//
//Revision 1.1  2004/02/13 15:32:18  kostas
//Added Doug's fitting programs
//

/*! Definition of the polynomium function */


#ifndef POLYNOMIUM
#define POLYNOMIUM
//#define FUNCTION_NAME Polynomium
#include "covfit/Function.h"

namespace CovFit {

class Polynomium: public Function{

public:
  
  Polynomium():Function(){}
  Polynomium(int n):Function(n){}
  virtual ~Polynomium(){}
  
  void f_init(FILE *fp)
  {
    /* Nothing to be done */
  }
  
  double f(double x, double *b)
  {
    double z;
    int i;
    z=0;
    for(i=npar-1;i>=0;i--) z = b[i]+x*z;
    return(z);
  }
  
  double df(double x, double *b, int i)
  {
    if(i==0)return(1.0);
    else return( pow(x,(double)i) );
  }
  
  double ddf(double x, double *b, int i, int j)
  {
    return(0.0);
  }
  
};

} // namespace CovFit

#endif
