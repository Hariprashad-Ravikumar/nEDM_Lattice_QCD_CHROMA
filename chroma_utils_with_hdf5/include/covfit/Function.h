//$Id: Function.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: Function.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:18  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2005/04/27 20:00:24  edwards
//Moved from include dir. Overhauled and now in own CovFit namespace.
//
//Revision 1.2  2004/02/13 15:52:25  kostas
//Comments
//
#ifndef FUNCTION_H
#define FUNCTION_H

#include <cstdio>
#include <cmath>

namespace CovFit {

//! The function base class
class Function{
protected:
  int npar ;
public:
  Function(int N):npar(N){}
  Function():npar(-1){}
  virtual ~Function(){}
  virtual void f_init(FILE *fp)=0;
  virtual double f(double x, double *b)=0;
  virtual double df(double x, double *b, int i)=0;
  virtual double ddf(double x, double *b, int i, int j)=0;
  int Npar(){return npar;}
  void setNpar(int n){npar=n;}
};

} // namespace CovFit

#endif
