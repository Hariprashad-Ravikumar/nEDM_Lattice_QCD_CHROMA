//$Id: fitter.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: fitter.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.2  2008/11/30 02:49:34  edwards
//Put in some judicious const char* in place of char* to make the compiler
//happy.
//
//Revision 1.1  2008/11/12 21:37:18  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.2  2005/11/02 15:52:40  kostas
//removed a not needed include malloc.h
//
//Revision 1.1  2005/04/27 20:00:24  edwards
//Moved from include dir. Overhauled and now in own CovFit namespace.
//
//Revision 1.6  2004/05/20 17:33:38  kostas
//fixed params set through FitterArgs
//
//Revision 1.5  2004/05/20 16:32:23  kostas
//improved the  interface to the fitters and added xml parsing stuff
//for spectrum calculations
//
//Revision 1.4  2004/02/26 20:39:42  kostas
//added functionality to the fitters so that they can be called
//fro c++ programs
//
//Revision 1.3  2004/02/23 02:11:27  edwards
//Changed to include  Function.h instead of function.h
//
//Revision 1.2  2004/02/13 15:52:25  kostas
//Comments
//

//
// This is just a c++ wrap around the c fitting codes provided by
// Doug Toussaint. A description of the algorithms used can be
// found in Doug's TASI lectures.
//
#ifndef FITTER_H
#define FITTER_H

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>
#include "xml_array.h"
#include "covfit/Function.h"
#include "covfit/covarmat.h"


using namespace std;

namespace CovFit {

class FitterArgs{
public:
  double eps ;
  int maxit ;
  int display ;
  double xlow ; // 
  double xhigh; // fit x in xlow<= x<=xhigh
  int npar ;   // number of parameters
  Array<int> fixpar ; // number of fixed parameters
  FitterArgs():eps(1.0e-4),
	       maxit(10000),
	       display(1),
	       xlow(0),
	       xhigh(1),
	       npar(10){}
  ~FitterArgs(){}
} ;



class BaseFitter{
private:
protected:

  int ndata ;  // number of data points 
  double eps ;
  int maxit ;
  int display ;
  double xlow ; // 
  double xhigh; // fit x in xlow<= x<=xhigh
  int npar ;   // number of parameters
  int fixpar ; // number of fixed parameters
  double *xdata ;
  double *ydata ;
  double *sigma ;
  double  *par ; // fit parameters
  double *epar ; // error in fit parameters
  int *parlist ; // flag for fixed parameters
  int ndof ; // number of degrees of freedom
  double chi_square ; // chi-square
  double confidence ; // confidence if fit 
  bool singular_Hessian ;
  bool convergence ;
  bool fit_failed ;
  Function *F ;

  void subst(double **a, double *b,double *x, int *p, int *q, int n) ;
  int factor(double **a, int *p, int *q, int n) ;
  double norm(double *x, int n);
  double dot(double *x, double *y, int n);
  void lineq(double *mat, double *vec, double *ans, int dim) ;
  void matinv(double *x,double *y,int dim) ;

 public:
  BaseFitter(Function& f, const char* parf):npar(-1),ndata(0),singular_Hessian(false)
  { 
    F=&f ;
    xdata = (double *)malloc(sizeof(double));
    ydata = (double *)malloc(sizeof(double));
    sigma = (double *)malloc(sizeof(double));

    ReadParams(parf) ;
  }

  BaseFitter(Function& f, FitterArgs& Arg ):eps(Arg.eps),
					    maxit(Arg.maxit),
					    display(Arg.display),
					    xlow(Arg.xlow),
					    xhigh(Arg.xhigh),
					    npar(Arg.npar),
					    ndata(0),singular_Hessian(false)
  {
    F=&f ;
    npar = F->Npar() ;
    par = new double[npar] ;
    epar = new double[npar] ;
    parlist = new int[npar];
    xdata = (double *)malloc(sizeof(double));
    ydata = (double *)malloc(sizeof(double));
    sigma = (double *)malloc(sizeof(double));
    setFixedParams(Arg.fixpar) ;
  }


  virtual ~BaseFitter(){
    delete [] par ;
    delete [] epar ;
    delete [] parlist ;
    free(xdata);
    free(ydata);
    free(sigma);
  }

  bool  FailedFit(){return fit_failed||singular_Hessian ;}
  bool  Convergence(){return convergence;}
  Double ChiSq(){return chi_square;}
  Double Confidence(){return confidence ;}
  int Ndof(){return ndof;}

  void setRange(const Double& x,const Double& xx){ xlow = x; xhigh = xx ;}

  void ReadParams(const char* file2);
  virtual void dphi(double *b, double *grad)=0;
  virtual void ddphi(double *b, double **a)=0;
  virtual double phi(double *b)=0;
  //virtual void ReadData(char *fp)=0;
  virtual void Fit()=0;

  void getParams(Array<Double>& p, Array<Double>& e);

  void setParams(const Array<Double>& p );

  void setFixedParams(const Array<int>& f ) ;

  double gammq(double a, double x) ;
  void gser(double *gamser,double a,double x,double *gln) ;
  void gcf(double *gammcf,double a, double x,double *gln);
  double gammln(double xx) ;

  double min();
  void graph(char *fitgraph);
};

class Fitter : public BaseFitter {
private:
protected:
  double *weight ;
public:
  Fitter(Function& f, const char* pf, const char* df):BaseFitter(f,pf)
  {
    weight = (double *)malloc(sizeof(double));
    ReadData(df); 
    ndof = ndata-npar+fixpar;
  }

  Fitter(Function& f, const char* df):BaseFitter(f,"stdin")
  {
    weight = (double *)malloc(sizeof(double));
    ReadData(df); 
    ndof = ndata-npar+fixpar;
  }
  
  Fitter(Function& f,FitterArgs& A, const char* df):BaseFitter(f,A)
  {
    weight = (double *)malloc(sizeof(double));
    ReadData(df); 
    ndof = ndata-npar+fixpar; 
  }
  Fitter(Function& f,FitterArgs& A):BaseFitter(f,A)
  {
    weight = (double *)malloc(sizeof(double));
  }

  virtual ~Fitter(){
    free(weight);
  }

  void ReadData(const char* f) ;
  void setData(int N, double *x, double *y, double *w);
  void setData(Array<Double>& x, Array<Double>& y, Array<Double>& ey) ;

  void dphi(double *b, double *grad);
  void ddphi(double *b, double **a);
  double phi(double *b);
 
  void Fit();
};

class coFitter : public BaseFitter {
private:
  bool wd_alloc ;
  double *wdata1 ;
  double *wdata2 ;
  double *wdata3 ;
  double *wdata4 ;

  void dumpmat(const char *label,const char *format,double *matrix) ;

protected:
  int nblocks ; // number of data blocks after blocking
  double *covar ;
  double *covarinv ;
public:  
  coFitter(Function& f,char *pf):covar(NULL),wd_alloc(false),BaseFitter(f,pf){}
  coFitter(Function& f,char *pf, CovarMat& CM):wd_alloc(false),BaseFitter(f,pf)
  {
    SetCovarMat(CM); 
  }

  coFitter(Function& f, CovarMat& CM):wd_alloc(false),BaseFitter(f,"stdin")
  {
    SetCovarMat(CM); 
  }
  
  coFitter(Function& f,FitterArgs& A,CovarMat& CM):wd_alloc(false),
						   BaseFitter(f,A)
  {
    SetCovarMat(CM);
  }

  coFitter(Function& f,FitterArgs& A):covar(NULL),
				      wd_alloc(false),
				      BaseFitter(f,A){}

  virtual ~coFitter(){
    if(wd_alloc){
      delete [] wdata1 ;
      delete [] wdata2 ;
      delete [] wdata3 ;
      delete [] wdata4 ;
    }
    if(covar!=NULL){
      free(covar   );
      free(covarinv);
    }
  }

  double phi(double *b); 
  void dphi(double *b, double *grad);
  void ddphi(double *b, double **a);
 
  void Fit();

  void SetCovarMat(CovarMat& CM);

};

} // namespace CovFit

#endif
