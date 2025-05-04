//$Id: covarmat.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: covarmat.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:18  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.5  2008/09/26 15:45:33  edwards
//Fixed some more missing includes to make the intel compiler happy.
//
//Revision 1.4  2005/11/30 16:30:34  kostas
//fixed memory leak bug in CovarMat (when Construct was called)
//
//Revision 1.3  2005/06/03 01:49:12  kostas
//jackknifed covariance matrix
//
//Revision 1.2  2005/06/02 21:35:34  kostas
//trivial modification
//
//Revision 1.1  2005/04/27 20:00:24  edwards
//Moved from include dir. Overhauled and now in own CovFit namespace.
//
//Revision 1.5  2004/05/19 22:30:14  kostas
//fixed yet another bug
//
//Revision 1.4  2004/05/19 21:53:12  kostas
//fixed bug
//
//Revision 1.3  2004/05/19 21:41:44  kostas
//added some Arrays in CovarMat
//
//Revision 1.2  2004/02/13 15:52:25  kostas
//Comments
//

//
// This is just a c++ wrap around the c fitting codes provided by
// Doug Toussaint. A description of the algorithms used can be
// found in Doug's TASI lectures.
//
#ifndef COVARMAT_H
#define COVARMAT_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <string.h>
#include "xml_array.h"
#include "covfit/proplist.h"

namespace CovFit {

class CovarMat{
  double *ylist ;
  double *covar ;
  double *xdata ;
  double *ydata ;
  int npoints   ;
  int blocksize ;
  int nblocks   ;
  int nprops    ;
 public:
  CovarMat():npoints(0),blocksize(-1),nblocks(0),nprops(-1){
    xdata = (double *)malloc(sizeof(double));
    ydata = (double *)malloc(sizeof(double));
    ylist = (double *)malloc(sizeof(double));
    covar = (double *)malloc(sizeof(double));
  }
  CovarMat(int np, int bl):npoints(np),blocksize(bl)
    {
      xdata = (double *)malloc(npoints*sizeof(double));
      ydata = (double *)malloc(npoints*sizeof(double));
      ylist = (double *)malloc(        sizeof(double));
      covar = (double *)malloc(npoints*npoints*sizeof(double));
      nprops  = 0 ;
      nblocks = 0 ;
    }

  void Construct(int np, int bl){
    npoints = np ;
    blocksize = bl ;
    
    xdata = (double *)realloc(xdata,         npoints*sizeof(double));
    ydata = (double *)realloc(ydata,         npoints*sizeof(double));
    covar = (double *)realloc(covar, npoints*npoints*sizeof(double));
    nprops  = 0 ;
    nblocks = 0 ;
  }


  /*!
    read a list of propagators, compute average propagator and
    covariance matrix.
    Propagator is list of numbers
    x0   value
    x1   value
    x2   value
    ...
    x_max value
  */
  void ReadDataList(char *file);
  void SetDataList(double *xd, double *ll,int length ) ;
  void SetDataList(const Array<Double>& xd, const PropList& ll) ;
  void CalcCovarMat();
  void CalcJackCovarMat(); //covariance matrix for a jackknifed observable
  void WriteCovarMat(char *file);
  void WriteCovarMat(const std::string& file){
    char tt[1000] ;
    strcpy(tt,file.data()) ;
    WriteCovarMat(tt);
  }
  void ReadCovarMat(char *file);
  int Nblocks(){return nblocks;}
  int Npoints(){return npoints;}
  double& Xdata(int i){return xdata[i];}
  double& Ydata(int i){return ydata[i];}
  double& operator[](int i){return *(covar+i);}
  double& operator()(int i, int j){return *(covar+i + j*npoints);}

  ~CovarMat(){
    free(xdata);
    free(ydata);
    free(covar);
    free(ylist); 
  } ;

};

} // namespace CovFit
#endif
