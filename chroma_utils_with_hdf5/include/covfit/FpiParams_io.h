//$Id: FpiParams_io.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: FpiParams_io.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:18  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.3  2005/07/29 14:35:16  edwards
//Changed include of files now down in "io" subdir.
//
//Revision 1.2  2005/05/13 16:46:48  kostas
//added fpi
//
//Revision 1.1  2005/05/13 16:41:47  kostas
//added fpi stuff
//
//

#include <iostream>
#include <cstdio>
#include <string>

#include "io/adat_io.h"
#include "FitParams_io.h"

#ifndef __READFPIPARAMS___
#define __READFPIPARAMS___

namespace CovFit {

  // Namespace composition
  using XMLArray::Array;
  using namespace ADATXML;

  //using namespace ADATUtil;
  //using namespace std;
  //using namespace CovFit ;

  class FpiParams_t{
  public:
    Fitter_t          fit ;
    CovarMat_t        cov ;
    std::string filenameSP ;
    std::string filenameSS ;
    
    Array<Double> quark_mass ;
    int mint ;
    int maxt ;
    Array<Double> fit_paramsSP ;
    Array<Double> fit_paramsSS ;
    
    Double mres ;
    int block ;
  } ;


void read(XMLReader& xml, const std::string& path, FpiParams_t& fpi) ;
void write(XMLWriter& xml, const std::string& path, const FpiParams_t& fpi) ;

} // namespace CovFit

#endif


