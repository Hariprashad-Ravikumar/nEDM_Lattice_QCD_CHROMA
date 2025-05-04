//$Id: FitParams_io.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: FitParams_io.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:18  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.4  2007/03/03 04:02:43  kostas
//fixed spectrum to accept periodic functions
//
//Revision 1.3  2006/12/08 05:23:06  kostas
//modified spectrum to admit arbitrary function
//
//Revision 1.2  2005/07/29 14:35:16  edwards
//Changed include of files now down in "io" subdir.
//
//Revision 1.1  2005/04/27 20:00:24  edwards
//Moved from include dir. Overhauled and now in own CovFit namespace.
//
//Revision 1.3  2004/10/18 18:16:22  kostas
//converted the particle database to an xmlarray
//
//Revision 1.2  2004/05/20 16:33:35  kostas
//cvs tags
//

#include <iostream>
#include <cstdio>
#include <string>

#include "xml_array.h"
#include "io/adat_xmlio.h"
#include "covfit/proplist.h"   // this is only to get Double

#ifndef __READFITPARAMS___
#define __READFITPARAMS___

namespace CovFit {

  // Namespace composition
  using XMLArray::Array;
  using namespace ADATXML;


  enum ParticleType
  { 
    MESON,
    BARYON
  } ;
  enum SourceType
  {
    SRC_TYPE_POINT_SOURCE, 
    SRC_TYPE_WALL_SOURCE, 
    SRC_TYPE_SHELL_SOURCE, 
    SRC_TYPE_BNDST_SOURCE, 
    SRC_TYPE_POINT_AND_BNDST_SOURCE, 
    SRC_TYPE_SHELL_AND_BNDST_SOURCE, 
    SRC_TYPE_POINT_AND_SHELL_AND_BNDST_SOURCE
  };

  enum SinkType
  {
    SNK_TYPE_POINT_SINK, 
    SNK_TYPE_WALL_SINK, 
    SNK_TYPE_POINT_AND_WALL_SINK, 
    SNK_TYPE_SHELL_SINK, 
    SNK_TYPE_POINT_AND_SHELL_SINK, 
    SNK_TYPE_BNDST_SINK, 
    SNK_TYPE_POINT_AND_BNDST_SINK, 
    SNK_TYPE_SHELL_AND_BNDST_SINK, 
    SNK_TYPE_POINT_AND_SHELL_AND_BNDST_SINK
  };


  class Fitter_t{
  public:
    int Niter ;
    Double Toler ;
  };

  class CovarMat_t{
  public:
    int block ;
    bool fold ;
  } ;


  class FitterParam{
  public:
    std::string fitfunc ;
    int min_dist ;
    int max_dist ;
    Array<Double> fit_params ;
    Array<int> fixedpar ;
    Array<int> mass_param ;
    double period ;
  } ;

  void read(XMLReader& xml, const std::string& path, FitterParam& fp) ;
  void write(XMLWriter& xml, const std::string& path, const FitterParam& fp) ;

 
  class State_t{
  public:
    std::string name ;
    std::string filename ;

    SourceType   source ;
    SinkType     sink ;
    ParticleType type ;

    Array<Double> quark_mass ;
    int stopper ;
    Array<int> range ;
    FitterParam fit ;
  } ;

  class FitParams_t{
  public:
    Fitter_t          fit ;
    CovarMat_t        cov ;
    Array<State_t>    sta ;
  } ;

  void read(XMLReader& xml, const std::string& path, SourceType& param) ;
  void read(XMLReader& xml, const std::string& path, SinkType& param) ;

  void read(XMLReader& xml, const std::string& path, FitParams_t& fp) ;
  void read(XMLReader& xml, const std::string& path, Fitter_t& ft) ;
  void read(XMLReader& xml, const std::string& path, CovarMat_t& cv) ;
  void read(XMLReader& xml, const std::string& path, State_t& st) ;
  void write(XMLWriter& xml, const std::string& path, const FitParams_t& fp) ;
  void write(XMLWriter& xml, const std::string& path, const Fitter_t& ft) ;
  void write(XMLWriter& xml, const std::string& path, const CovarMat_t& cv) ;
  void write(XMLWriter& xml, const std::string& path, const State_t& st) ;

} // namespace CovFit

#endif


