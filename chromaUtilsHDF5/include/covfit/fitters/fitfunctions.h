//$Id: fitfunctions.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: fitfunctions.h,v $
//Revision 2.0  2008/12/05 04:43:33  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.11  2007/05/10 02:16:06  kostas
//added the fit function I need for periodic two meson spectroscopy
//
//Revision 1.10  2007/03/08 05:15:20  kostas
//fixed the problems alternating exps. had with spectrum
//
//Revision 1.9  2007/03/03 04:02:43  kostas
//fixed spectrum to accept periodic functions
//
//Revision 1.8  2006/12/27 04:23:31  kostas
//Added periodic alternating exponential fit function
//
//Revision 1.7  2006/12/06 04:17:58  kostas
//added function "factory" ... Robert sorry for the abuse of terminology...
//
//Revision 1.6  2006/08/02 21:47:36  kostas
//added alternating exponential fits
//
//Revision 1.5  2006/07/12 00:20:26  kostas
//added some things about gA fits.
//
//Revision 1.4  2006/05/17 21:14:18  kostas
//added genneric chiPT log fit function
//
//Revision 1.3  2006/04/24 20:56:08  kostas
//Added the strange functions that apear in HBchiPT
//
//Revision 1.2  2006/03/15 19:55:07  kostas
//added some fit fitfunctions
//
//Revision 1.1  2005/11/23 03:18:51  kostas
//moved the fitter header files into the main include tree
//
//Revision 1.5  2005/10/03 16:04:02  kostas
//fixed bugs
//
//Revision 1.4  2005/03/04 15:03:43  kostas
//added anti-periodic exponentials
//
//Revision 1.3  2004/09/26 19:23:58  kostas
//added antiperiodic exponential fit for baryons (not sure if it works)
//
//Revision 1.2  2004/04/29 21:06:26  kostas
//Added periodic exponetials
//
//Revision 1.1  2004/03/17 22:30:51  kostas
//modified how to compile fitfunctions
//
/* fit  function */

#include <string>
#include <iostream>

#include "polynomium.h"
#include "expon.h"
#include "expon_altern.h"
#include "exponper_altern.h"
#include "exponper_altern_onemass.h"
#include "exponper.h"
#include "exponper_subtr.h"
#include "exponaper.h"
#include "exponaper_bar.h"
#include "chiPTfp.h"
#include "chiPTlog.h"
//#include "Lookup.h"
#include "BaryonFfunc.h"
#include "BaryonSfunc.h"
#include "gAJfunc.h"
#include "gAKfunc.h"
#include "boundmom.h"

namespace CovFit {
  Function* CreateFunction(std::string& name, int N) ;
  Function* CreateFunction(std::string& name, double p, int N) ;
}


