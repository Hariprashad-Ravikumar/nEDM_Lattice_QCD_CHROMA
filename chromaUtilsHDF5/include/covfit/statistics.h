#include <cmath>
#include "covfit/proplist.h"

#ifndef STATISTICS_H
#define STATISTICS_H

namespace CovFit {

  template<typename T> T mean(const Array<T>& data){
    T me ;
    me = data[0] - data[0] ; // trick in case T is not a simple type
    for(int i(0);i<data.size();i++){
      me += data[i] ;
    }
    T norm = me ;
    norm = double(1.0/double(data.size())) ;
    me *= norm ;
    
    return me ;
  }
  
  template<typename T> T var(const Array<T>& data){
    T me = mean(data) ;
    T va ;
    T vm ;
    va = me - me; // trick in case T is not a simple type
    vm = va ; 
    for(int i(0);i<data.size();i++){
      T dd(data[i]) ;
      dd -= me ;
      vm += dd ;
      va += dd*dd ;
    }
    T norm = me ;
    norm = double(data.size()-1) ;
    va /= norm ;
    vm *= vm ;
    norm = double(data.size()*(data.size()-1)) ;
    vm /= norm ;
    va -= vm ;
    // page 607 numerical recipes : 
    //          var = (<(x-mean)^2> - <x-mean><x-mean>/N)/(N-1) 
    
    return va ;
  }
  
  template<typename T> T stdev(const Array<T>& data){
    return sqrt(var(data));
  } 
  
  template<typename T> T err(const Array<T>& data){
    return stdev(data)/sqrt(Double(data.size())) ;
  } 
  
  template<typename T> T jackerr(const Array<T>& data){
    return Double(data.size())*stdev(data)/sqrt(double(data.size()+1)); 
  }

  template<typename T> T booterr(const Array<T>& data){
    return stdev(data)/sqrt(double(data.size()-1)); 
  }

}  // namespace CovFit

#endif
