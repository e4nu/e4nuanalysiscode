/**                                     
 * \info This script contains general information on radiative correction constants
 **/

#ifndef _RadCONST_I_H_
#define _RadCONST_I_H_
#include "conf/TargetI.h"

namespace e4nu {
  namespace conf {
    const double tH = 0.027760 ;
    const double tHe = 0.005772;
    const double tC = 0.004183;
    const double tFe = 0.008532 ;
    double GetThickness( double tgt ) { 
      if ( tgt == conf::kPdgHe3 ) return tHe ; 
      if ( tgt == conf::kPdgC12 ) return tC ; 
      if ( tgt == conf::kPdgFe56 ) return tFe ; 
      return tH ; 
    }
  }
}
#endif
