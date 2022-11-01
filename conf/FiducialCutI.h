/**                                                                                                                                                                                                                   
 * \info This script contains general information on ficucial cut parameters
 **/

#include <iostream>

#ifndef _FIDUCIALCUT_I_H_
#define _FIDUCIALCUT_I_H_

namespace e4nu {
  namespace FiducialCutI {
    const unsigned int GetTorusCurrent( const double E ) {
      if ( E < 2 ) return 750 ;
      else if ( E > 2 && E < 5 ) return 2250 ;
      return 0 ; 
    }
  }
}

#endif
