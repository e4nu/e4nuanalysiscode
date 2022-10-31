/**                                                                                                                                                                                                                   
 * \info This script contains general information on ficucial cut parameters
 **/

#include <iostream>

#ifndef _FIDUCIALCUT_I_H_
#define _FIDUCIALCUT_I_H_

namespace e4nu {
  namespace FiducialCutI {
    const unsigned int GetTorusCurrent( const unsigned double E ) {
      if ( E < 2 ) return 750 ;
      else ( E > 2 and E < 5 ) return 2250 ;
      return 0 ; 
    }
  }
}

#endif _FIDUCIALCUT_I_H
