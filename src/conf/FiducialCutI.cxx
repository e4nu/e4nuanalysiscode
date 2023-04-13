/**                                                                                                                                                                                                                   * \info This script contains general information on ficucial cut parameters
 **/

#include <iostream>
#include "conf/FiducialCutI.h"

using namespace e4nu ;

unsigned int conf::GetTorusCurrent( const double E ) {
  if ( E < 2 ) return 750 ;
  else if ( E > 2 && E < 5 ) return 2250 ;
  return 0 ; 
}

