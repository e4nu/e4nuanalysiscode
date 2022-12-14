/** 
 * \info This script contains general information on the acceptance maps
 **/

#ifndef _ACCEPTANCEMAPS_I_H_
#define _ACCEPTANCEMAPS_I_H_

#include <iostream>
#include <TFile.h>
#include <map>

namespace e4nu {
  namespace conf { 
    
    std::string GetAcceptanceFile( const int particle, const unsigned int target, const double E ) ;
    std::map<int,TFile*> GetAcceptanceFileMap( const unsigned int target, const double E ) ;
  }
}
#endif 
