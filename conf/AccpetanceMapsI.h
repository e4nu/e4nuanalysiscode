/** 
 * \info This script contains general information on the acceptance maps
 **/

#include <string>
#include <iostream>

#include "conf/TargetI.h"
#include "conf/ParticleI.h"

#ifndef _ACCEPTANCEMAPS_I_H_
#define _ACCEPTANCEMAPS_I_H_

namespace e4nu {
  namespace AcceptanceMapsI { 
    
    std::string GetAcceptanceFile( const int paticle, const unsigned int target, const double E, std::string local_path ) {
      std::string base_dir = local_path + "data/AcceptanceMaps/CLAS6/";
      std::string file = base_dir + "e2a_maps_" ;
      if( target == TargetI::kPdgHe3 )	file += "3He_" ;
      else if ( target == TargetI::kPdgHe4 ) file += "4He_E_2_261.root" ;
      else if ( target == TargetI::kPdgC12 ) file += "12C_E_1_161.root";
      else return ;

      if ( E == 1.161 ) file += "1_161";
      else if ( E == 2.261 ) file+= "2_261" ; 
      else if ( E == 4.461 ) file+= "4_461" ; 
      else return ;

      if( particle == ParticleI::kPdgProton ) file += "_p.root";
      else if ( particle ==ParticleI::kPdgPiP ) file += "_pip.root" ;
      else if ( particle ==ParticleI::kPdgPiM ) file += "_pim.root" ;
      else file += ".root" ;
      
      return file ; 
    }

  }
}
#endif _ACCEPTANCEMAPS_I_H_
