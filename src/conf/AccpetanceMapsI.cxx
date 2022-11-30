/** 
 * \info This script contains general information on the acceptance maps
 **/

#include <string>
#include <iostream>

#include "conf/TargetI.h"
#include "conf/ParticleI.h"
#include "conf/AccpetanceMapsI.h"

using namespace e4nu ; 
    
std::string conf::GetAcceptanceFile( const int particle, const unsigned int target, const double E, const std::string local_path ) {
  std::string base_dir = local_path + "data/AcceptanceMaps/CLAS6/";
  std::string file = base_dir + "e2a_maps_" ;
  if( target == conf::kPdgHe3 )file += "3He_" ;
  else if ( target == conf::kPdgHe4 ) file += "4He_E_2_261.root" ;
  else if ( target == conf::kPdgC12 ) file += "12C_E_1_161.root";
  
  if ( E == 1.161 ) file += "1_161";
  else if ( E == 2.261 ) file+= "2_261" ; 
  else if ( E == 4.461 ) file+= "4_461" ; 
  
  if( particle == conf::kPdgProton ) file += "_p.root";
  else if ( particle == conf::kPdgPiP ) file += "_pip.root" ;
  else if ( particle == conf::kPdgPiM ) file += "_pim.root" ;
  else file += ".root" ;
  
  return file ; 
}

