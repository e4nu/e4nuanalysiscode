/** 
 * \info This script contains general information on the acceptance maps
 **/

#include <string>
#include <iostream>
#include <cstdlib>
#include "conf/TargetI.h"
#include "conf/ParticleI.h"
#include "conf/AccpetanceMapsI.h"

using namespace e4nu ; 

std::string conf::GetAcceptanceFile( const int particle, const unsigned int target, const double E ) {

  static const  std::string local_path = std::getenv("E4NUANALYSIS");

  static const std::string base_dir = local_path + "/data/AcceptanceMaps/CLAS6/";
  std::string file = base_dir + "e2a_maps_" ;
  if( target == conf::kPdgHe3 )file += "3He_" ;
  else if ( target == conf::kPdgHe4 ) file += "4He_" ;
  else if ( target == conf::kPdgC12 ) file += "12C_";
  
  if ( E == 1.161 ) file += "E_1_161";
  else if ( E == 2.261 ) file+= "E_2_261" ; 
  else if ( E == 4.461 ) file+= "E_4_461" ; 
  
  if( particle == conf::kPdgProton ) file += "_p.root";
  else if ( particle == conf::kPdgPiP ) file += "_pip.root" ;
  else if ( particle == conf::kPdgPiM ) file += "_pim.root" ;
  else file += ".root" ;
  
  return file ; 
}

std::map<int,TFile*> conf::GetAcceptanceFileMap( const unsigned int target, const double EBeam ) {
  std::map<int,TFile*> acc_map ; 
  acc_map[conf::kPdgProton] = TFile::Open( conf::GetAcceptanceFile( conf::kPdgProton, target, EBeam ).c_str(), "READ" ) ;
  acc_map[conf::kPdgPiP] = TFile::Open( conf::GetAcceptanceFile( conf::kPdgPiP, target, EBeam ).c_str(), "READ" ) ;
  acc_map[conf::kPdgPiM] = TFile::Open( conf::GetAcceptanceFile( conf::kPdgPiM, target, EBeam ).c_str(), "READ" ) ;
  acc_map[conf::kPdgElectron] = TFile::Open( conf::GetAcceptanceFile( conf::kPdgElectron, target, EBeam ).c_str(), "READ" ) ; 
  return acc_map ; 
}
