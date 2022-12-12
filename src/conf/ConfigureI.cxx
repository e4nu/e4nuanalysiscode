/** 
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "conf/ConfigureI.h"

namespace e4nu { 

  namespace conf {

    ConfigureI::ConfigureI( ) {;}

    ConfigureI::ConfigureI( const double EBeam, const unsigned int TargetPdg ) { 
      kEBeam = EBeam ; 
      kTargetPdg = TargetPdg ;
      PrintConfiguration();
    }
    
    ConfigureI::~ConfigureI() {
      kTopology_map.clear() ; 
    }


    ConfigureI::ConfigureI( const std::string input_file ) {
      std::cout << " Configuring analysis from " << input_file << " ..." << std::endl;
      std::ifstream file (input_file.c_str()) ;
      
      std::vector<std::string> param, value ; 
      std::string parameter_name, line, parameter_value ; 
      if( file.is_open() ) { 
	while ( getline (file,line) ) {
	  std::istringstream line_as_stream( line ) ; 
	  line_as_stream >> parameter_name >> parameter_value ; 
	  param.push_back(parameter_name) ; 
	  value.push_back(parameter_value) ; 
	}
      }
      
      for ( unsigned int i = 0 ; i < param.size(); ++i ) {
	if ( param[i] == "UseAllSectors" ) {
	  if( value[i] == "true" ) kUseAllSectors = true ; 
	  else kUseAllSectors = false ; 
	} else if ( param[i] == "ApplyFiducial" ){
	  if( value[i] == "false" ) kApplyFiducial = false ; 
	  else kApplyFiducial = true ; 
	} else if ( param[i] == "ApplyAccWeights" ){
	  if ( value[i] == "false" ) kApplyAccWeights = false ; 
	  else kApplyAccWeights = true ; 
	} else if ( param[i] == "ApplyReso" ) {
	  if ( value[i] == "false" ) kApplyReso = false ;
	  else kApplyReso = true ; 
	} else if ( param[i] == "ApplyPhiOpeningAngle" ) {
	  if( value[i] == "true" ) kApplyPhiOpeningAngle = true ; 
	  else kApplyPhiOpeningAngle = false ; 
	} else if ( param[i] == "UsePhiThetaBand" ) {
	  if( value[i] == "true" ) kUsePhiThetaBand = true ;
	  else kUsePhiThetaBand = false ; 
	} else if ( param[i] == "ApplyThetaSlice" ) {
	  if ( value[i] == "true" ) kApplyThetaSlice = true ; 
	  else kApplyThetaSlice = false ; 
	} else if ( param[i] == "ApplyGoodSectorPhiSlice" ) { 
	  if ( value[i] == "true" ) kApplyGoodSectorPhiSlice = true ; 
	  else kApplyGoodSectorPhiSlice = false ;
	} else if ( param[i] == "IsData" ) {
	  if( value[i] == "true" ) kIsData = true ; 
	  else kIsData = false ; 
	} else if( param[i] == "ApplyQ2Cut" ) { 
	  if( value[i] == "true" ) kQ2Cut = true ; 
	  else kQ2Cut = false ; 
	} else if( param[i] == "ApplyWCut" ) { 
	  if( value[i] == "true" ) kWCut = true ; 
	  else kWCut = false ; 	
	} else if( param[i] == "IsElectronData" ) { 
	  if( value[i] == "true" ) kIsElectron = true ; 
	  else kIsElectron = false ; 
	}else if ( param[i] == "offset" ) koffset = std::stod( value[i] ) ; 
	else if ( param[i] == "EBeam" ) kEBeam = std::stod( value[i] ) ; 
	else if ( param[i] == "TargetPdg" ) kTargetPdg = (unsigned int) std::stoi( value[i] ) ; 
	else if ( param[i] == "Toplogy") {
	  std::string delimeter = ":" ;
	  std::string separator = "," ;
	  std::string element, m_element ;
	  std::istringstream particle_list( value[i] ) ;
	  while( getline( particle_list, element, ',' ) ) {
	      std::istringstream part_list_mult ( element ) ;
	      std::vector<std::string> temp ; 
	      while ( getline( part_list_mult, m_element, ':' ) ) {
		temp.push_back( m_element ) ; 
	      }
	      if( temp.size() == 2 ) {
		std::pair<int, unsigned int> pair ( std::stoi( temp[0]), (unsigned int) std::stoi( temp[1] ) ) ; // Pdg, multiplicity
		kTopology_map.insert( pair ) ;
	      }
	  }
	}
      }
      PrintConfiguration() ;
    }

    void ConfigureI::PrintConfiguration(void) const { 
      std::cout << "*********************************************************************" << std::endl;
      std::cout << "*                         E4NU ANALYSIS CONF                       **" << std::endl;
      std::cout << "*********************************************************************" << std::endl;
      std::cout << "UseAllSectors:" << kUseAllSectors << std::endl;
      std::cout << "ApplyFiducial:" << kApplyFiducial << std::endl;
      std::cout << "ApplyReso:" << kApplyReso << std::endl;
      std::cout << "ApplyQ2Cut: " << kQ2Cut << std::endl;
      std::cout << "ApplyWCut: " << kWCut << std::endl;
      std::cout << "ApplyPhiOpeningAngle:" << kApplyPhiOpeningAngle << std::endl;
      std::cout << "UsePhiThetaBand:"<< kUsePhiThetaBand << std::endl;
      std::cout << "ApplyThetaSlice:"<< kApplyThetaSlice << std::endl;
      std::cout << "ApplyGoodSectorPhiSlice:"<<kApplyGoodSectorPhiSlice<<std::endl;
      std::cout << "EBeam: " << kEBeam << " GeV " << std::endl;
      std::cout << "Target Pdg : " << kTargetPdg << std::endl;
      if ( kIsData ) std::cout << "\nIsData" << std::endl;
      else std::cout << "Is MC Data " << std::endl;
      if( kIsElectron) std::cout << " Electron scattering data" <<std::endl;
      std::cout << "Topology: \n " << std::endl;
      for ( auto it = kTopology_map.begin() ; it != kTopology_map.end() ; ++it ) { 
	std::cout << "    " << it->first << ", multiplicity " << it->second << std::endl;
      }
      std::cout << "*********************************************************************" << std::endl;
    }

  }
}
