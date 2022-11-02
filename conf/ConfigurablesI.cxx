/**
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "conf/ConfigurablesI.h"

namespace e4nu { 

  namespace conf {

    e4nuConfigurationI::e4nuConfigurationI( ) {;}
    e4nuConfigurationI::e4nuConfigurationI( std::string input_file ) {
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
	if ( param[i] == "UseAllSectors" && value[i] == "true" ) kUseAllSectors = true ; 
	else if ( param[i] == "ApplyFiducial" && value[i] == "false" ) kApplyFiducial = false ; 
	else if ( param[i] == "ApplyAccWeights" && value[i] == "false" ) kApplyAccWeights = false ; 
	else if ( param[i] == "ApplyReso" && value[i] == "false" ) kApplyReso = false ;
	else if ( param[i] == "ApplyPhiOpeningAngle" && value[i] == "true" ) kApplyPhiOpeningAngle = true ; 
	else if ( param[i] == "UsePhiThetaBand" && value[i] == "true" ) kUsePhiThetaBand = true ;
	else if ( param[i] == "ApplyThetaSlice" && value[i] == "true" ) kApplyThetaSlice = true ; 
	else if ( param[i] == "ApplyGoodSectorPhiSlice" && value[i] == "true" ) kApplyGoodSectorPhiSlice = true ; 
	else if ( param[i] == "offset" ) koffset = std::stod( value[i] ) ; 
	else if ( param[i] == "EBeam" ) kEBeam = std::stod( value[i] ) ; 
	else if ( param[i] == "TargetPdg" ) kTargetPdg = (unsigned int) std::stoi( value[i] ) ; 
	else if ( param[i] == "IsData" ) kIsData = true ; 

      }

    }
  }
}
