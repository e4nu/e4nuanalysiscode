/**
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date Jun 2025
 **/
#include "conf/CLAS6RunRequirements.h"

using namespace e4nu ;

std::vector<unsigned int> conf::GetInvalidRuns( const double Beam, const std::string Target ) {
  std::vector<unsigned int> invalid_runs ; 

  if( Target == "C12" ) { 
    if( Beam == 1.161 ) { invalid_runs = {18294,18297,18298,18306,18307}; }
    if( Beam == 2.261 ) { invalid_runs = {18096, 18131}; }
    if( Beam == 4.461 ) { invalid_runs = {18519} ; }
  } 

  // Add here requirements for other targets !
  return invalid_runs;
}
