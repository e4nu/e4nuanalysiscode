// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/CLASAnalysisI.h"

using namespace e4nu ; 

CLASAnalysisI::CLASAnalysisI() {
  this->Initialize() ;
}

CLASAnalysisI::~CLASAnalysisI() {;}

bool CLASAnalysisI::LoadData( const std::string file ) {

  return kIsDataLoaded ; 
}

bool CLASAnalysisI::LoadData( const std::string file, const unsigned int nmax ) {

  return kIsDataLoaded ; 
}

EventI * CLASAnalysisI::GetEvent( const unsigned int event_id ) {
  return nullptr ; 
}

unsigned int CLASAnalysisI::GetNEvents( void ) const {
  return 0;
}
void CLASAnalysisI::Initialize() { 

}
