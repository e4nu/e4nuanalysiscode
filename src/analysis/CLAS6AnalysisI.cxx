// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/CLAS6AnalysisI.h"

using namespace e4nu ; 

CLAS6AnalysisI::CLAS6AnalysisI() {
  this->Initialize() ;
}

CLAS6AnalysisI::~CLAS6AnalysisI() {;}

bool CLAS6AnalysisI::LoadData( void ) {

  return kIsDataLoaded ; 
}

EventI * CLAS6AnalysisI::GetEvent( const unsigned int event_id ) {
  return nullptr ; 
}

unsigned int CLAS6AnalysisI::GetNEvents( void ) const {
  return 0;
}
void CLAS6AnalysisI::Initialize() { 

}
