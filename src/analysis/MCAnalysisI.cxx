// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/MCAnalysisI.h"

using namespace e4nu ; 

MCAnalysisI::MCAnalysisI() {
  this->Initialize() ;
}

MCAnalysisI::~MCAnalysisI() {;}

bool MCAnalysisI::LoadData( const std::string file ) {
  if( ! fIsDataLoaded ) fData = new MCEventHolder(file);
  else return fIsDataLoaded ; 

  fIsDataLoaded = true ;

  return fIsDataLoaded ; 
}

bool MCAnalysisI::LoadData( const std::string file, const unsigned int nmax ) {

  if( ! fIsDataLoaded ) fData = new MCEventHolder( file, nmax ) ;
  else return fIsDataLoaded ; 

  fIsDataLoaded = true ;

  return fIsDataLoaded ; 
}

void MCAnalysisI::GetEvent( const unsigned int event_id ) const {
  fData -> LoadEvent(event_id) ; 
}

void MCAnalysisI::Initialize() { 

}
