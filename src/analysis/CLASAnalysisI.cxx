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
  //if( ! fIsDataLoaded ) fData = new CLASEventHolder(file);
  //else return fIsDataLoaded ; 

  fIsDataLoaded = true ;

  return fIsDataLoaded ; 
}

bool CLASAnalysisI::LoadData( const std::string file, const unsigned int nmax ) {
  /*
  if( ! fIsDataLoaded ) fData = new CLASEventHolder( file, nmax ) ;
  else return fIsDataLoaded ; 
  */
  fIsDataLoaded = true ;

  return fIsDataLoaded ; 
}

void CLASAnalysisI::GetEvent( const unsigned int event_id ) const {
  //  fData -> LoadEvent(event_id) ; 
}

void CLASAnalysisI::Initialize() { 

}
