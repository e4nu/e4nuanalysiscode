// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/E4NuAnalysis.h"

using namespace e4nu ; 

E4NuAnalysis::E4NuAnalysis() {;}

E4NuAnalysis::E4NuAnalysis( const std::string conf_file ) : AnalysisI(conf_file), MCAnalysisI(), CLASAnalysisI() {;}

E4NuAnalysis::E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) : AnalysisI(EBeam, TargetPdg), MCAnalysisI(), CLASAnalysisI() {;}

E4NuAnalysis::~E4NuAnalysis() {;}

bool E4NuAnalysis::LoadData( const std::string file ) {
  if( IsData() ) return CLASAnalysisI::LoadData(file);
  return MCAnalysisI::LoadData(file) ; 
}

bool E4NuAnalysis::LoadData( const std::string file, const unsigned int nmax ) {

  if( IsData() ) return CLASAnalysisI::LoadData(file,nmax);
  return MCAnalysisI::LoadData( file, nmax ) ; 
}

void E4NuAnalysis::GetEvent( const unsigned int event_id ) const {
  if( IsData() ) CLASAnalysisI::GetEvent( event_id ) ; 
  return MCAnalysisI::GetEvent( event_id ) ; 
}
