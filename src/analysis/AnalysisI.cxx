// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/AnalysisI.h"

using namespace e4nu ; 
using namespace e4nu::conf ; 

AnalysisI::AnalysisI() : ConfigureI() { 
  this->Initialize() ;
}

AnalysisI::AnalysisI( const std::string conf_file ) : ConfigureI( conf_file ) { this->Initialize() ; }

AnalysisI::AnalysisI( const double EBeam, const unsigned int TargetPdg ) : ConfigureI( EBeam, TargetPdg) { this->Initialize() ; }

AnalysisI::~AnalysisI() { ; } 

bool AnalysisI::Analyse(void) const {

  return true ; 
}

void AnalysisI::Initialize() { 
  fIsDataLoaded = false ;
}
