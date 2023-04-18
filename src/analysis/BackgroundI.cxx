/**
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/BackgroundI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/Utils.h"

using namespace e4nu; 

BackgroundI::BackgroundI( ) {
  this->Initialize();
}

BackgroundI::BackgroundI( const std::string input_file ) : ConfigureI( input_file ) {
  this->Initialize();
}

BackgroundI::BackgroundI( const double EBeam, const unsigned int TargetPdg ) : ConfigureI( EBeam, TargetPdg ) {
  this->Initialize();
}    

BackgroundI::~BackgroundI() {
  delete kRotation ;
}

void BackgroundI::Initialize(void){
  if( kIsConfigured && ApplyFiducial() ) {
    kRotation = new Subtraction();
    kRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), GetFiducialCut() );
    kRotation->ResetQVector(); 
  }
}

unsigned int BackgroundI::GetMinParticleMultiplicity( int pdg ) const {
  std::map<int,unsigned int> Topology = GetTopology();
  return Topology[pdg] ;
}
