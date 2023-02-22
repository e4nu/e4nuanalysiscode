/**
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/BackgroundI.h"
#include "conf/FiducialCutI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/Utils.h"

using namespace e4nu; 

BackgroundI::BackgroundI( ) {

  if( kIsConfigured ) kIsConfigured = InitializeFiducial() ; 

  if( kIsConfigured && ApplyFiducial() ) {
    fRotation = new Subtraction();
    fRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), fFiducialCut);
    fRotation->ResetQVector(); //Resets q vector to (0,0,0) 
  }
}

BackgroundI::BackgroundI( const std::string input_file ) : ConfigureI( input_file ) {

  if( kIsConfigured ) kIsConfigured = InitializeFiducial() ; 

  if( kIsConfigured && ApplyFiducial() ) {
    fRotation = new Subtraction();
    fRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), fFiducialCut);
    fRotation->ResetQVector(); //Resets q vector to (0,0,0)
  }
}

BackgroundI::BackgroundI( const double EBeam, const unsigned int TargetPdg ) : ConfigureI( EBeam, TargetPdg ) {

  if( kIsConfigured ) kIsConfigured = InitializeFiducial() ; 

  if( kIsConfigured && ApplyFiducial() ) {
    fRotation = new Subtraction();
    fRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), fFiducialCut);
    fRotation->ResetQVector(); //Resets q vector to (0,0,0)
  }
}    

BackgroundI::~BackgroundI() {
  delete fRotation ;
  delete fFiducialCut ;
}

bool BackgroundI::InitializeFiducial(void) {
  double EBeam = GetConfiguredEBeam() ; 
  unsigned int Target = GetConfiguredTarget() ;

  if( ApplyFiducial() ) {
    // Initialize fiducial for this run
    fFiducialCut = new Fiducial() ; 
    fFiducialCut -> InitPiMinusFit( EBeam ) ; 
    fFiducialCut -> InitEClimits(); 
    fFiducialCut -> up_lim1_ec -> Eval(60) ;
    fFiducialCut -> SetConstants( conf::GetTorusCurrent( EBeam ), Target , EBeam ) ;
    fFiducialCut -> SetFiducialCutParameters( EBeam ) ;
    if( !fFiducialCut ) return false ; 
  } else { return true ; }

  return true ; 
}


