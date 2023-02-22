/** 
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/AnalysisI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/DetectorUtils.h"
#include "utils/Utils.h"

using namespace e4nu; 

AnalysisI::AnalysisI( ) { ; }

AnalysisI::AnalysisI( const std::string input_file ) : BackgroundI( input_file ) { if( kIsConfigured ) kIsConfigured = InitializeFiducial() ;  }

AnalysisI::AnalysisI( const double EBeam, const unsigned int TargetPdg ) : BackgroundI( EBeam, TargetPdg ) {;}    

AnalysisI::~AnalysisI() {;}

bool AnalysisI::Analyse( EventI * event ) {

  TLorentzVector in_mom = event -> GetInLepton4Mom() ;
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ;
  double EBeam = in_mom.E() ; 
  
  if( !utils::IsValidSector( out_mom.Phi(), EBeam, UseAllSectors() ) ) return false ;   

  if( out_mom.Theta() * 180 / TMath::Pi() < GetElectronMinTheta( out_mom ) ) return false ;
  ++kNEventsAfterEThetaCut ;     

  if( ApplyOutElectronCut() ){
    if( out_mom.P() < conf::GetMinMomentumCut( conf::kPdgElectron, EBeam ) ) return false ; 
  ++kNEventsAfterEMomCut ;     
  }

  if( ApplyThetaSlice() ) {
    if( out_mom.Theta() * 180./TMath::Pi() < conf::kMinEThetaSlice ) return false ; 
    if( out_mom.Theta() * 180./TMath::Pi() > conf::kMaxEThetaSlice ) return false ; 
    ++kNEventsAfterThetaCut ; 
  }

  if( ApplyPhiOpeningAngle() ) {
    double phi = out_mom.Phi() ; 
    if ( ! IsData() ) phi += TMath::Pi() ; 
    if ( ! conf::ValidPhiOpeningAngle( phi ) ) return false ;  
  }
  ++kNEventsAfterPhiOpeningAngleCut ; 

  if( ApplyGoodSectorPhiSlice() ) {
    double phi = out_mom.Phi() ; 
    if ( ! IsData() ) phi += TMath::Pi() ; 
    if ( ! conf::GoodSectorPhiSlice( phi ) ) return false ; 
  }
  ++kNEventsAfterPhiCut ; 

  double reco_Q2 = utils::GetRecoQ2( out_mom, EBeam ) ; 
  double W_var = utils::GetRecoW( out_mom, EBeam ) ;

  if( ApplyQ2Cut() ) {
    double MaxQ2 = 0 ; 
    if( conf::GetQ2Cut( MaxQ2, EBeam ) ) {
      if( reco_Q2 < MaxQ2 ) return false ; 
    }
    ++kNEventsAfterQ2Cut ; 
  }

  if( ApplyWCut() ) {
    double MinW = 0 ; 
    if( conf::GetWCut( MinW, EBeam ) ) {
      if( W_var > MinW ) return false ; 
    }
    ++kNEventsAfterWCut ; 
  }

  return true ; 
}


double AnalysisI::GetElectronMinTheta( TLorentzVector emom ) {
  return kElectronFit ->Eval(emom.P()) ; 
}

bool AnalysisI::Finalise(void) const{
  std::cout << " Events after electron momentum cut = " << kNEventsAfterEMomCut << std::endl;
  std::cout << " Events after electron theta cut = " << kNEventsAfterEThetaCut << std::endl;
  if( ApplyQ2Cut() ) std::cout << " Events after Q2 cut = " << kNEventsAfterQ2Cut << std::endl;
  if( ApplyWCut() ) std::cout << " Events after W cut = "<< kNEventsAfterWCut <<std::endl;
  if( ApplyThetaSlice() ) std::cout << " Events after theta cut = " << kNEventsAfterThetaCut << std::endl;
  if( ApplyPhiOpeningAngle() ) std::cout << " Events after phi opening angle cut = " << kNEventsAfterPhiOpeningAngleCut << std::endl;
  if( ApplyGoodSectorPhiSlice() ) std::cout << " Events after phi cut = " << kNEventsAfterPhiCut << std::endl;
  return true ; 
}
