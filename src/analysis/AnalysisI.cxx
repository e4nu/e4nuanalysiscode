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
#include "utils/Utils.h"

using namespace e4nu; 

AnalysisI::AnalysisI( ) {;}

AnalysisI::AnalysisI( const std::string input_file ) : BackgroundI( input_file ) { ; }

AnalysisI::AnalysisI( const double EBeam, const unsigned int TargetPdg ) : BackgroundI( EBeam, TargetPdg ) {;}    

AnalysisI::~AnalysisI() {
  
}

bool AnalysisI::Analyse( EventI * event ) {

  TLorentzVector in_mom = event -> GetInLepton4Mom() ;
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ;
  double EBeam = in_mom.E() ; 

  if( out_mom.Theta() * 180 / TMath::Pi() < GetElectronMinTheta( out_mom ) ) return false ;
  ++fNEventsAfterEThetaCut ;     

  if( ApplyOutElectronCut() ){
    if( out_mom.P() < conf::GetMinMomentumCut( conf::kPdgElectron, EBeam ) ) return false ; 
  }
  ++fNEventsAfterEMomCut ;     

  if( ApplyThetaSlice() ) {
    if( out_mom.Theta() * 180./TMath::Pi() < conf::kMinEThetaSlice ) return false ; 
    if( out_mom.Theta() * 180./TMath::Pi() > conf::kMaxEThetaSlice ) return false ; 
    ++fNEventsAfterThetaCut ; 
  }

  if( ApplyPhiOpeningAngle() ) {
    double phi = out_mom.Phi() ; 
    if ( ! IsData() ) phi += TMath::Pi() ; 
    if ( ! conf::ValidPhiOpeningAngle( phi ) ) return false ;  
  }
  ++fNEventsAfterPhiOpeningAngleCut ; 

  if( ApplyGoodSectorPhiSlice() ) {
    double phi = out_mom.Phi() ; 
    if ( ! IsData() ) phi += TMath::Pi() ; 
    if ( ! conf::GoodSectorPhiSlice( phi ) ) return false ; 
  }
  ++fNEventsAfterPhiCut ; 

  double reco_Q2 = utils::GetRecoQ2( out_mom, EBeam ) ; 
  double W_var = utils::GetRecoW( out_mom, EBeam ) ;

  if( ApplyQ2Cut() ) {
    double MaxQ2 = 0 ; 
    if( conf::GetQ2Cut( MaxQ2, EBeam ) ) {
      if( reco_Q2 < MaxQ2 ) return false ; 
    }
    ++fNEventsAfterQ2Cut ; 
  }

  if( ApplyWCut() ) {
    double MinW = 0 ; 
    if( conf::GetWCut( MinW, EBeam ) ) {
      if( W_var > MinW ) return false ; 
    }
    ++fNEventsAfterWCut ; 
  }

  return true ; 
}


double AnalysisI::GetElectronMinTheta( TLorentzVector emom ) {
  return fElectronFit ->Eval(emom.P()) ; 
}

bool AnalysisI::Finalise(void) const{
  std::cout << " Events after electron momentum cut = " << fNEventsAfterEMomCut << std::endl;
  std::cout << " Events after electron theta cut = " << fNEventsAfterEThetaCut << std::endl;
  if( ApplyQ2Cut() ) std::cout << " Events after Q2 cut = " << fNEventsAfterQ2Cut << std::endl;
  if( ApplyWCut() ) std::cout << " Events after W cut = "<< fNEventsAfterWCut <<std::endl;
  if( ApplyThetaSlice() ) std::cout << " Events after theta cut = " << fNEventsAfterThetaCut << std::endl;
  if( ApplyPhiOpeningAngle() ) std::cout << " Events after phi opening angle cut = " << fNEventsAfterPhiOpeningAngleCut << std::endl;
  if( ApplyGoodSectorPhiSlice() ) std::cout << " Events after phi cut = " << fNEventsAfterPhiCut << std::endl;
  return true ; 
}
