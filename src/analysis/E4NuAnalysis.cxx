// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/E4NuAnalysis.h"
#include "conf/ParticleI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/Utils.h"

using namespace e4nu ; 

E4NuAnalysis::E4NuAnalysis() {;}

E4NuAnalysis::E4NuAnalysis( const std::string conf_file ) : AnalysisI(conf_file), MCAnalysisI(), CLASAnalysisI() {;}

E4NuAnalysis::E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) : AnalysisI(EBeam, TargetPdg), MCAnalysisI(), CLASAnalysisI() {;}

E4NuAnalysis::~E4NuAnalysis() {;}

bool E4NuAnalysis::LoadData( const std::string file ) {
  //  if( IsData() ) return CLASAnalysisI::LoadData(file);
  return MCAnalysisI::LoadData(file) ; 
}

bool E4NuAnalysis::LoadData( const std::string file, const unsigned int nmax ) {

  //if( IsData() ) return CLASAnalysisI::LoadData(file,nmax);
  return MCAnalysisI::LoadData( file, nmax ) ; 
}

EventI * E4NuAnalysis::GetValidEvent( const unsigned int event_id ) {
  //if( IsData() ) CLASAnalysisI::GetEvent( event_id ) ; 
  return MCAnalysisI::GetValidEvent( event_id ) ; 
}

unsigned int E4NuAnalysis::GetNEvents( void ) const {
  return MCAnalysisI::GetNEvents() ;
}

bool E4NuAnalysis::Analyse(void) {
  unsigned int total_nevents = GetNEvents() ;
  // Loop over events
  for( unsigned int i = 0 ; i < total_nevents ; ++i ) {
    //Print percentage
    double progress = i / (double) total_nevents ; 
    if( i==0 || i % 100000 == 0 ) utils::PrintProgressBar(progress);

    std::unique_ptr<EventI> event = std::unique_ptr<EventI>((EventI*) MCAnalysisI::GetValidEvent(i) ); 
    if( ! event ) continue ;

    TLorentzVector in_mom = event -> GetInLepton4Mom() ;
    TLorentzVector out_mom = event -> GetOutLepton4Mom() ;
    double EBeam = in_mom.E() ; 

    if( out_mom.P() < conf::GetMinMomentumCut( conf::kPdgElectron, EBeam ) ) continue ; 
    ++fNEventsAfterEMomCut ;     
    
    double reco_Q2 = utils::GetRecoQ2( out_mom, EBeam ) ; 
    double W_var = utils::GetRecoW( out_mom, EBeam ) ;

    if( ApplyQ2Cut() ) {
      double MaxQ2 = 0 ; 
      if( conf::GetQ2Cut( MaxQ2, EBeam ) ) {
	if( reco_Q2 < MaxQ2 ) continue ; 
      }
      ++fNEventsAfterQ2Cut ; 
    }

    if( ApplyWCut() ) {
      double MinW = 0 ; 
      if( conf::GetWCut( MinW, EBeam ) ) {
	if( W_var > MinW ) continue ; 
      }
      ++fNEventsAfterWCut ; 
    }

    if( ApplyThetaSlice() ) {
      if( out_mom.Theta() * 180./TMath::Pi() < conf::kMinEThetaSlice ) continue ; 
      if( out_mom.Theta() * 180./TMath::Pi() > conf::kMaxEThetaSlice ) continue ; 
      ++fNEventsAfterThetaCut ; 
    }

    if( ApplyPhiOpeningAngle() ) {
      if ( ! conf::ValidPhiOpeningAngle( out_mom.Phi() ) ) continue ;  
      ++fNEventsAfterPhiOpeningAngleCut ; 
    }

    if( ApplyGoodSectorPhiSlice() ) {
      if ( ! conf::GoodSectorPhiSlice( out_mom.Phi() ) ) continue ; 
      ++fNEventsAfterPhiCut ; 
    }

  }

  return true ; 
}

bool E4NuAnalysis::Finalise( const std::string out_file ) {
  //if( IsData() )
  bool is_ok = MCAnalysisI::Finalise( out_file ) ; 

  std::cout << " Events after electron momentum cut = " << fNEventsAfterEMomCut << std::endl;
  if( ApplyQ2Cut() ) std::cout << " Events after Q2 cut = " << fNEventsAfterQ2Cut << std::endl;
  if( ApplyWCut() ) std::cout << " Events after W cut = "<< fNEventsAfterWCut <<std::endl;
  if( ApplyThetaSlice() ) std::cout << " Events after theta cut = " << fNEventsAfterThetaCut << std::endl;
  if( ApplyPhiOpeningAngle() ) std::cout << " Events after phi opening angle cut = " << fNEventsAfterPhiOpeningAngleCut << std::endl;
  if( ApplyGoodSectorPhiSlice() ) std::cout << " Events after phi cut = " << fNEventsAfterPhiCut << std::endl;

  return is_ok ; 
}
