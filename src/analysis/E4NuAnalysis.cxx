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
  // Loop over events
  for( unsigned int i = 0 ; i < GetNEvents() ; ++i ) {
    EventI * event = MCAnalysisI::GetValidEvent(i) ; 
    if( ! event ) continue ; 
    
    double weight = event->GetWeight() ; 
    if ( weight < 0 || weight > 10 ) continue ; 

    TLorentzVector in_mom = event -> GetInLepton4Mom() ;
    TLorentzVector out_mom = event -> GetOutLepton4Mom() ;
    double EBeam = in_mom.E() ; 

    if( out_mom.E() < conf::GetMinMomentumCut( conf::kPdgElectron, EBeam ) ) continue ; 
    
    double reco_Q2 = utils::GetRecoQ2( out_mom, EBeam ) ; 
    double W_var = utils::GetRecoW( out_mom, EBeam ) ;

    if( ApplyQ2Cut() ) {
      double MaxQ2 = 0 ; 
      if( conf::GetQ2Cut( MaxQ2, EBeam ) ) {
	if( reco_Q2 < MaxQ2 ) continue ; 
      }
    }

    if( ApplyWCut() ) {
      double MinW = 0 ; 
      if( conf::GetWCut( MinW, EBeam ) ) {
	if( W_var > MinW ) continue ; 
      }
    }

    if( ApplyThetaSlice() ) {
      if( out_mom.Theta() < conf::kMinEThetaSlice * 180./TMath::Pi() ) continue ; 
      if( out_mom.Theta() > conf::kMaxEThetaSlice * 180./TMath::Pi() ) continue ; 
    }

    if( ApplyPhiOpeningAngle() ) {
      if ( ! conf::ValidPhiOpeningAngle( out_mom.Phi() ) ) continue ;  
    }

    if( ApplyGoodSectorPhiSlice() ) {
      if ( ! conf::GoodSectorPhiSlice( out_mom.Phi() ) ) continue ; 
    }


    //    conf::ValidPhiOpeningAngle( 10. , true ) ; 
    /*
      double GetMinMomentumCut( const int particle_pdg, const double EBeam ) ;
      bool ValidPhiOpeningAngle( double, const bool apply = true ) ;
  bool GoodSectorPhiSlice( double phi , const bool apply = true ) ;
  bool GetQ2Cut( double & Q2cut, const double Ebeam, const bool apply_Q2cut = true ) ;
  bool GetWCut( double & WCut, const double Ebeam ) ;

     */
  }

  return true ; 
}
