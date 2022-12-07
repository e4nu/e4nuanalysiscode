// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/E4NuAnalysis.h"
#include "conf/ParticleI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AnalysisCutsI.h"

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

    //    std::cout << conf::GetAcceptanceFile( 11, GetConfiguredTarget(), in_mom.E(), "/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster" ) ; 
    std::cout << conf::kPdgElectron << std::endl;

    if( ApplyThetaSlice() ) {
      if( in_mom.Theta() < conf::kMinEThetaSlice * 180./TMath::Pi() ) continue ; 
      if( in_mom.Theta() > conf::kMaxEThetaSlice * 180./TMath::Pi() ) continue ; 
    }

    //    if( in_mom.E() < conf::GetMinMomentumCut( 11/*conf::kPdgElectron*/, in_mom.E() ) ) continue ; 

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
