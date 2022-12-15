// _______________________________________________
/*
 * Event Interface base class
 * 
 */
#include <iostream>
#include "physics/EventI.h"
#include "conf/ParticleI.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"

using namespace e4nu ; 

EventI::EventI() { 
  this->Initialize() ;
}

EventI::~EventI() {
  this->Clear();
}

void EventI::SetOutLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fOutLepton.SetPxPyPzE( px, py, pz, E ) ; 
  return ; 
}

void EventI::SetInLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fInLepton.SetPxPyPzE( px, py, pz, E ) ; 
  return ; 
} 

void EventI::SetFinalParticle( const int pdg, const double E, const double px, const double py, const double pz ) {
  TLorentzVector mom;
  mom.SetPxPyPzE( px, py, pz, E ) ; 
						
  if( fFinalParticles.find(pdg) == fFinalParticles.end() ) {
    std::vector<TLorentzVector> vct ;
    vct.push_back(mom); 
    fFinalParticles.insert( std::pair<int,std::vector<TLorentzVector>>(pdg, vct) ) ; 
  } else {
    fFinalParticles[pdg].push_back( mom ) ; 
  }
}

double EventI::GetObservable( std::string observable ) {
  unsigned int target = fTargetPdg ; 
  double EBeam = GetInLepton4Mom().E();
  TLorentzVector ef4mom = GetOutLepton4Mom() ;
  TLorentzVector p4mom ; 
  bool event_wproton = false ; 
  if( fFinalParticles.find( conf::kPdgProton ) != fFinalParticles.end() ) { 
    event_wproton = true ;
    double max_mom = 0 ; 
    for ( unsigned int i = 0 ; i < fFinalParticles[conf::kPdgProton].size() ; ++i ) {
      if( p4mom.P() > max_mom ) p4mom = fFinalParticles[conf::kPdgProton][i] ;
    }
  }

  if( observable == "ECal" ) {
    if ( event_wproton == false ) return 0 ; 
    return utils::GetECal( ef4mom, p4mom, target ) ; 
  } else if ( observable == "QELRecoEnu" ) {
    return utils::GetQELRecoEnu( ef4mom, target ) ;
  } else if ( observable == "EnergyTransfer" ) {
    return utils::GetEnergyTransfer( ef4mom, EBeam ) ; 
  } else if ( observable == "RecoQ2" ) {
    return utils::GetRecoQ2( ef4mom, EBeam ) ; 
  } else if ( observable == "RecoXBJK" ) {
    return utils::GetRecoXBJK( ef4mom, EBeam ) ;
  } else if ( observable == "RecoW" ) {
    return utils::GetRecoW(ef4mom,EBeam ); 
  } else if ( observable == "DeltaAlphaT" ) {
    if ( event_wproton == false ) return 0; 
    return utils::DeltaAlphaT( ef4mom.Vect(), p4mom.Vect() , EBeam ) ; 
  } else if ( observable == "DeltaPhiT" ) {
    if ( event_wproton == false ) return 0; 
    return utils::DeltaPhiT( ef4mom.Vect(), p4mom.Vect() , EBeam ) ;
  } else if ( observable == "EFMom" ) { 
    return ef4mom.P() ;
  }
  std::cout << observable << " is NOT defined " << std::endl;
  return 0 ; 
}

void EventI::Initialize() { 
  fFinalParticles.clear() ; 
  fIsMC = false ; 
  fEventID = 0 ; 
  fWeight = 0 ; 
  fEventID = 0 ; 
  fTargetPdg = 0 ; 
  fInLeptPdg = 11 ; 
  fOutLeptPdg = 11 ; 
  fNP = 0 ; 
  fNN = 0 ; 
  fNPiP = 0 ; 
  fNPiM = 0 ; 
  fNPi0 = 0 ; 
  fNKM = 0 ; 
  fNKP = 0 ; 
  fNK0 = 0 ; 
  fNEM = 0 ; 
  fNOther = 0 ;

  fInLepton.SetPxPyPzE( 0,0,0,0 ) ;
  fOutLepton.SetPxPyPzE( 0,0,0,0 ) ;

}

void EventI::Clear() { 

  fFinalParticles.clear() ; 

}
