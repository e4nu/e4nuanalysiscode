// _______________________________________________
/*
 * Event Interface base class
 * 
 */
#include <iostream>
#include "physics/EventI.h"
#include "conf/ParticleI.h"
#include "utils/ParticleUtils.h"

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
