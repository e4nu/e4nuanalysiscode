// _______________________________________________
/*
 * MC EVENT STRUCTURE
 * 
 */

#include "physics/MCEvent.h"
#include "utils/ParticleUtils.h"
#include "conf/ParticleI.h"

using namespace e4nu ; 

MCEvent::MCEvent(): EventI() { 
  fIsMC = true ; 
}

MCEvent::~MCEvent() {;}

void MCEvent::SetOutLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fOutLepton.SetPxPyPzE( px, py, pz, E ) ; 

  double phi = fOutLepton.Phi() + TMath::Pi() ; // The GENIE Coordinate system is flipped with respect to CLAS
  fOutLepton.SetPhi( phi ) ; 

  return ; 
}

void MCEvent::SetInLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fInLepton.SetPxPyPzE( px, py, pz, E ) ; 
  double phi = fInLepton.Phi() + TMath::Pi() ; // The GENIE Coordinate system is flipped with respect to CLAS
  fInLepton.SetPhi( phi ) ; 

  return ; 
} 

void MCEvent::SetFinalParticle( const int pdg, const double E, const double px, const double py, const double pz ) {
  TLorentzVector mom;
  mom.SetPxPyPzE( px, py, pz, E ) ; 

  double phi = mom.Phi() + TMath::Pi() ; // The GENIE Coordinate system is flipped with respect to CLAS                                                                                                           
  mom.SetPhi( phi ) ;
						
  if( fFinalParticles.find(pdg) == fFinalParticles.end() ) {
    std::vector<TLorentzVector> vct ;
    vct.push_back(mom); 
    fFinalParticles.insert( std::pair<int,std::vector<TLorentzVector>>(pdg, vct) ) ; 
  } else {
    fFinalParticles[pdg].push_back( mom ) ; 
  }
}

