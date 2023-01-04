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
  fIsEM = false ; 
  fIsCC = false ; 
  fIsNC = false ; 
  fIsQEL = false ;
  fIsRES = false ; 
  fIsMEC = false ; 
  fIsDIS = false ;
  
  fTrueQ2s = 0 ; 
  fTrueWs = 0 ; 
  fTruexs = 0 ; 
  fTrueys = 0 ; 
  fTrueQ2 = 0 ; 
  fTrueW = 0 ; 
  fTruex = 0 ; 
  fTruey = 0 ; 
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
  TLorentzVector mom( px, py, pz, E ) ; 

  double phi = mom.Phi() + TMath::Pi() ; // The GENIE Coordinate system is flipped with respect to CLAS                                                                                                           
  mom.SetPhi( phi ) ;
						
  if( fFinalParticles.find(pdg) == fFinalParticles.end() ) {
    std::vector<TLorentzVector> vct (0);
    vct.push_back(mom); 
    fFinalParticles.insert( std::pair<int,std::vector<TLorentzVector>>(pdg, vct) ) ; 
  } else {
    fFinalParticles[pdg].push_back( mom ) ; 
  }
}

double MCEvent::GetMottXSecWeight(void) { 
  fMottXSecWght = 1./utils::GetMottXSecScale( GetOutLepton4Mom(), GetInLepton4Mom().E(), fIsEM ) ; 
  return fMottXSecWght ; 
}
