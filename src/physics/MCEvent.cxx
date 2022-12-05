// _______________________________________________
/*
 * MC EVENT STRUCTURE
 * 
 */

#include "physics/MCEvent.h"

using namespace e4nu ; 

MCEvent::MCEvent(): EventI() { 
  fIsMC = true ; 
}

MCEvent::~MCEvent() {

}

void MCEvent::SetOutLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  
  EventI::SetOutLeptonKinematics(E,px,py,pz);

  TLorentzVector temp = GetOutLepton4Mom() ;
  double phi = temp.Phi() + TMath::Pi() ; // The GENIE Coordinate system is flipped with respect to CLAS
  temp.SetPhi( phi ) ; 
  EventI::SetOutLeptonKinematics( temp ) ; 

  return ; 
}

void MCEvent::SetInLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  EventI::SetInLeptonKinematics(E,px,py,pz);
  
  TLorentzVector temp = GetInLepton4Mom() ;
  double phi = temp.Phi() + TMath::Pi() ; // The GENIE Coordinate system is flipped with respect to CLAS
  temp.SetPhi( phi ) ; 
  EventI::SetInLeptonKinematics( temp ) ; 

  return ; 
} 
