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
  // FLIP 
  EventI::SetOutLeptonKinematics( temp ) ; 
  return ; 
}

void MCEvent::SetInLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  EventI::SetInLeptonKinematics(E,px,py,pz);
  TLorentzVector temp = GetInLepton4Mom() ;
  // FLIP 
  EventI::SetInLeptonKinematics( temp ) ; 

  return ; 
} 
