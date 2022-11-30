// _______________________________________________
/*
 * Event Interface base class
 * 
 */
#include <iostream>
#include "physics/EventI.h"

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

void EventI::Initialize() { 

}

void EventI::Clear() { 

}
