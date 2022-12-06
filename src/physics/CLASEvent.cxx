// _______________________________________________
/*
 * CLAS EVENT STRUCTURE
 * 
 */

#include "physics/CLASEvent.h"
#include "utils/ParticleUtils.h"
#include "conf/ParticleI.h"

using namespace e4nu ; 

CLASEvent::CLASEvent(): EventI() { 
  fIsMC = false ; 
}

CLASEvent::~CLASEvent() {

}

void CLASEvent::SetOutLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  // To IMPLEMENT

  return ; 
}

void CLASEvent::SetInLeptonKinematics( const double E, const double px, const double py, const double pz ) {

  return ; 
} 
