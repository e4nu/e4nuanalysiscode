// _______________________________________________
/*
 * CLAS6 EVENT STRUCTURE
 * 
 */

#include "physics/CLAS6Event.h"
#include "utils/ParticleUtils.h"
#include "conf/ParticleI.h"

using namespace e4nu ; 

CLAS6Event::CLAS6Event(): EventI() { 
  fIsMC = false ; 
}

CLAS6Event::~CLAS6Event() {

}

void CLAS6Event::SetOutLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  // To IMPLEMENT

  return ; 
}

void CLAS6Event::SetInLeptonKinematics( const double E, const double px, const double py, const double pz ) {

  return ; 
} 
