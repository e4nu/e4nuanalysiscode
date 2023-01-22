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

double MCEvent::GetMottXSecWeight(void) { 
  fMottXSecWght = 1./utils::GetMottXSecScale( GetOutLepton4Mom(), GetInLepton4Mom().E(), fIsEM ) ; 
  return fMottXSecWght ; 
}
