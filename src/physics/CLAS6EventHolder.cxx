// _______________________________________________
/*
 * CLAS6EventHolder implementation 
 */

#include "physics/CLAS6EventHolder.h"
#include "physics/CLAS6Event.h"

using namespace e4nu ; 

CLAS6EventHolder::CLAS6EventHolder(): EventHolderI() { this->Initialize(); } 

CLAS6EventHolder::~CLAS6EventHolder() {
  //  this->Clear();
}

CLAS6EventHolder::CLAS6EventHolder( const std::string file, const unsigned first_event, const unsigned int nmaxevents ): EventHolderI( file, first_event , nmaxevents ) { 
  this->Initialize() ; 
}

CLAS6EventHolder::CLAS6EventHolder( const std::vector<std::string> files ): EventHolderI( files ) { 
  this->Initialize() ; 
}


bool CLAS6EventHolder::LoadBranch(void) {
  if( !fEventHolderChain ) return false ; 

  //  fEventHolderChain->SetBranchAddress("iev", &iev, &b_iev);
  
  return true ;
}

EventI * CLAS6EventHolder::GetEvent(const unsigned int event_id) {

  if ( event_id > (unsigned int) fMaxEvents ) return nullptr ; 

  //  static 
  CLAS6Event * event = new CLAS6Event() ; 
  fEventHolderChain->GetEntry( event_id ) ; 
  
  //event -> SetEventID( iev ) ;
  
  return event ; 
}

void CLAS6EventHolder::Initialize() { 
  this->LoadBranch();
}

void CLAS6EventHolder::Clear() { 

}
