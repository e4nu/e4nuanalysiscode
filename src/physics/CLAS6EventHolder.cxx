// _______________________________________________
/*
 * CLAS6EventHolder implementation 
 */

#include "physics/CLAS6EventHolder.h"

using namespace e4nu ; 

CLAS6EventHolder::CLAS6EventHolder(): EventHolderI() { this->Initialize(); } 

CLAS6EventHolder::~CLAS6EventHolder() {
  this->Clear();
}

CLAS6EventHolder::CLAS6EventHolder( const std::string file, const int first_event, const int nmaxevents ): EventHolderI( file, first_event, nmaxevents ) { 
  this->Initialize() ; 
}

CLAS6EventHolder::CLAS6EventHolder( const std::vector<std::string> files ): EventHolderI( files ) { 
  this->Initialize() ; 
}


bool CLAS6EventHolder::LoadBranch(void) {
  if( !fEventHolderChain ) return false ; 

  return true ;
}

unsigned int CLAS6EventHolder::GetNEventsChain(void) { 
  if(!fEventHolderChain) return false ; 
 
  return fEventHolderChain->GetEntries();
}

bool CLAS6EventHolder::LoadEvent( const unsigned int event_id ) {

  return true ; 
}

bool CLAS6EventHolder::LoadAllEvents(void) {
  
  return true ; 
}

void CLAS6EventHolder::Initialize() { 
  this->LoadBranch();
}

void CLAS6EventHolder::Clear() { 
  //  b_iev = nullptr ;   
}
