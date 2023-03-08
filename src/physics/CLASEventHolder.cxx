// _______________________________________________
/*
 * CLASEventHolder implementation 
 */

#include "physics/CLASEventHolder.h"

using namespace e4nu ; 

CLASEventHolder::CLASEventHolder(): EventHolderI() { this->Initialize(); } 

CLASEventHolder::~CLASEventHolder() {
  this->Clear();
}

CLASEventHolder::CLASEventHolder( const std::string file, const int first_event, const int nmaxevents ): EventHolderI( file, first_event, nmaxevents ) { 
  this->Initialize() ; 
}

CLASEventHolder::CLASEventHolder( const std::vector<std::string> files ): EventHolderI( files ) { 
  this->Initialize() ; 
}


bool CLASEventHolder::LoadBranch(void) {
  if( !fEventHolderChain ) return false ; 

  return true ;
}

unsigned int CLASEventHolder::GetNEventsChain(void) { 
  if(!fEventHolderChain) return false ; 
 
  return fEventHolderChain->GetEntries();
}

bool CLASEventHolder::LoadEvent( const unsigned int event_id ) {

  return true ; 
}

bool CLASEventHolder::LoadAllEvents(void) {
  
  return true ; 
}

void CLASEventHolder::Initialize() { 
  this->LoadBranch();
}

void CLASEventHolder::Clear() { 
  //  b_iev = nullptr ;   
}
