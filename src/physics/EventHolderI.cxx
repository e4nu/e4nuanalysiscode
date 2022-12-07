// _______________________________________________
/*
 * EventHolder Interface base class
 * 
 */
#include <iostream>
#include "physics/EventHolderI.h"

using namespace e4nu ; 

EventHolderI::EventHolderI() { 
  this->Initialize() ;
}

EventHolderI::~EventHolderI() {
  this->Clear();
}

EventHolderI::EventHolderI( const std::string file ) { 
  this->Initialize() ; 
  std::cout<< "Loading "<< file << " ... \n" ;
  if( this->LoadMembers( file ) ) fIsConfigured = true ; 
  fMaxEvents = fEventHolderChain ->GetEntries() ;
}

EventHolderI::EventHolderI( const std::string file, const int nmaxevents ) { 
  this->Initialize() ; 
  if( this->LoadMembers( file ) ) { 
    fIsConfigured = true ; 
    if( nmaxevents > fEventHolderChain -> GetEntries() || nmaxevents < 0 ) fMaxEvents = fEventHolderChain ->GetEntries() ;
    else fMaxEvents = nmaxevents ; 

    std::cout<< "Loading "<< fMaxEvents << " from " << file << " ... \n" ;
  }
  else fIsConfigured = false ; 
}

EventHolderI::EventHolderI( const std::vector<std::string> files ) { 
  this->Initialize() ; 
  for ( unsigned int file_id = 0 ; file_id < files.size() ; ++file_id ) {
    fIsConfigured *= this->LoadMembers( files[file_id] ) ; 
  }
  if( fIsConfigured ) { 
    fMaxEvents = fEventHolderChain -> GetEntries() ; 
  }
}

bool EventHolderI::LoadMembers( const std::string file ) {
  if ( ! fEventHolderChain -> Add( file.c_str() ) ) return false ; 
  return true ; 
} 

void EventHolderI::Initialize() { 
  fEventHolderChain = new TChain("gst","e4nu_analysis") ; 
  fIsConfigured = true ; 
  fMaxEvents = -1 ; 
}

void EventHolderI::Clear() { 
  fEventHolderChain = nullptr ;
  fMaxEvents = -1 ;

}
