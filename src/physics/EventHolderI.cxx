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

EventHolderI::EventHolderI( const std::string file, const unsigned int first_event, const unsigned int nmaxevents ) { 
  this->Initialize() ; 
  if( this->LoadMembers( file ) ) { 
    fIsConfigured = true ; 
    if( nmaxevents > fEventHolderChain -> GetEntries() || nmaxevents == 0 ) fMaxEvents = fEventHolderChain ->GetEntries() ;
    else fMaxEvents = nmaxevents ; 
    fFirstEvent = first_event ;
    std::cout<< "Loading "<< fMaxEvents << " from " << file ;
    if( fFirstEvent != 0 ) std::cout << " Starting from event " << fFirstEvent ;
    std::cout << " ... \n" ;
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
  fEventHolderChain = std::unique_ptr<TChain>( new TChain("gst","e4nu_analysis") ); 
  fIsConfigured = true ; 
  fMaxEvents = -1 ; 
}

void EventHolderI::Clear() { 
  fEventHolderChain = nullptr ;
  fMaxEvents = 0 ;
  fFirstEvent = 0 ; 
}
