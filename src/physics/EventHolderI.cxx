// _______________________________________________
/*
 * EventHolder Interface base class
 * 
 */
#include <iostream>
#include "physics/EventHolderI.h"
#include <filesystem>
#include <sys/types.h>
#include <dirent.h>

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
    if( !LoadMembers( files[file_id] ) ) fIsConfigured=false ; 
  }
  if( fIsConfigured ) { 
    fMaxEvents = fEventHolderChain -> GetEntries() ; 
  }
}

bool EventHolderI::LoadMembers( const std::string file ) {
  if( file.find(".root") != std::string::npos) {
    std::cout << " Loading single file from " << file << std::endl;
    if ( ! fEventHolderChain -> Add( file.c_str() ) ) return false ;
  } else { 
    DIR* dirp = opendir(file.c_str());
    struct dirent * dp;
    std::cout << " Loading files from " << file << std::endl;
    int count = 0 ;
    while ((dp = readdir(dirp)) != NULL) {
      std::string s1 = dp->d_name;
      if (s1.find("gst.root") != std::string::npos) {	
	std::cout << " Adding " << (file+s1)<< std::endl;
	++count ; 
	if ( ! fEventHolderChain -> Add( (file+s1).c_str() ) ) {
	  std::cout << " Failded to add " << file+s1 << std::endl;
	  closedir(dirp);
	  return false ;
	}
      }
    }
    closedir(dirp);
  }

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
