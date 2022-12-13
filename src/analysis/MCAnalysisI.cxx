// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/MCAnalysisI.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"
#include "conf/ParticleI.h"
#include "utils/DetectorUtils.h"
#include "conf/FiducialCutI.h"

using namespace e4nu ; 

MCAnalysisI::MCAnalysisI() {
  this->Initialize() ;
}

MCAnalysisI::~MCAnalysisI() {;}

bool MCAnalysisI::LoadData( const std::string file ) {
  if( ! fIsDataLoaded ) fData = new MCEventHolder(file);
  else return fIsDataLoaded ; 

  fIsDataLoaded = true ;

  return fIsDataLoaded ; 
}

bool MCAnalysisI::LoadData( const std::string file, const unsigned int nmax ) {

  if( ! fIsDataLoaded ) fData = new MCEventHolder( file, nmax ) ;
  else return fIsDataLoaded ; 

  fIsDataLoaded = true ;

  return fIsDataLoaded ; 
}

EventI * MCAnalysisI::GetEvent( const unsigned int event_id ) {
  return fData -> GetEvent(event_id) ; 
}

EventI * MCAnalysisI::GetValidEvent( const unsigned int event_id ) {

  MCEvent * event = (MCEvent*) fData -> GetEvent(event_id) ; 
  if( !event ) return nullptr ; 

  TLorentzVector in_mom = event -> GetInLepton4Mom() ; 
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ; 

  // Check run is correct
  double EBeam = GetConfiguredEBeam() ; 
  if ( in_mom.E() != EBeam ) return nullptr ; 
  if ( (unsigned int) event -> GetTargetPdg() != GetConfiguredTarget() ) return nullptr ; 

  // Check weight is physical
  double wght = event->GetWeight() ; 
  if ( wght < 0 || wght > 10 ) return nullptr ; 
  
  // Apply Fiducial volume cuts
  if( ApplyFiducial() ) {
    if (! kFiducialCut -> EFiducialCut(EBeam, out_mom.Vect() ) ) return nullptr ;
  }

  // Apply smaring to particles
  if( ApplyReso() ) {
    this -> SmearParticles( event ) ; 
  }

  // Get Topology Definition
  std::map<int,unsigned int> Topology = GetTopology(); 
  std::map<int,std::vector<TLorentzVector>> part_map = event -> GetFinalParticles4Mom() ;

  // Signal event
  bool is_signal = false ; 
  double acc_wght = 1 ; 

  for( auto it = Topology.begin() ; it != Topology.end() ; ++it ) {
    if( it->first == conf::kPdgElectron ) {
      // Apply acceptance for electrons
      acc_wght *= utils::GetAcceptanceMapWeight( it->first, out_mom, GetConfiguredTarget(), EBeam ) ; 
      if( fabs( acc_wght) != acc_wght ) return nullptr ; 
      is_signal = true ; 
    } else {
      // If we have an exclusive event, apply acceptance weight as well for signal event
      if( part_map.count( it->first ) && part_map[it->first].size() == it->second ) {
	is_signal = true ; 
	for( unsigned int i = 0 ; i < it->second ; ++i ) {
	  acc_wght *= utils::GetAcceptanceMapWeight( it->first, part_map[it->first][i], GetConfiguredTarget(), EBeam ) ; 
	  if( fabs( acc_wght) != acc_wght ) return nullptr ; 
	} 
      } else {
	is_signal = false ; 
      }
    }
    if( is_signal ) wght *= acc_wght ; 

    // Background event
    // BACKGROUND SUBSTRACTION METHOD HERE
  }
   
  // Apply Mottxsec weight
  //
  if( IsElectronData() ) {
    wght *= utils::GetXSecScale(out_mom, EBeam, true ) ; 
  }

  // Set Final weight
  event->SetWeight( wght ) ; 
  
  return event ; 
}

void MCAnalysisI::SmearParticles( MCEvent * event ) {
  
  TLorentzVector in_mom = event -> GetInLepton4Mom() ; 
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ; 
  std::map<int,std::vector<TLorentzVector>> part_map = event -> GetFinalParticles4Mom() ;

  double EBeam = GetConfiguredEBeam() ; 

  utils::ApplyResolution( conf::kPdgElectron , in_mom, EBeam ) ; 
  event -> EventI::SetOutLeptonKinematics( in_mom ) ; 
  
  utils::ApplyResolution( conf::kPdgElectron , out_mom, EBeam ) ; 
  event -> EventI::SetOutLeptonKinematics( out_mom ) ; 

  for( std::map<int,std::vector<TLorentzVector>>::iterator it = part_map.begin() ; it != part_map.end() ; ++it ) {
    std::vector<TLorentzVector> vtemp ; 
    for( unsigned int i = 0 ; i < (it->second).size() ; ++i ) { 
      TLorentzVector temp = (it->second)[i] ; 
      utils::ApplyResolution( it->first, temp, EBeam ) ;
      vtemp.push_back(temp) ; 
    }
    part_map[it->first] = vtemp ; 
  }
  event -> EventI::SetFinalParticlesKinematics( part_map ) ; 
} 

unsigned int MCAnalysisI::GetNEvents( void ) const {
  return (unsigned int) fData ->GetNEvents() ; 
}
void MCAnalysisI::Initialize() { 
  // Initialize fiducial for this run
  kFiducialCut = new Fiducial() ;
  double EBeam = GetConfiguredEBeam() ; 

  if( ApplyFiducial() ) {
    kFiducialCut -> InitPiMinusFit( EBeam ) ; 
    kFiducialCut -> InitEClimits(); 
    kFiducialCut -> up_lim1_ec -> Eval(60) ;
    kFiducialCut -> SetConstants( conf::GetTorusCurrent( EBeam ), GetConfiguredTarget() , EBeam ) ;
    SetFiducial( kFiducialCut -> SetFiducialCutParameters( EBeam ) ) ;
    if( ApplyFiducial() ) std::cout << " Succesfully setup Fiducial volume parameters " << std::endl;
    else std::cout << " Turning off fiducial volume cut for this run " << std::endl;
 }
}
