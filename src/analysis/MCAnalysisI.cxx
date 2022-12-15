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
#include "conf/AccpetanceMapsI.h"

using namespace e4nu ; 

MCAnalysisI::MCAnalysisI() {
  kAcceptanceMap.clear();
  kAccMap.clear();
  kGenMap.clear();
  this->Initialize() ;
}

MCAnalysisI::~MCAnalysisI() {
  delete fData;

  kAcceptanceMap.clear();
  kAccMap.clear();
  kGenMap.clear();
}

bool MCAnalysisI::LoadData( const std::string file ) {
  double nevents = GetNEventsToRun() ; 
  double first_event = GetFirstEventToRun() ; 

  if( ! fIsDataLoaded ) { 
    fData = new MCEventHolder( file, first_event, nevents ) ;
    fIsDataLoaded = true ;
  }
  return fIsDataLoaded ; 
}

EventI * MCAnalysisI::GetEvent( const unsigned int event_id ) {
  return fData -> GetEvent(event_id) ; 
}

EventI * MCAnalysisI::GetValidEvent( const unsigned int event_id ) {
  
  MCEvent * event = (MCEvent*) fData -> GetEvent(event_id) ; 
  if( !event ) return nullptr ; 

  ++fEventsBeforeCuts ;

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
    ++fNEventsAfterFiducial ;
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
      // Apply acceptance for electron
      if( ApplyAccWeights() ) {
	if( kAccMap[it->first] && kGenMap[it->first] ) acc_wght *= utils::GetAcceptanceMapWeight( *kAccMap[it->first], *kGenMap[it->first], out_mom ) ; 
	if( fabs( acc_wght ) != acc_wght ) return nullptr ; 
      }
      is_signal = true ; 
    } else if( part_map.count( it->first ) && part_map[it->first].size() == it->second ) {
      // If we have an exclusive event, apply acceptance weight as well for signal event
      is_signal = true ; 
      for( unsigned int i = 0 ; i < it->second ; ++i ) {
	if( ApplyAccWeights() ) {
	  // Do I need particle id ? 
	  if( kAccMap[it->first] && kGenMap[it->first] ) acc_wght *= utils::GetAcceptanceMapWeight( *kAccMap[it->first], *kGenMap[it->first], part_map[it->first][i] ) ; 
	  if( fabs( acc_wght) != acc_wght ) return nullptr ; 
	}
      } 
    } else {
      is_signal = false ; 
    }

    if( is_signal && ApplyAccWeights() ) wght *= acc_wght ; 
      
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

  fData = nullptr ; 

  // Get run configurables
  double EBeam = GetConfiguredEBeam() ; 
  unsigned int Target = GetConfiguredTarget() ;

  if( ApplyFiducial() ) {
    // Initialize fiducial for this run
    kFiducialCut = std::unique_ptr<Fiducial>( new Fiducial() ) ; 
    kFiducialCut -> InitPiMinusFit( EBeam ) ; 
    kFiducialCut -> InitEClimits(); 
    kFiducialCut -> up_lim1_ec -> Eval(60) ;
    kFiducialCut -> SetConstants( conf::GetTorusCurrent( EBeam ), Target , EBeam ) ;
    kFiducialCut -> SetFiducialCutParameters( EBeam ) ;
 }

  // Initialize acceptance map histograms from file
  if( ApplyAccWeights() ) { 
    kAcceptanceMap = conf::GetAcceptanceFileMap2( Target, EBeam ) ; 

    // THESE ONE BELOW CAUSE A SMALL MEMORY LEAK - INVESTIGATE
    kAccMap[conf::kPdgElectron] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgElectron] -> Get("Accepted Particles") ) ) ;
    kAccMap[conf::kPdgProton] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgProton] -> Get("Accepted Particles") ) );
    kAccMap[conf::kPdgPiP] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiP] -> Get("Accepted Particles") ) );
    kAccMap[conf::kPdgPiM] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiM] -> Get("Accepted Particles") ) );

    kGenMap[conf::kPdgElectron] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgElectron] -> Get("Generated Particles") ) );
    kGenMap[conf::kPdgProton] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgProton] -> Get("Generated Particles") ) ) ;
    kGenMap[conf::kPdgPiP] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiP] -> Get("Generated Particles") ) ) ;
    kGenMap[conf::kPdgPiM] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiM] -> Get("Generated Particles") ) ) ;
    
    for( auto it = kAccMap.begin() ; it != kAccMap.end(); ++it ) {
      kAccMap[it->first]->SetDirectory(nullptr);
      kGenMap[it->first]->SetDirectory(nullptr);
    }
  }  

}

bool MCAnalysisI::Finalise( const std::string out_file ) {

  std::cout << " Total Number of Events Processed = " << fEventsBeforeCuts << std::endl;
  std::cout << " Events after fiducial cuts = " << fNEventsAfterFiducial << std::endl;

  return true ; 
}
