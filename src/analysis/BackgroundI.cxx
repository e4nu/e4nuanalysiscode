/**
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/BackgroundI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/Utils.h"

using namespace e4nu; 

BackgroundI::BackgroundI( ) {
  if( kIsConfigured && ApplyFiducial() ) {
    kRotation = new Subtraction();
    kRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), GetFiducialCut() );
    kRotation->ResetQVector(); 
  }
}

BackgroundI::BackgroundI( const std::string input_file ) : ConfigureI( input_file ) {
  if( kIsConfigured && ApplyFiducial() ) {
    kRotation = new Subtraction();
    kRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), GetFiducialCut() );
    kRotation->ResetQVector(); 
  }
}

BackgroundI::BackgroundI( const double EBeam, const unsigned int TargetPdg ) : ConfigureI( EBeam, TargetPdg ) {
  if( kIsConfigured && ApplyFiducial() ) {
    kRotation = new Subtraction();
    kRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), GetFiducialCut() );
    kRotation->ResetQVector(); 
  }
}    

BackgroundI::~BackgroundI() {
  delete kRotation ;
}

bool BackgroundI::ChechMinParticleMultiplicity( int pdg, unsigned int part_mult ) const {
  std::map<int,unsigned int> Topology = GetTopology();
  if( Topology.find(pdg) != Topology.end() ) { 
    if ( Topology[pdg] == 0 ) return true ; 
    if ( part_mult < Topology[pdg] ) return false ; 
  }
  return true ; 
}

bool BackgroundI::NewBackgroundSubstraction( std::map<int,std::vector<MCEvent>> & event_holder ) { 
  if( !ApplyFiducial()  ) return true ; 
  if( !GetSubtractBkg() ) return true ;
  Fiducial * fiducial = GetFiducialCut() ; 

  unsigned int max_mult = GetMaxBkgMult(); // Max multiplicity specified in conf file
  unsigned int min_mult = GetMinBkgMult(); // Signal multiplicity
  std::map<int,unsigned int> Topology = GetTopology();
  unsigned int m = max_mult ;
  
  while ( m > min_mult ) {
    if( event_holder.find(m) != event_holder.end() ) {
      std::cout<< " Substracting background events with with multiplicity " << m << ". The total number of bkg events is: " << event_holder[m].size() <<std::endl; 

      for( unsigned int event_id = 0 ; event_id < event_holder[m].size() ; ++event_id ) { 
	// Add counter for same multiplicity 
	double N_all = 0 ; 
	std::map<std::map<std::vector<int>,std::vector<int>>, double> probability_count ; // size of pdg_vector is multiplicity
                                                                                          // probability_counts is the number of events with that specific topology and id list 
	// Create map id
	std::map<std::vector<int>,std::vector<int>> new_topology ;

	// Start rotations
	for ( unsigned int rot_id = 0 ; rot_id < GetNRotations() ; ++rot_id ) { 
	  // Set rotation around q3 vector
	  TVector3 VectorRecoQ = event_holder[m][event_id].GetRecoq3() ;	
	  double rotation_angle = gRandom->Uniform(0,2*TMath::Pi());
	  
	  // The following mapps contain the output of the fiducial cut for each particle in an event
	  std::map<int,std::vector<bool>> is_particle_contained ; 

	  // Rotate all Hadrons
	  std::map<int,std::vector<TLorentzVector>> rot_particles = event_holder[m][event_id].GetFinalParticles4Mom() ;
	  std::map<int,std::vector<TLorentzVector>> rot_particles_uncorr = event_holder[m][event_id].GetFinalParticlesUnCorr4Mom() ;       
	  unsigned int rot_event_mult = 0 ; // rotated event multiplicity
	  std::vector<int> part_pdg_list, part_id_list ; 
	  for( auto it = rot_particles.begin() ; it != rot_particles.end() ; ++it ) {
	    int part_pdg = it->first ; 
	    if( Topology.find( part_pdg ) == Topology.end() ) continue ; // Skip particles which are not in signal definition 
	 
	    for ( unsigned int part_id = 0 ; part_id < (it->second).size() ; ++part_id ) {
	      TVector3 part_vect = (it->second)[part_id].Vect() ;
	      part_vect.Rotate(rotation_angle,VectorRecoQ);
	      
	      // Check which particles are in fiducial	      
	      if ( is_particle_contained.find(part_pdg) != is_particle_contained.end() ) {
                is_particle_contained[part_pdg].push_back( fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_vect ) ) ; 
              } else {
		std::vector<bool> temp = { fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_vect ) } ; 
                is_particle_contained[part_pdg] = temp ;
              }
	      
	      // Calculate rotated event multiplicity
	      if( is_particle_contained[part_pdg][part_id] ) {
		++rot_event_mult ; 
		part_pdg_list.push_back( part_pdg ) ; 
		part_id_list.push_back( part_id ) ; 
	      } 
	    }
	  }

	  // If multiplicity < minimum multiplicity, remove
	  if( rot_event_mult < min_mult ) continue ; 
	  
	  // Check if particle multiplicity is above signal particle multiplicity 
	  bool is_signal_bkg = true ; 
	  for( unsigned int k = 0 ; k < part_pdg_list.size() ; ++k ) { 
	    unsigned int count_particle = 0 ; 
	    for( unsigned int l = 0 ; l < is_particle_contained[part_pdg_list[k]].size() ; ++l ) { 
	      if( is_particle_contained[part_pdg_list[k]][l] ) ++count_particle ; 
	    }
	 
	    is_signal_bkg *= ChechMinParticleMultiplicity( part_pdg_list[k], count_particle ) ;
	  }

	  if( ! is_signal_bkg ) {
	    continue ; 
	  }

	  // If multiplicity is the same as the original event multiplicity,
	  if( rot_event_mult == m ) {
	    ++N_all ; 
	    continue ; 
	  }

	  // For other case, store combination of particles which contribute
	  // And add entry in corresponding map 
	  new_topology[part_pdg_list] = part_id_list ; 
	  probability_count[new_topology] += 1 ; 

	  // Clean maps for next iteration 
	  is_particle_contained.clear() ;
	}// Close rotation loop

	// Skip if denominator is 0
	if( N_all == 0 ) continue ; 

	// Store event particles with correct weight and multiplicty
	double event_wgt = event_holder[m][event_id].GetEventWeight() ;
	std::map<int,std::vector<TLorentzVector>> particles = event_holder[m][event_id].GetFinalParticles4Mom() ;
	std::map<int,std::vector<TLorentzVector>> particles_uncorr = event_holder[m][event_id].GetFinalParticlesUnCorr4Mom() ;
	for( auto it = probability_count.begin() ; it != probability_count.end() ; ++it ) { 
	  for( auto it_key = it->first.begin() ; it_key != it->first.end() ; ++it_key ) { 
	    std::map<int,std::vector<TLorentzVector>> temp_corr_mom ;
	    std::map<int,std::vector<TLorentzVector>> temp_uncorr_mom ;
	    int new_multiplicity = 0 ; 
	    for( unsigned int k = 0 ; k < (it_key->first).size() ; ++k ) { 
	      int particle_pdg = (it_key->first)[k] ; 
	      int particle_id = (it_key->second)[k] ; 
	      
	      temp_corr_mom[particle_pdg].push_back( particles[particle_pdg][particle_id] ) ; 
	      temp_uncorr_mom[particle_pdg].push_back( particles_uncorr[particle_pdg][particle_id] ) ; 
	      ++new_multiplicity ; 
	    }

	    event_holder[m][event_id].SetFinalParticlesKinematics( temp_corr_mom ) ; 
	    event_holder[m][event_id].SetFinalParticlesUnCorrKinematics( temp_uncorr_mom ) ; 

	    double probability = - (it->second) * event_wgt / N_all ; 
	    event_holder[m][event_id].SetEventWeight( probability ) ; 

	    if ( event_holder.find(new_multiplicity) != event_holder.end() ) {
	      event_holder[new_multiplicity].push_back( event_holder[m][event_id] ) ; 
	    } else {
	      std::vector<e4nu::MCEvent> temp = { event_holder[m][event_id] } ;
	      event_holder[new_multiplicity] = temp ; 
	    }
	  }
	}
	
	new_topology.clear();
      } // Close event loop 
    }
    --m; 
  }

  delete fiducial ; 

  return true ; 
}


bool BackgroundI::AcceptanceCorrection( std::map<int,std::vector<MCEvent>> & event_holder ) { 

  // We need to correct for signal events that are reconstructed outside of the fiducial
  if( !ApplyFiducial()  ) return true ; 
  if( !GetSubtractBkg() ) return true ;
  Fiducial * fiducial = GetFiducialCut() ; 

  unsigned int min_mult = GetMinBkgMult(); // Signal multiplicity
  std::map<int,unsigned int> Topology = GetTopology();
  std::vector<MCEvent> signal_events = event_holder[min_mult] ; 

  unsigned int n_truesignal = signal_events.size() ;
  for( unsigned int i = 0 ; i < n_truesignal ; ++i ) { 

    long N_signal_detected = 0 ; 
    long N_signal_undetected = 0 ; 
    
    // Start rotations
    for ( unsigned int rot_id = 0 ; rot_id < GetNRotations() ; ++rot_id ) { 
      // Set rotation around q3 vector
      TVector3 VectorRecoQ = signal_events[i].GetRecoq3() ;	
      double rotation_angle = gRandom->Uniform(0,2*TMath::Pi());
	  
      std::map<int,std::vector<TLorentzVector>> rot_particles = signal_events[i].GetFinalParticles4Mom() ;
      std::map<int,std::vector<TLorentzVector>> rot_particles_uncorr = signal_events[i].GetFinalParticlesUnCorr4Mom() ;
      
      // Rotate all particles 
      for( auto it = rot_particles.begin() ; it != rot_particles.end() ; ++it ) {
	int part_pdg = it->first ; 
	bool is_contained = true ; 
	if( Topology.find( part_pdg ) == Topology.end() ) continue ; // Skip particles which are not in signal definition 
	
	for ( unsigned int part_id = 0 ; part_id < (it->second).size() ; ++part_id ) {
	  TVector3 part_vect = (it->second)[part_id].Vect() ;
	  part_vect.Rotate(rotation_angle,VectorRecoQ);
	      
	  // Check which particles are in fiducial	      
	  is_contained *= fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_vect ) ;
	}
	if( is_contained ) ++N_signal_detected ; 
	else ++N_signal_undetected ; 
      }
    }

      // Add missing signal events
      MCEvent temp_event = signal_events[i] ;
      double event_wgt = temp_event.GetEventWeight() ;
      temp_event.SetEventWeight( + event_wgt * N_signal_undetected / N_signal_detected ) ; 
      signal_events.push_back(temp_event);

  }
  // Store correction
  event_holder[min_mult] = signal_events ; 
  
  delete fiducial ;
  return true ; 
}
