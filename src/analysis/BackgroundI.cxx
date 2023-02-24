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

bool BackgroundI::NewBackgroundSubstraction( std::map<int,std::vector<MCEvent>> & event_holder ) { 
  if( !ApplyFiducial()  ) return true ; 
  if( !GetSubtractBkg() ) return true ;
  Fiducial * fiducial = GetFiducialCut() ; 

  unsigned int max_mult = GetMaxBkgMult(); // Max multiplicity specified in conf file
  unsigned int min_mult = GetMinBkgMult(); // Signal multiplicity
  std::map<int,unsigned int> Topology = GetTopology();
  unsigned int m = max_mult ;
  
  while ( m >= min_mult ) {
    if( event_holder.find(m) != event_holder.end() ) {
      std::cout<< " Substracting background events with with multiplicity " << m << ". The total number of bkg events is: " << event_holder[m].size() <<std::endl; 

      for( unsigned int event_id = 0 ; event_id < event_holder[m].size() ; ++event_id ) { 
	// Add counter for same multiplicity 
	unsigned int N_tot = 0 ; 
	std::map<std::map<std::vector<int>,std::vector<int>>, unsigned int> counter ; // map< map<pdg_vector, part_id vector> , vector<cout> > 
	                                                                              // size of pdg_vector is multiplicity
	                                                                              // counters is the number of events with that specific topology and id list 
	// The following mapps contain the output of the fiducial cut for each particle in an event
	std::map<int,std::vector<bool>> is_particle_contained ; 
	// Create map id
	std::map<std::vector<int>,std::vector<int>> new_topology ;
	std::vector<int> part_pdg_list, part_id_list ;	
	// Start rotations
	for ( unsigned int rot_id = 0 ; rot_id < GetNRotations() ; ++rot_id ) { 
	  // Set rotation around q3 vector
	  TVector3 VectorRecoQ = event_holder[m][event_id].GetRecoq3() ;	
	  double rotation_angle = gRandom->Uniform(0,2*TMath::Pi());
	  
	  // Rotate all Hadrons
	  std::map<int,std::vector<TLorentzVector>> rot_particles = event_holder[m][event_id].GetFinalParticles4Mom() ;
	  std::map<int,std::vector<TLorentzVector>> rot_particles_uncorr = event_holder[m][event_id].GetFinalParticlesUnCorr4Mom() ;       
	  unsigned int rot_event_mult = 0 ; // rotated event multiplicity
	  for( auto it = rot_particles.begin() ; it != rot_particles.end() ; ++it ) { 
	    for ( unsigned int part_id = 0 ; (it->second).size() ; ++part_id ) {
	      int part_pdg = it->first ; 
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
	  
	  // If multiplicity is the same as the original event multiplicity,
	  if( rot_event_mult == m ) ++N_tot ; 
      
	  // For other case, store combination of particles which contribute
	  // And add entry in corresponding map 
	  new_topology[part_pdg_list] = part_id_list ; 
	  counter[new_topology] += 1 ; 
	
	}// Close rotation loop

	// Here to add to correct multiplicity with correct event

	// Clean maps for next iteration 
	is_particle_contained.clear() ;
	new_topology.clear();
      } // Close event loop 
    }
    --m; 
  }


  return true ; 


}
