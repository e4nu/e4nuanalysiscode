/**
 * \info These parameters are configurable 
 * the default values are set here
 **/

#ifndef _BACKGROUND_I_H_
#define _BACKGROUND_I_H_

#include <vector>
#include <map>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "analysis/ConfigureI.h"
#include "physics/EventI.h"
#include "physics/MCEvent.h"
#include "utils/Subtraction.h"

namespace e4nu { 

  class BackgroundI : public ConfigureI {

  public: 
    // Default constructor
    BackgroundI() ;
    BackgroundI( const std::string input_file ) ;
    BackgroundI( const double EBeam, const unsigned int TargetPdg ) ;

    unsigned int GetMinParticleMultiplicity( int pdg ) const ;
    
    // These template class guarantees the same substraction method for data and MC
    // Definition must be in header file to avoid linking issuess

    template <class T>
      bool NewBackgroundSubstraction( std::map<int,std::vector<T*>> & event_holder ) { 
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
	    // Start rotations
	    for ( unsigned int rot_id = 0 ; rot_id < GetNRotations() ; ++rot_id ) { 
	      // Set rotation around q3 vector
	      TVector3 VectorRecoQ = event_holder[m][event_id]->GetRecoq3() ;	
	      double rotation_angle = gRandom->Uniform(0,2*TMath::Pi());
	  
	      // Rotate all Hadrons
	      std::map<int,std::vector<TLorentzVector>> rot_particles = event_holder[m][event_id]->GetFinalParticles4Mom() ;
	      std::map<int,std::vector<TLorentzVector>> rot_particles_uncorr = event_holder[m][event_id]->GetFinalParticlesUnCorr4Mom() ;       
	      unsigned int rot_event_mult = 0 ; // rotated event multiplicity
	      std::vector<int> part_pdg_list, part_id_list ; 

	      for( auto it = rot_particles.begin() ; it != rot_particles.end() ; ++it ) {
		int part_pdg = it->first ; 
		if( Topology.find( part_pdg ) == Topology.end() ) continue ; // Skip particles which are not in signal definition 

		for ( unsigned int part_id = 0 ; part_id < (it->second).size() ; ++part_id ) {
		  TVector3 part_vect = (it->second)[part_id].Vect() ;
		  part_vect.Rotate(rotation_angle,VectorRecoQ);
	      
		  // Check which particles are in fiducial
		  bool is_particle_contained = fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_vect ) ; 
	      
		  // Calculate rotated event multiplicity
		  if( is_particle_contained ) {
		    ++rot_event_mult ; 
		    part_pdg_list.push_back( part_pdg ) ; 
		    part_id_list.push_back( part_id ) ; 
		  } 
		} // end loop over a specific pdg
	      } // end pdg key loop

	      // If multiplicity < minimum multiplicity, remove
	      if( rot_event_mult < min_mult ) {
		continue ; 
	      }

	      // Check if particle multiplicity is above signal particle multiplicity 
	      bool is_signal_bkg = true ; 
	      for( auto it = Topology.begin(); it!=Topology.end();++it){
		if( it->first == conf::kPdgElectron ) continue ; 
		if( it->second == 0 ) continue ; // If min mult is 0, continue
		unsigned int count_p = 0 ; 
		for( unsigned int idp = 0 ; idp < part_pdg_list.size() ; ++idp ){
		  if( part_pdg_list[idp] == it->first ) ++count_p ;  
		}
		if( count_p < Topology[it->first] ) {
		  is_signal_bkg =false ;
		  break; 
		}
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
	      std::map<std::vector<int>,std::vector<int>> new_topology ;
	      new_topology[part_pdg_list] = part_id_list ; 
	      probability_count[new_topology] += 1 ; 
	 
	    }// Close rotation loop

	    // Skip if denominator is 0
	    if( N_all == 0 ) continue ; 

	    // Store event particles with correct weight and multiplicty
	    T * temp_event = event_holder[m][event_id] ;
	    double event_wgt = temp_event->GetEventWeight() ;
	    std::map<int,std::vector<TLorentzVector>> particles = temp_event->GetFinalParticles4Mom() ;
	    std::map<int,std::vector<TLorentzVector>> particles_uncorr = temp_event->GetFinalParticlesUnCorr4Mom() ;
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

		temp_event->SetFinalParticlesKinematics( temp_corr_mom ) ; 
		temp_event->SetFinalParticlesUnCorrKinematics( temp_uncorr_mom ) ; 

		double probability = - (it->second) * event_wgt / N_all ; 
		temp_event->SetEventWeight( probability ) ; 

		if ( event_holder.find(new_multiplicity) != event_holder.end() ) {
		  event_holder[new_multiplicity].push_back( temp_event ) ; 
		} else {
		  std::vector<T*> temp = { temp_event } ;
		  event_holder[new_multiplicity] = temp ; 
		}
	      }
	    }
	
	  } // Close event loop 
	}
	--m; 
      }

      return true ; 
    }

    template <class T>
      bool AcceptanceCorrection( std::map<int,std::vector<T*>> & event_holder ) { 

      // We need to correct for signal events that are reconstructed outside of the fiducial
      if( !ApplyFiducial()  ) return true ; 
      if( !GetSubtractBkg() ) return true ;
      Fiducial * fiducial = GetFiducialCut() ; 

      unsigned int min_mult = GetMinBkgMult(); // Signal multiplicity
      std::map<int,unsigned int> Topology = GetTopology();
      std::vector<T*> signal_events = event_holder[min_mult] ; 

      unsigned int n_truesignal = signal_events.size() ;
      for( unsigned int i = 0 ; i < n_truesignal ; ++i ) { 

	long N_signal_detected = 0 ; 
	long N_signal_undetected = 0 ; 
    
	// Start rotations
	for ( unsigned int rot_id = 0 ; rot_id < GetNRotations() ; ++rot_id ) { 
	  // Set rotation around q3 vector
	  TVector3 VectorRecoQ = signal_events[i]->GetRecoq3() ;	
	  double rotation_angle = gRandom->Uniform(0,2*TMath::Pi());
	  
	  std::map<int,std::vector<TLorentzVector>> rot_particles = signal_events[i]->GetFinalParticles4Mom() ;
	  std::map<int,std::vector<TLorentzVector>> rot_particles_uncorr = signal_events[i]->GetFinalParticlesUnCorr4Mom() ;
      
	  // Rotate all particles 
	  bool is_contained = true ; 
	  for( auto it = rot_particles.begin() ; it != rot_particles.end() ; ++it ) {
	    int part_pdg = it->first ; 
	    if( Topology.find( part_pdg ) == Topology.end() ) continue ; // Skip particles which are not in signal definition 
	
	    for ( unsigned int part_id = 0 ; part_id < (it->second).size() ; ++part_id ) {
	      TVector3 part_vect = (it->second)[part_id].Vect() ;
	      part_vect.Rotate(rotation_angle,VectorRecoQ);
	      
	      // Check which particles are in fiducial	      
	      is_contained *= fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_vect ) ;
	    }
	  }
	  if( is_contained ) ++N_signal_detected ; 
	  else ++N_signal_undetected ; 
	}
	if( N_signal_detected == 0 ) continue ; 
	// Add missing signal events
	T * temp_event = signal_events[i] ;
	double event_wgt = temp_event->GetEventWeight() ;
	temp_event->SetEventWeight( + event_wgt * N_signal_undetected / N_signal_detected ) ; 
	temp_event->SetIsUndetectedSignal(true);
	signal_events.push_back(temp_event);
      }
      // Store correction
      event_holder[min_mult] = signal_events ; 
  
      return true ; 
    }



    template <class T> 
      bool SubtractBackground( std::map<int,std::vector<T>> & event_holder ) {
      if( !ApplyFiducial()  ) return true ; 
      if( !GetSubtractBkg() ) return true ;

      unsigned int max_mult = GetMaxBkgMult(); // Max multiplicity specified in conf file
      unsigned int min_mult = GetMinBkgMult(); // Signal multiplicity
      std::map<int,unsigned int> Topology = GetTopology();
      unsigned int m = max_mult ;
      while ( m >= min_mult ) {
	if( event_holder.find(m) != event_holder.end() ) {
	  std::cout<< " Number of events with multiplicity " << m << " = " << event_holder[m].size() <<std::endl; 
	}
	--m; 
      }

      std::map<int,std::vector<TLorentzVector>> particles ;
      std::map<int,std::vector<TLorentzVector>> particles_uncorr ; 
      TLorentzVector V4_el ;
      unsigned int bkg_mult = min_mult + 1 ;
      if( bkg_mult > max_mult ) return true ; 

      // Create new map for event particles
      std::map<int,std::vector<TLorentzVector>> t_particles ;
      std::map<int,std::vector<TLorentzVector>> t_particles_uncorr ;
      if ( event_holder.find(bkg_mult) != event_holder.end() ) {
	for ( unsigned int i = 0 ; i < event_holder[bkg_mult].size() ; ++i ) {
	  particles = event_holder[bkg_mult][i].GetFinalParticles4Mom();
	  particles_uncorr = event_holder[bkg_mult][i].GetFinalParticlesUnCorr4Mom(); // This map needs to change with the cuts as well...
	  V4_el = event_holder[bkg_mult][i].GetOutLepton4Mom();
	  double event_wgt = event_holder[bkg_mult][i].GetEventWeight() ;

	  unsigned int n_pions =  particles[conf::kPdgPiP].size() + particles[conf::kPdgPiM].size() + particles[conf::kPdgPhoton].size() ; 

	  // Set rotation around q3 vector
	  kRotation->SetQVector( event_holder[bkg_mult][i].GetRecoq3() );	 

	  // remove multiplicity 2 contribution to signal...
	  // 2p0pi . 1p0pi ( multiplicity 2 . multiplicity 1 ) 
	  if( particles[conf::kPdgProton].size() == 2 && n_pions == 0 ) { 
	    double E_tot_2p[2]={0};
	    double p_perp_tot_2p[2]={0};
	    double N_prot_both = 0;
	    double P_N_2p[2]={0};

	    TVector3 V3_2prot_corr[2];
	    for ( unsigned k = 0 ; k < 2 ; ++k ) {
	      V3_2prot_corr[k] = particles[conf::kPdgProton][k].Vect() ; 
	    }
	
	    TVector3 V3_2prot_uncorr[2];
	    for ( unsigned k = 0 ; k < 2 ; ++k ) {
	      V3_2prot_uncorr[k] = particles_uncorr[conf::kPdgProton][k].Vect() ; 
	    }

	    kRotation->prot2_rot_func( V3_2prot_corr, V3_2prot_uncorr, V4_el, E_tot_2p, p_perp_tot_2p, P_N_2p , &N_prot_both);

	    for( unsigned int j = 0 ; j < bkg_mult ; ++j ) {
	      std::vector<TLorentzVector> corr_mom = { particles[conf::kPdgProton][j] } ;
	      std::vector<TLorentzVector> uncorr_mom = { particles_uncorr[conf::kPdgProton][j] } ;
	      t_particles[conf::kPdgProton] = corr_mom ;
	      t_particles_uncorr[conf::kPdgProton] = uncorr_mom ;
	  
	      // Store background event with the detected proton
	      event_holder[bkg_mult][i].SetFinalParticlesKinematics( t_particles ) ; 
	      event_holder[bkg_mult][i].SetFinalParticlesUnCorrKinematics( t_particles_uncorr ) ; 
	      event_holder[bkg_mult][i].SetEventWeight( -P_N_2p[j] * event_wgt ) ; 

	      if ( event_holder.find(min_mult) != event_holder.end() ) {
		event_holder[min_mult].push_back( event_holder[bkg_mult][i] ) ; 
	      } else {
		std::vector<T> temp = { event_holder[bkg_mult][i] } ;
		event_holder[min_mult] = temp ; 
	      }
	  
	    }
	  } // close 2p0pi topology check 
	  
	  if( particles[conf::kPdgProton].size() == 1 && n_pions == 1 ) { 
	    TVector3 V3_prot_uncorr[1] ; 
	    for ( unsigned k = 0 ; k < 1 ; ++k ) {
              V3_prot_uncorr[k] = particles_uncorr[conf::kPdgProton][k].Vect() ;
            }

	    TVector3 V3_pi_corr[1] ; 
	    double pi_charge[1] ; 
	    if ( particles[conf::kPdgPiP].size() == 1 ) {
	      V3_pi_corr[0] = particles[conf::kPdgPiP][0].Vect() ;  
	      pi_charge[0] = 1; 
	    } else if ( particles[conf::kPdgPiM].size() == 1 ) {
	      V3_pi_corr[0] = particles[conf::kPdgPiM][0].Vect() ;
	      pi_charge[0] = -1 ;
	    } else if ( particles[conf::kPdgPhoton].size() == 1 ) {
	      V3_pi_corr[0] = particles[conf::kPdgPhoton][0].Vect() ;  
	      pi_charge[0] = 0 ; 
	    }

	    double N_piphot_det;
	    double N_piphot_undet;

	    kRotation->prot1_pi1_rot_func(V3_prot_uncorr[0],V3_pi_corr[0], pi_charge[0], &N_piphot_det,&N_piphot_undet);

	    std::vector<TLorentzVector> corr_mom = { particles[conf::kPdgProton][0] } ;
	    std::vector<TLorentzVector> uncorr_mom = { particles_uncorr[conf::kPdgProton][0] } ;
	    t_particles[conf::kPdgProton] = corr_mom ;
	    t_particles_uncorr[conf::kPdgProton] = uncorr_mom ;
	  
	    // Store background event with the detected proton
	    event_holder[bkg_mult][i].SetFinalParticlesKinematics( t_particles ) ; 
	    event_holder[bkg_mult][i].SetFinalParticlesUnCorrKinematics( t_particles_uncorr ) ; 
	    event_holder[bkg_mult][i].SetEventWeight( -N_piphot_undet/N_piphot_det * event_wgt ) ; 
	    
	    if ( event_holder.find(min_mult) != event_holder.end() ) {
	      event_holder[min_mult].push_back( event_holder[bkg_mult][i] ) ; 
	    } else {
	      std::vector<T> temp = { event_holder[bkg_mult][i] } ;
	      event_holder[min_mult] = temp ; 
	    } 
	  }
	  
	  particles.clear() ;
	  particles_uncorr.clear() ;
	} // close loop over bkg m 2 events
      } /// outside of m 2 bkg

  
      bkg_mult = min_mult + 2 ;
      if( bkg_mult > max_mult ) return true ; 
      // remove multiplicity 3 contribution to signal...
      if ( event_holder.find(bkg_mult) != event_holder.end() ) {
	for ( unsigned int i = 0 ; i < event_holder[bkg_mult].size() ; ++i ) {
	  particles = event_holder[bkg_mult][i].GetFinalParticles4Mom();
	  particles_uncorr = event_holder[bkg_mult][i].GetFinalParticlesUnCorr4Mom();
	  V4_el = event_holder[bkg_mult][i].GetOutLepton4Mom();
	  kRotation->SetQVector( event_holder[bkg_mult][i].GetRecoq3() );	
	  double event_wgt = event_holder[bkg_mult][i].GetEventWeight() ;
      
	  // 2p1pi 
	  unsigned int n_pions = particles[conf::kPdgPiP].size() + particles[conf::kPdgPiM].size() + particles[conf::kPdgPhoton].size() ;  
	  if( particles[conf::kPdgProton].size() == 2 && n_pions == 1 ) {
	
	    TVector3 V3_2prot_corr[2], V3_2prot_uncorr[2];
	    for ( unsigned k = 0 ; k < 2 ; ++k ) {
	      V3_2prot_corr[k] = particles[conf::kPdgProton][k].Vect() ; 
	      V3_2prot_uncorr[k] = particles_uncorr[conf::kPdgProton][k].Vect() ; 
	    }
	
	    TVector3 V3_1pi_corr(0,0,0) ; 
	    int pi_charge = 0 ; 
	    if( particles[conf::kPdgPiP].size() == 1 ) {
	      V3_1pi_corr = particles[conf::kPdgPiP][0].Vect() ; 
	      pi_charge = 1 ; 
	    } else if( particles[conf::kPdgPiM].size() == 1 ) {
	      V3_1pi_corr = particles[conf::kPdgPiM][0].Vect() ;
	      pi_charge = -1 ; 
	    } else if( particles[conf::kPdgPhoton].size() == 1 ) V3_1pi_corr = particles[conf::kPdgPhoton][0].Vect() ; 
	
	    double Ecal_2p1pi_to2p0pi[2] = {0};
	    double p_miss_perp_2p1pi_to2p0pi[2]={0}; 
	    double P_2p1pito2p0pi[2] = {0};
	    double P_2p1pito1p1pi[2] = {0};
	    double P_2p1pito1p0pi[2] = {0};
	    double Ptot = 0;
	  
	    kRotation->prot2_pi1_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_1pi_corr, pi_charge, 
					  V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,
					  P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

	    // 2p 1pi . 2p 0pi
	    // Store one proton only so ECal calculation is consistent with master
	    for( unsigned int j = 0 ; j < particles[conf::kPdgProton].size() ; ++j ) {
	      std::vector<TLorentzVector> corr_mom = { particles[conf::kPdgProton][j] } ;
	      std::vector<TLorentzVector> uncorr_mom = { particles_uncorr[conf::kPdgProton][j] } ;
	      t_particles[conf::kPdgProton] = corr_mom ;
	      t_particles_uncorr[conf::kPdgProton] = uncorr_mom ;
	      // Store information
	      event_holder[bkg_mult][i].SetFinalParticlesKinematics( t_particles ) ; 
	      event_holder[bkg_mult][i].SetFinalParticlesUnCorrKinematics( t_particles_uncorr ) ; 
	      event_holder[bkg_mult][i].SetEventWeight( P_2p1pito2p0pi[j] * event_wgt ) ; 
	  
	      if ( event_holder.find(min_mult) != event_holder.end() ) {
		event_holder[min_mult].push_back( event_holder[bkg_mult][i] ) ; 
	      } else {
		std::vector<T> temp = { event_holder[bkg_mult][i] } ;
		event_holder[min_mult] = temp ; 
	      }
	    }

	    // 2p 1pi . 1p 1pi
	    // Store one proton only for a consistent ECal calculation
	    for( unsigned int j = 0 ; j < particles[conf::kPdgProton].size() ; ++j ) {
	      std::vector<TLorentzVector> corr_mom = { particles[conf::kPdgProton][j] } ;
	      std::vector<TLorentzVector> uncorr_mom = { particles_uncorr[conf::kPdgProton][j] } ;
	      t_particles[conf::kPdgProton] = corr_mom ;
	      t_particles_uncorr[conf::kPdgProton] = uncorr_mom ;
	      // Store information
	      event_holder[bkg_mult][i].SetFinalParticlesKinematics( t_particles ) ; 
	      event_holder[bkg_mult][i].SetFinalParticlesUnCorrKinematics( t_particles_uncorr ) ; 
	      event_holder[bkg_mult][i].SetEventWeight( P_2p1pito1p1pi[j] * event_wgt ) ; 
	      if ( event_holder.find(min_mult) != event_holder.end() ) {
		event_holder[min_mult].push_back( event_holder[bkg_mult][i] ) ; 
	      } else {
		std::vector<T> temp = { event_holder[bkg_mult][i] } ;
		event_holder[min_mult] = temp ; 
	      }
	    }
	    // 2p 1pi . 1p 0pi
	    // Store one proton only so ECal calculation is consistent with master
	    for( unsigned int j = 0 ; j < particles[conf::kPdgProton].size() ; ++j ) {
	      std::vector<TLorentzVector> corr_mom = { particles[conf::kPdgProton][j] } ;
	      std::vector<TLorentzVector> uncorr_mom = { particles_uncorr[conf::kPdgProton][j] } ;
	      t_particles[conf::kPdgProton] = corr_mom ;
	      t_particles_uncorr[conf::kPdgProton] = uncorr_mom ;
	      // Store information
	      event_holder[bkg_mult][i].SetFinalParticlesKinematics( t_particles ) ; 
	      event_holder[bkg_mult][i].SetFinalParticlesUnCorrKinematics( t_particles_uncorr ) ; 
	      event_holder[bkg_mult][i].SetEventWeight( -P_2p1pito1p0pi[j] * event_wgt ) ; 
	      if ( event_holder.find(min_mult) != event_holder.end() ) {
		event_holder[min_mult].push_back( event_holder[bkg_mult][i] ) ; 
	      } else {
		std::vector<T> temp = { event_holder[bkg_mult][i] } ;
		event_holder[min_mult] = temp ; 
	      }
	    }
	  } // outside topology 2p1pi check 
    
      
	  // 3p0pi .2p . 1p
	  if( particles[conf::kPdgProton].size() == 3 && n_pions == 0 ) {
	    const int N_3p = 3;
	    const int N_comb = 3; // N_2p!/(N_2p!*(N_3p-N_2p)!
	    const int N_2p = 2;
	    TVector3 V3_prot_uncorr[N_3p],V3_prot_corr[N_3p];
	    for ( unsigned k = 0 ; k < N_3p ; ++k ) {
	      V3_prot_corr[k] = particles[conf::kPdgProton][k].Vect() ; 
	      V3_prot_uncorr[k] = particles_uncorr[conf::kPdgProton][k].Vect() ; 
	    }

	    // Store all possible combinations
	    TLorentzVector TLV_prot_corr[N_comb][N_2p] ; 
	    TLorentzVector TLV_prot_uncorr[N_comb][N_2p] ; 
	    unsigned int count = 0 ; 
	    for( unsigned int i1 = 0; i1 < N_comb; ++i1) {
	      for( unsigned int i2 = 0; i2 < N_comb; ++i2) {
		if( ( i1 == i2 ) || ( i1 > i2 ) ) continue ; 
		TLV_prot_corr[count][0] = particles[conf::kPdgProton][i1];
		TLV_prot_corr[count][1] = particles[conf::kPdgProton][i2];
		++count ; 
	      }
	    }
  
	    double P_3pto1p[N_3p];
	    double N_p1[N_3p]={0};
	    double N_p_three=0;
	    double E_cal_3pto1p[3]={0};
	    double p_miss_perp_3pto1p[3]={0};
	    double E_cal_3pto2p[3][N_2p]={0};
	    double p_miss_perp_3pto2p[3][N_2p]={0};
	    double P_3pto2p[3][N_2p]={0};

	    kRotation->prot3_rot_func( V3_prot_corr,V3_prot_uncorr,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, E_cal_3pto1p,p_miss_perp_3pto1p,&N_p_three);

	    //3p to 2p.1p
	    for(int c = 0; c < N_comb; c++) { // Loop over number of combinations
	      for(int j = 0; j < N_2p; j++) { // Loop over two protons
		std::vector<TLorentzVector> corr_mom = { TLV_prot_corr[c][j] };
		std::vector<TLorentzVector> uncorr_mom = {  TLV_prot_uncorr[c][j] };
		t_particles[conf::kPdgProton] = corr_mom ;
		t_particles_uncorr[conf::kPdgProton] = uncorr_mom; 
		// Store information
		event_holder[bkg_mult][i].SetFinalParticlesKinematics( t_particles ) ; 
		event_holder[bkg_mult][i].SetFinalParticlesUnCorrKinematics( t_particles_uncorr ) ; 
		event_holder[bkg_mult][i].SetEventWeight( P_3pto2p[c][j] * event_wgt ) ; 
		if ( event_holder.find(min_mult) != event_holder.end() ) {
		  event_holder[min_mult].push_back( event_holder[bkg_mult][i] ) ; 
		} else {
		  std::vector<T> temp = { event_holder[bkg_mult][i] } ;
		  event_holder[min_mult] = temp ; 
		}
	      }
	    }
	    // 3p . 1p 
	    for(int j = 0; j < N_3p; j++) {
	      P_3pto1p[j]= N_p1[j] / N_p_three;

	      std::vector<TLorentzVector> corr_pmom = { particles[conf::kPdgProton][j] } ; 
	      std::vector<TLorentzVector> uncorr_pmom = {particles_uncorr[conf::kPdgProton][j] }; 
	      t_particles[conf::kPdgProton] = corr_pmom ;
	      t_particles_uncorr[conf::kPdgProton] = uncorr_pmom ;
	      // Store information
	      event_holder[bkg_mult][i].SetFinalParticlesKinematics( t_particles ) ; 
	      event_holder[bkg_mult][i].SetFinalParticlesUnCorrKinematics( t_particles_uncorr ) ; 
	      event_holder[bkg_mult][i].SetEventWeight( -P_3pto1p[j] * event_wgt ) ; 
	      if ( event_holder.find(min_mult) != event_holder.end() ) {
		event_holder[min_mult].push_back( event_holder[bkg_mult][i] ) ; 
	      } else {
		std::vector<T> temp = { event_holder[bkg_mult][i] } ;
		event_holder[min_mult] = temp ; 
	      }
	    }
	  } // outside topology 3p 
	} //
      }
  
      bkg_mult = min_mult + 3 ;
      if( bkg_mult > max_mult ) return true ; 
      // remove multiplicity 4 contribution to signal...
      if ( event_holder.find(bkg_mult) != event_holder.end() ) {
	for ( unsigned int i = 0 ; i < event_holder[bkg_mult].size() ; ++i ) {
	  particles = event_holder[bkg_mult][i].GetFinalParticles4Mom();
	  particles_uncorr = event_holder[bkg_mult][i].GetFinalParticlesUnCorr4Mom();
	  V4_el = event_holder[bkg_mult][i].GetOutLepton4Mom();
	  kRotation->SetQVector( event_holder[bkg_mult][i].GetRecoq3() );	
	  double event_wgt = event_holder[bkg_mult][i].GetEventWeight() ;
	  unsigned int n_pions = particles[conf::kPdgPiP].size() + particles[conf::kPdgPiM].size() + particles[conf::kPdgPhoton].size() ;  

	  // 2p2pi 
	  if( particles[conf::kPdgProton].size() == 2 && n_pions == 2 ) { 
	    const int N_2pi=2;
	    double Ecal_2p2pi[2];
	    double p_miss_perp_2p2pi[2];
	    double Ptot_2p[2]={0};
	    double pion_acc_ratio[N_2pi] = {1};

	    TVector3 V3_2prot_corr[2];
	    TVector3 V3_2prot_uncorr[2];
	    TVector3 V3_2pi_corr[N_2pi];

	    for ( unsigned k = 0 ; k < 2 ; ++k ) {
	      V3_2prot_corr[k] = particles[conf::kPdgProton][k].Vect() ; 
	      V3_2prot_uncorr[k] = particles_uncorr[conf::kPdgProton][k].Vect() ; 
	    }
	    /*
	      for ( unsigned k = 0 ; k < 2 ; ++k ) {
	      V3_2pi_corr[k] = particles[conf::kPdgPiP][k].Vect() ; 
	      }
	      kRotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi_corr,charge_pi, V4_el, Ecal_2p2pi,p_miss_perp_2p2pi,Ptot_2p);
	    */
	    // THIS IS NOT COMPLETE - MUST ADD ALL PIONS
	  }
	  // 3p1pi
	  if( particles[conf::kPdgProton].size() == 3 && n_pions == 1 ) {
	    double P_tot_3p[3]={0};
	    double Ecal_3p1pi[3]={0};
	    double p_miss_perp_3p1pi[3]={0};
	    TVector3 V3_prot_uncorr[3],V3_prot_corr[3];
	    for ( unsigned k = 0 ; k < 3 ; ++k ) {
	      V3_prot_corr[k] = particles[conf::kPdgProton][k].Vect() ; 
	      V3_prot_uncorr[k] = particles_uncorr[conf::kPdgProton][k].Vect() ; 
	    }

	    TVector3 V3_pi_corr(0,0,0) ; 
	    int pi_charge = 0 ; 
	    if( particles[conf::kPdgPiP].size() == 1 ) {
	      V3_pi_corr = particles[conf::kPdgPiP][0].Vect() ; 
	      pi_charge = 1 ; 
	    } else if( particles[conf::kPdgPiM].size() == 1 ) {
	      V3_pi_corr = particles[conf::kPdgPiM][0].Vect() ;
	      pi_charge = -1 ; 
	    } else if( particles[conf::kPdgPhoton].size() == 1 ) V3_pi_corr = particles[conf::kPdgPhoton][0].Vect() ; 

	    kRotation->prot3_pi1_rot_func(V3_prot_corr, V3_prot_uncorr, V3_pi_corr, pi_charge , V4_el, Ecal_3p1pi, p_miss_perp_3p1pi, P_tot_3p);

	    for(int j = 0; j < 3; j++) {
	      std::vector<TLorentzVector> corr_pmom = { particles[conf::kPdgProton][j] } ; 
	      std::vector<TLorentzVector> uncorr_pmom = {particles_uncorr[conf::kPdgProton][j] }; 
	      t_particles[conf::kPdgProton] = corr_pmom ;
	      t_particles_uncorr[conf::kPdgProton] = uncorr_pmom ;
	      // Store information
	      event_holder[bkg_mult][i].SetFinalParticlesKinematics( t_particles ) ; 
	      event_holder[bkg_mult][i].SetFinalParticlesUnCorrKinematics( t_particles_uncorr ) ; 
	      event_holder[bkg_mult][i].SetEventWeight( P_tot_3p[j] * event_wgt ) ; 
	      if ( event_holder.find(min_mult) != event_holder.end() ) {
		event_holder[min_mult].push_back( event_holder[bkg_mult][i] ) ; 
	      } else {
		std::vector<T> temp = { event_holder[bkg_mult][i] } ;
		event_holder[min_mult] = temp ; 
	      }
	    }
	  }
	}
      }
      return true ; 
    }
      
  protected:
    virtual ~BackgroundI();
    Subtraction * kRotation = nullptr ;

  };
}

#endif
