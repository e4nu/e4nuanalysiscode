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
#include "utils/Subtraction.h"
#include "conf/ParticleI.h"

namespace e4nu {

  class BackgroundI : public ConfigureI {

  public:
    // Default constructor
    BackgroundI() ;
    BackgroundI( const std::string input_file ) ;
    BackgroundI( const double EBeam, const unsigned int TargetPdg ) ;
    void Initialize(void);

    // Some useful functions:
    unsigned int GetMinParticleMultiplicity( int pdg ) const ;

    // Definition of Background mehtod
    // These template class guarantees the same substraction method for data and MC
    // Definition must be in header file to avoid linking issuess
    template <class T>
    bool BackgroundSubstraction( std::map<int,std::vector<T>> & event_holder ) {
      bool apply_fiducial = ApplyFiducial();

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

            // We want to correct for the acceptance so we have a "perfectly detected" event independent of the acceptance
            ApplyAcceptanceCorrection( event_holder[m][event_id], true ) ;

            // Start rotations
            for ( unsigned int rot_id = 0 ; rot_id < GetNRotations() ; ++rot_id ) {
              // Set rotation around q3 vector
              TVector3 VectorRecoQ = event_holder[m][event_id].GetRecoq3() ;
              double rotation_angle = gRandom->Uniform(0,2*TMath::Pi());

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
                  bool is_particle_contained = fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_vect, IsData(), apply_fiducial ) ;

                  TLorentzVector part_shift = (it->second)[part_id] ;
                  if( fFidAngleShift != 0 ) {
                    part_shift.SetPhi( (it->second)[part_id].Phi() + fFidAngleShift * TMath::Pi() / 180. ) ;
                    is_particle_contained *= fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_shift.Vect(), IsData(), apply_fiducial ) ;
                    part_shift.SetPhi( (it->second)[part_id].Phi() - fFidAngleShift * TMath::Pi() / 180. ) ;
                    is_particle_contained *= fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_shift.Vect(), IsData(), apply_fiducial ) ;
                  }

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
            for( auto it = probability_count.begin() ; it != probability_count.end() ; ++it ) {
              for( auto it_key = it->first.begin() ; it_key != it->first.end() ; ++it_key ) {
                T temp_event = event_holder[m][event_id] ;
                double event_wgt = temp_event.GetEventWeight() ;
                std::map<int,std::vector<TLorentzVector>> particles = temp_event.GetFinalParticles4Mom() ;
                std::map<int,std::vector<TLorentzVector>> particles_uncorr = temp_event.GetFinalParticlesUnCorr4Mom() ;

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

                temp_event.SetFinalParticlesKinematics( temp_corr_mom ) ;
                temp_event.SetFinalParticlesUnCorrKinematics( temp_uncorr_mom ) ;

                double probability = - (it->second) * event_wgt / N_all ;
                temp_event.SetEventWeight( probability ) ;

                // Store analysis record after background substraction (4) :
                temp_event.StoreAnalysisRecord(kid_bkgcorr+m); // Id is the bkg id (4) + original multiplicity.
                // For m = signal_multiplicity, id = 3+signal_mult

                // Account for acceptance
                ApplyAcceptanceCorrection( temp_event ) ;

                // Add in map
                if ( event_holder.find(new_multiplicity) != event_holder.end() ) {
                  event_holder[new_multiplicity].push_back( temp_event ) ;
                } else {
                  std::vector<T> temp = { temp_event } ;
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
    bool HadronsAcceptanceCorrection( std::map<int,std::vector<T>> & event_holder ) {

      // We need to correct for signal events that are reconstructed outside of the fiducial
      bool apply_fiducial = ApplyFiducial();
      if( !GetSubtractBkg() || ! apply_fiducial ) return true ;

      Fiducial * fiducial = GetFiducialCut() ;
      if( !fiducial ) return false ;
      std::cout << " Applying Acceptance Correction to hadrons ... " << std::endl;

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
          TVector3 VectorRecoQ = signal_events[i].GetRecoq3() ;
          double rotation_angle = gRandom->Uniform(0,2*TMath::Pi());

          std::map<int,std::vector<TLorentzVector>> rot_particles = signal_events[i].GetFinalParticles4Mom() ;
          std::map<int,std::vector<TLorentzVector>> rot_particles_uncorr = signal_events[i].GetFinalParticlesUnCorr4Mom() ;

          // Rotate all particles
          bool is_contained = true ;
          for( auto it = rot_particles.begin() ; it != rot_particles.end() ; ++it ) {
            int part_pdg = it->first ;
            if( Topology.find( part_pdg ) == Topology.end() ) continue ; // Skip particles which are not in signal definition

            for ( unsigned int part_id = 0 ; part_id < (it->second).size() ; ++part_id ) {
              TVector3 part_vect = (it->second)[part_id].Vect() ;
              part_vect.Rotate(rotation_angle,VectorRecoQ);

              // Check which particles are in fiducial
              is_contained = fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_vect, IsData(), apply_fiducial ) ;

              TLorentzVector part_shift = (it->second)[part_id] ;
              if( fFidAngleShift != 0 ) {
                part_shift.SetPhi( (it->second)[part_id].Phi() + fFidAngleShift * TMath::Pi() / 180. ) ;
                is_contained *= fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_shift.Vect(), IsData(), apply_fiducial ) ;
                part_shift.SetPhi( (it->second)[part_id].Phi() - fFidAngleShift * TMath::Pi() / 180. ) ;
                is_contained *= fiducial->FiducialCut( part_pdg, GetConfiguredEBeam(), part_shift.Vect(), IsData(), apply_fiducial ) ;
              }

              if( !is_contained ) break ;
            }
            if( !is_contained ) break ;
          }
          if( is_contained ) ++N_signal_detected ;
          else ++N_signal_undetected ;
        }
        if( N_signal_detected == 0 ) continue ;
        // Add missing signal events
        T temp_event = signal_events[i] ;
        double event_wgt = temp_event.GetEventWeight() ;
        temp_event.SetEventWeight( + event_wgt * N_signal_undetected / N_signal_detected ) ;
        signal_events.push_back(temp_event);

        // Store analysis record after acceptance correction (3) :
        temp_event.StoreAnalysisRecord(kid_acc);

      }
      // Store correction
      event_holder[min_mult] = signal_events ;

      return true ;
    }


    template <class T>
    bool ElectronAcceptanceCorrection( std::map<int,std::vector<T>> & event_holder ) {

      // We need to correct for signal events that are reconstructed outside of the fiducial
      bool apply_fiducial = ApplyFiducial();
      if( !GetSubtractBkg() || apply_fiducial ) return true ;

      Fiducial * fiducial = GetFiducialCut() ;
      if( !fiducial ) return false ;
      std::cout << " Applying Electron Acceptance Correction ... " << std::endl;

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
          TVector3 BeamVector (0,0,1);
          double rotation_angle = gRandom->Uniform(0,2*TMath::Pi());

          TVector3 emom = signal_events[i].GetOutLepton4Mom().Vect() ;

          // Rotate all particles
          emom.Rotate(rotation_angle, BeamVector);

          bool is_contained = fiducial->FiducialCut( conf::kPdgElectron, GetConfiguredEBeam(), emom, IsData(), apply_fiducial ) ;

          TLorentzVector part_shift = signal_events[i].GetOutLepton4Mom();
          if( fFidAngleShift != 0 ) {
            part_shift.SetPhi( emom.Phi() + fFidAngleShift * TMath::Pi() / 180. ) ;
            is_contained *= fiducial->FiducialCut( conf::kPdgElectron, GetConfiguredEBeam(), part_shift.Vect(), IsData(), apply_fiducial ) ;
            part_shift.SetPhi( emom.Phi() - fFidAngleShift * TMath::Pi() / 180. ) ;
            is_contained *= fiducial->FiducialCut( conf::kPdgElectron, GetConfiguredEBeam(), part_shift.Vect(), IsData(), apply_fiducial ) ;
          }

          if( is_contained ) ++N_signal_detected ;
          else ++N_signal_undetected ;
        }
        if( N_signal_detected == 0 ) continue ;
        // Add missing signal events
        T temp_event = signal_events[i] ;

        double event_wgt = temp_event.GetEventWeight() ;
        temp_event.SetEventWeight( + event_wgt * N_signal_undetected / N_signal_detected ) ;
        signal_events.push_back(temp_event);

        // Store analysis record after acceptance correction (3) :
        temp_event.StoreAnalysisRecord(kid_acc);
      }
      // Store correction
      event_holder[min_mult] = signal_events ;

      return true ;
    }

  protected:
    virtual ~BackgroundI();
    Subtraction * kRotation = nullptr ;

  };
}

#endif
