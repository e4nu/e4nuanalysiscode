/**
 * This file contains the particle threshold maps for every particle 
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _PARTICLE_THRESHOLDS_H_
#define _PARTICLE_THRESHOLDS_H_

namespace e4nu { 
  namespace ParticleConfigurables 
    {
      unsigned double GetParticleResolucion( int particle_pdg, bool apply_resolution = true ) ; 
      unsigned double GetECalOffset( const unsigned int target_pdg ) ;
      unsigned double GetBindingEnergy( const unsigned int target_pdg ) ; 
      unsigned double GetTargetMass( const unsigned int target_pdg ) ; 
      unsigned double GetResidualTargetMass( const unsigned int target_pdg ) ; 
      unsigned int GetTargetNProtons( const unsigned int target_pdg ); 
      unsigned int GetTargetNProtons( const unsigned int target_pdg );
    }
}

#endif _PARTICLE_THRESHOLDS_H_
