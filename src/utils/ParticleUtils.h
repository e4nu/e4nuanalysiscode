/**
 * This file contains utils specific for particles 
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _PARTICLE_UTILS_H_
#define _PARTICLE_UTILS_H_
#include <iostream>
#include <string> 
#include "TLorentzVector.h"

namespace e4nu { 
  namespace utils
    {
      void ApplyResolution( const int pdg, TLorentzVector & mom, const double EBeam ) ; 
      double GetParticleResolucion( const int particle_pdg, const double EBeam ) ; 
      double GetParticleMass( const int pdg ) ; 
      int GetParticleCharge( const int pdg ) ;
      std::string PdgToString( const int pdg ) ;
    }
}

#endif 
