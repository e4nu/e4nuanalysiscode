/**
 * This file contains utils specific for particles 
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _PARTICLE_UTILS_H_
#define _PARTICLE_UTILS_H_

namespace e4nu { 
  namespace ParticleUtils
    {
      unsigned double GetParticleResolucion( int particle_pdg, bool apply_resolution = true ) ; 
     }
}

#endif _PARTICLE_UTILS_H_
