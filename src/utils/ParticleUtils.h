/**
 * This file contains utils specific for particles 
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _PARTICLE_UTILS_H_
#define _PARTICLE_UTILS_H_

namespace e4nu { 
  namespace utils
    {
      double GetParticleResolucion( const int particle_pdg, const double EBeam, const bool apply_resolution = true ) ; 
     }
}

#endif 
