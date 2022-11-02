/** 
 * \info This script contains general information on particle properties
 **/

#ifndef _PARTICLE_I_H_
#define _PARTICLE_I_H_

namespace e4nu {
  namespace conf { 
   
    // Pdg codes 
    const int kPdgProton = 2212 ;
    const int kPdgNeutron = 2112 ;
    const int kPdgPiP = 211 ; 
    const int kPdgPiM = -211 ; 
    const int kPdgPi0 = 111 ; 
    const int kPdgElectron = 11 ; 
    const int kPdgPositron = -11 ; 

    // Mass
    const double kMassProton = 0.9382720813 ; 
    const double kMassNeutron = 0.939565 ; 
    const double kMassPiM = 0.13957 ; 
    const double kMassPim = 0.139570 ; 
    const double kMassPi0 = 0.139570 ;
    const double kMassElectron = 0.000510998 ; 

    // Detector resolution for each particle 
    const double kProtonRes = 0.01 ; 
    const double kElectronRes = 0.005 ; 
    const double kPionRes = 0.007 ; 
  }
}
#endif 
