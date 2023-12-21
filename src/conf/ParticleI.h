/** 
 * \info This script contains general information on particle properties
 **/

#ifndef _PARTICLE_I_H_
#define _PARTICLE_I_H_

namespace e4nu {
  namespace conf { 
   
    // Pdg codes 
    static const int kPdgProton = 2212 ;
    static const int kPdgNeutron = 2112 ;
    static const int kPdgPiP = 211 ; 
    static const int kPdgPiM = -211 ; 
    static const int kPdgPi0 = 111 ; 
    static const int kPdgKP = 321 ;
    static const int kPdgKM = -321 ;
    static const int kPdgK0 = 311 ;
    static const int kPdgElectron = 11 ; 
    static const int kPdgPositron = -11 ; 
    static const int kPdgPhoton = 22 ; 
    // Mass
    static const double kProtonMass = 0.9382720813 ; 
    static const double kNeutronMass = 0.939565 ; 
    static const double kPiPMass = 0.13957 ; 
    static const double kPiMMass = 0.139570 ; 
    static const double kPi0Mass = 0.139570 ;
    static const double kElectronMass = 0.000510998 ; 
    static const double kNucleonMass    =  (kProtonMass+kNeutronMass)/2.;  //GeV

    // Detector resolution for each particle 
    static const double kProtonRes = 0.01 ; 
    static const double kElectronRes = 0.005 ; 
    static const double kPionRes = 0.007 ; 

    // Particle charge
    static const int kElectronCharge = -1 ; 
    static const int kProtonCharge = 1 ;
    static const int kPiPCharge = 1 ;
    static const int kPiMCharge = - 1 ; 
    static const int kPi0Charge = 0 ;
    static const int kNeutronCharge = 0 ; 
    static const int kPhotonCharge = 0 ; 
  }
}
#endif 
