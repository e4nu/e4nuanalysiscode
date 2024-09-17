#ifndef _ACCEPTANCE_UTILS_H_
#define _ACCEPTANCE_UTILS_H_

#include "plotting/PlottingUtils.h"

namespace e4nu {
  namespace plotting {

    // Define observables here
    extern double TotWeight, ECal, Recoq3, RecoW;
    extern double pfl, pfl_theta, pfl_phi;
    extern double proton_mom, proton_phi, proton_theta;
    extern double pim_mom, pim_theta, pim_phi;
    extern double pip_mom, pip_theta, pip_phi;
    extern double HadAlphaT, HadDeltaPT, HadDeltaPTx, HadDeltaPTy, HadDeltaPhiT;
    extern double AlphaT, DeltaPT, DeltaPhiT;
    extern double RecoXBJK, RecoEnergyTransfer, RecoQ2, HadSystemMass, RecoQELEnu;
    extern double MissingEnergy, MissingAngle, MissingMomentum;
    extern double InferedNucleonMom, HadronsAngle, Angleqvshad;
    extern double AdlerAngleThetaP, AdlerAnglePhiP, AdlerAngleThetaPi, AdlerAnglePhiPi;
    extern double RecoEvPion, RecoWPion, ElectronPT, PionPT;

    // This function computes the acceptance for a given observable
    // It uses the mc file location to find the _true and _truereco.root files
    std::string ComputeAcceptance(std::vector<std::string> mc_files, std::string observable, std::string title,
				  std::string input_MC_location, std::string output_location, std::string output_file_name,
				  std::string analysis_id = "default", bool store_root=false) ; 


    std::string ComputeRadCorr(std::vector<std::string> mc_files, std::string observable, std::string title,
			       std::string input_MC_location, std::string output_location,  std::string output_file_name, 
			       std::string analysis_id="default", bool store_root=false) ;
      
  }
}

#endif



