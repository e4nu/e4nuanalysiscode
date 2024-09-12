#ifndef _ACCEPTANCE_UTILS_H_
#define _ACCEPTANCE_UTILS_H_

#include "plotting/PlottingUtils.h"

namespace e4nu {
  namespace plotting {

    // Observable definition;
    double TotWeight;
    double ECal, Recoq3, RecoW;
    double pfl, pfl_theta, pfl_phi;
    double proton_mom, proton_phi, proton_theta;
    double pim_mom, pim_theta, pim_phi;
    double pip_mom, pip_theta, pip_phi;
    double HadAlphaT, HadDeltaPT, HadDeltaPTx, HadDeltaPTy, HadDeltaPhiT;
    double AlphaT, DeltaPT, DeltaPhiT;
    double RecoXBJK, RecoEnergyTransfer, RecoQ2, HadSystemMass, RecoQELEnu;
    double MissingEnergy, MissingAngle, MissingMomentum;
    double InferedNucleonMom;
    double HadronsAngle, Angleqvshad;
    double AdlerAngleThetaP, AdlerAnglePhiP, AdlerAngleThetaPi, AdlerAnglePhiPi;
    double RecoEvPion, RecoWPion;
    double ElectronPT, PionPT;
    long NEntries;

    double GetContent( const std::string observable ) ;

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



