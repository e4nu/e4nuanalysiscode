#ifndef _SYST_UTILS_H_
#define _SYST_UTILS_H_

#include "plotting/PlottingUtils.h"

namespace e4nu {

  namespace plotting {
    // Observables defined in plotting utils
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
  }

  namespace systematics {
    // This function calculates the corresponding xsec distributions for a given observable
    // No acceptance correction is used
    void ComputeHistSyst( std::vector<std::string> input_files, std::vector<std::string> tags, std::string observable, bool is_data,
			                    std::string input_location, std::string output_location, std::string analysis_id );
    void AddSystematic( TH1D & hist, const double rel_error, const std::string name ) ;
    TH1D * AddSystematic( TH1D & hist, const TH1D & hist_w_error ) ;
    TH1D * SectorVariationError( TH1D & hist, const std::vector<TH1D*> h_per_sector ) ;
  }
}

#endif
