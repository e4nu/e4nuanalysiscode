#ifndef _PLOTTING_UTILS_H_
#define _PLOTTING_UTILS_H_

#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"
#include <iomanip>
#include <map>

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

    double GetObservable( const std::string observable ) ;
    void NormalizeHist( TH1D * h, double normalization_factor );
    void CorrectData(TH1D* h, TH1D* acc);
    std::string GetAxisLabel( std::string observable, unsigned int id_axis );
    std::vector<double> GetUniformBinning( unsigned int nbins, double min, double max);
    std::vector<double> GetECalBinning( unsigned int nbins_tail, unsigned int nbins_peak, double min, double max, double EBeam);
    std::vector<double> GetBinning( std::string observable, double EBeam, std::string analysis_key="default" );
    std::vector<double> GetAdditionalBinning( std::string second_observable, double EBeam, std::string analysis_id="default" ) ;
    std::string GetAlternativeObs( std::string observable );
    std::string GetObsName( std::string observable );
    std::string GetUnit( std::string observable );
    double GetMaximum( std::vector<TH1D*> predictions);
    bool PlotZoomIn(std::string analysis_id="default");
    void StandardFormat( TH1D * prediction, std::string title, int color, int style, std::string observable, double y_max = 0, std::string y_axis_label ="");
    std::vector<std::string> SplitString(std::string s, char d=',' ) ;
    std::string GetArg(std::string op, int argc, char ** argv );
    bool ExistArg(std::string op, int argc, char ** argv );
  }
}

#endif
