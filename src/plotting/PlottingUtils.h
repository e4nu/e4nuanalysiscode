#ifndef _PLOTTING_UTILS_H_
#define _PLOTTING_UTILS_H_

#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TTree.h"
#include <iomanip>
#include <map>
#include "TGraph2D.h"

namespace e4nu {
  namespace plotting {

    // Define observables here
    extern double TotWeight, ECal, BeamE, Recoq3, RecoW;
    extern double EventWght, AccWght, MottXSecScale;
    extern double Efl, pfl, pfl_theta, pfl_phi, pfl_T ;
    extern double proton_mom, proton_phi, proton_theta;
    extern double pim_mom, pim_theta, pim_phi;
    extern double pip_mom, pip_theta, pip_phi;
    extern double HadAlphaT, HadDeltaPT, HadDeltaPTx, HadDeltaPTy, HadDeltaPhiT;
    extern double AlphaT, DeltaPT, DeltaPhiT;
    extern double RecoXBJK, RecoEnergyTransfer, RecoQ2, HadSystemMass, RecoQELEnu;
    extern double MissingEnergy, MissingAngle, MissingMomentum, MissingTransMomentum, CorrMissingEnergy, CorrMissingEnergy1, CorrMissingEnergy2, CorrMissingEnergy3;
    extern double InferedNucleonMom, HadronsAngle, Angleqvshad;
    extern double AdlerAngleThetaP, AdlerAnglePhiP, AdlerAngleThetaPi, AdlerAnglePhiPi;
    extern double RecoEvPion, RecoWPion, ElectronPT, PionPT;
    extern bool IsBkg;
    extern int ElectronSector, resid, InitialNEvents;
    extern bool QEL, RES, DIS, MEC;
    extern double MCNormalization, DataNormalization;
    extern long NEntries ;
    extern int TrueNProtons, TrueNNeutrons, TrueNPiP, TrueNPiM, TrueNPi0, TrueNCh;
    // Declaring an external variable to use in the code to store the Graph relating El', Ehad with Ebeam-Emiss.
    // This is done also for slices on pt (hence _1, _2, and _3 correspond to different pt slices, hardcoded).
    extern TGraph2D* graph_oscillations, *graph_oscillations_1, *graph_oscillations_2, *graph_oscillations_3 ;
    void SetAnalysisBranch( TTree * tree ) ;

    int ColorBlindPalette(int color_id ) ;
    double GetObservable( const std::string observable ) ;
    void NormalizeHist( TH1D * h, double normalization_factor );
    void NormalizeHist( TH2D * h, double normalization_factor );
    void CorrectData(TH1D* h, TH1D* acc);
    void CorrectData(TH2D* h, TH2D* acc);
    std::string GetAxisLabel( std::string observable, unsigned int id_axis );
    std::string GetAxisLabel(std::string observable_x, std::string observable_y, unsigned int id_axis);
    std::vector<double> GetUniformBinning( unsigned int nbins, double min, double max);
    std::vector<double> GetECalBinning( unsigned int nbins_tail, unsigned int nbins_peak, double min, double max, double EBeam);
    std::vector<double> GetBinning( std::string observable, double EBeam, std::string analysis_key="default" );
    std::vector<double> GetAdditionalBinning( std::string second_observable, double EBeam, std::string analysis_id="default" ) ;
    std::string GetAlternativeObs( std::string observable );
    std::string GetObsName( std::string observable );
    std::string GetUnit( std::string observable );
    double GetMaximum( std::vector<TH1D*> predictions );
    double GetMaximum( std::vector<TH2D*> predictions );
    double GetMinimum( std::vector<TH1D*> predictions );
    std::vector<double> GetEThetaRange( TTree & tree ) ;
    double GetEPhiRange( TTree & tree ) ;
    bool PlotZoomIn(std::string analysis_id="default");
    void StandardFormat( TH1D * prediction, std::string title, int color, int style, std::string observable, bool is_log = false, double y_max = 0, std::string y_axis_label ="" );
    void StandardFormat( TH2D * prediction, std::string title, int color, int style, std::string observable_x, std::string observable_y, double z_max = 0, std::string z_axis_label ="");
    std::vector<std::string> SplitString(std::string s, char d=',' ) ;
    std::string GetArg(std::string op, int argc, char ** argv );
    bool ExistArg(std::string op, int argc, char ** argv );
    // Functions for oscillation study
    void GetMissingEnergyGraph( const std::string mc_file, const bool diff_Ebeam = false );
    double ComputeMissingEnergy( const double event_efl, const double event_ehad, const unsigned int slice = 0 );
    int GetClosestBin( TH2D * hist, double cut_value, std::string axis );
  }
}

#endif
