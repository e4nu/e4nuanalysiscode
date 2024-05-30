/**
 * This class deals with common analysis features which are present in data and MC
 **/

#ifndef _ANALYSIS_I_H_
#define _ANALYSIS_I_H_

#include <vector>
#include <map>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "utils/Fiducial.h"
#include "analysis/BackgroundI.h"

namespace e4nu {

  class AnalysisI : public BackgroundI {

  public:
    // Default constructor
    AnalysisI() ;
    AnalysisI( const std::string input_file ) ;
    AnalysisI( const double EBeam, const unsigned int TargetPdg ) ;

    bool Analyse( Event & event ) ;
    void Initialize(void) ;
    bool Finalise(void) ;

  protected :
    void CookEvent( Event & event ) ;
    void PlotBkgInformation( const Event event ) ;
    void ApplyMomentumCut( Event & event ) ;
    bool ApplyFiducialCut( Event & event, bool apply_fiducial ) ; 
    double GetElectronMinTheta( TLorentzVector emom ) ;

    // Loading all observables
    void LoadFinalObservables( Event event );
    double GetObservable( std::string obs ); 

    // ID for Background historams
    unsigned int kid_signal, kid_tottruebkg, kid_totestbkg, kid_acccorr;
    unsigned int kid_2p0pitruebkg, kid_1p1pitruebkg, kid_2p1pitruebkg, kid_1p2pitruebkg ;
    unsigned int kid_2p0piestbkg, kid_1p1piestbkg, kid_2p1piestbkg, kid_1p2piestbkg ;

    virtual ~AnalysisI();

  private :
    TF1 * kElectronFit = nullptr ;

    // XSec value. if data it is set to 0
    double kXSec = 0 ;


  // Definition of final observables
  bool fCC, fNC, fEM, fQEL, fRES, fMEC, fDIS, fIsBkg ;
  int fID, fTargetPdg, fInLeptonPdg, fOutLeptonPdg, kresid ;
  unsigned int fTrueNProtons, fTrueNNeutrons, fTrueNPiP, fTrueNPiM, fTrueNPi0, fTrueNKP, fTrueNKM, fTrueNK0, fTrueNEM, fTrueNOther, fTopMult, fInitialNEvents; 
  unsigned int fElectronSector, fRecoNProtons, fRecoNNeutrons, fRecoNPiP, fRecoNPiM, fRecoNPi0, fRecoNKP, fRecoNKM, fRecoNK0, fRecoNEM;
  double fTotWeight, fAccWght, fEventWght, fBeamE, fTrueQ2s, fTrueWs,fTruexs, fTrueys,fTrueQ2,fTrueW,fTruex,fTruey ;
  double fEfl, fpfl, fpflx, fpfly, fpflz, fpfl_theta, fpfl_phi, fRecoQELEnu, fRecoEnergyTransfer, fRecoq3, fRecoQ2, fRecoXBJK, fRecoW, fMottXSecScale;
  double kHadronsAngle,kproton_E, kproton_mom, kproton_momx, kproton_momy, kproton_momz, kproton_theta, kproton_phi, kECal, kDiffECal, kAlphaT, kDeltaPT;
  double kDeltaPhiT, kHadAlphaT, kHadDeltaPT, kHadDeltaPTx, kHadDeltaPTy, kHadDeltaPhiT, kInferedNucleonMom ;
  double fpip_E, fpip_mom, fpip_momx, fpip_momy, fpip_momz, fpip_theta, fpip_phi;
  double fAdlerAngleThetaP, fAdlerAnglePhiP, fAdlerAngleThetaPi, fAdlerAnglePhiPi, fAngleqvshad, fHadSystemMass ;
  double fpim_E, fpim_mom, fpim_momx, fpim_momy, fpim_momz, fpim_theta, fpim_phi;
  double fMissingEnergy, fMissingMomentum, fMissingAngle, fConversionFactor, fTotalXSec, fMCNormalization;

  };
}

#endif
