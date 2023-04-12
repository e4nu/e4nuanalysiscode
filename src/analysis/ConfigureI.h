/**
 * \info 
 * The analysis parameters are configurable 
 * the default values are set here
 * These can also be configured provided a configuration file (input_file)
 * The class also deals with the setup of Fiducials, ElectronFits, and other 
 * cofnigurables which are needed for both data and MC analyses
 **/

#ifndef _CONFIGURE_I_H_
#define _CONFIGURE_I_H_

#include <vector>
#include <map>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "conf/FiducialCutI.h"
#include "utils/Fiducial.h"
#include <TRandom3.h>

namespace e4nu { 

  class ConfigureI {

  public: 
    // Default constructor
    ConfigureI() ;

    // Configure with input file 
    ConfigureI( const std::string input_file ) ;
    ConfigureI( const double EBeam, const unsigned int TargetPdg ) ;
      
    void PrintConfiguration(void) const ; 
    double GetConfiguredEBeam(void) const { return kEBeam ; }
    unsigned int GetConfiguredTarget(void) const { return kTargetPdg ; }
    bool UseAllSectors(void) const { return kUseAllSectors ; } 
    bool ApplyFiducial(void) const { return kApplyFiducial ; }
    bool ApplyEFiducial(void) const { return kApplyEFiducial ; }
    bool ApplyHadFiducial(void) const { return kApplyHadFiducial ; }
    bool ApplyCorrWeights(void) const { return kApplyCorrWeights ; }
    bool ApplyAccWeights(void) const { return kApplyAccWeights ; }
    bool ApplyMottScaling(void) const { return kApplyMottWeight ; }
    bool ApplyReso(void) const { return kApplyReso ; }
    bool ApplyPhiOpeningAngle(void) const { return kApplyPhiOpeningAngle ; }
    bool UsePhiThetaBand(void) const { return kUsePhiThetaBand ; } 
    bool ApplyThetaSlice(void) const { return kApplyThetaSlice ; } 
    bool ApplyGoodSectorPhiSlice(void) const { return kApplyGoodSectorPhiSlice ; }
    bool IsData(void) const { return kIsData ; } 
    bool GetOffSet(void) const { return koffset ; } 
    bool ApplyQ2Cut(void) const{ return kQ2Cut ; }
    bool ApplyWCut(void) const{ return kWCut ; }
    bool ApplyMomCut(void) const { return kApplyMomCut ; } 
    bool ApplyOutElectronCut(void) const { return kOutEMomCut ; }      
    bool IsElectronData(void) const { return kIsElectron ; }
    bool IsConfigured(void) const { return kIsConfigured ; }
    bool IsNoFSI(void) const { return kNoFSI ; }

    unsigned int GetNEventsToRun(void) const { return kNEvents ; } 
    unsigned int GetFirstEventToRun(void) const { return kFirstEvent ; } 
    unsigned int GetMaxBkgMult(void) const { return kMaxBkgMult ; }
    unsigned int GetMinBkgMult(void) const { return kMult_signal ; }
    unsigned int GetNRotations(void) const { return kNRotations ; } 
    bool GetSubtractBkg(void) const { return kSubtractBkg ; }
    bool GetDebugBkg(void) const { return kDebugBkg ; } 

    std::map<int,unsigned int> GetTopology(void) const{ return kTopology_map ; } 
    unsigned int GetNTopologyParticles(void) ;    
    Fiducial * GetFiducialCut(void) { return kFiducialCut ; } 

    double GetElectronMinTheta( TLorentzVector emom ) ;      
    
    // Histogram Configurables
    std::vector<std::string> GetObservablesTag(void) const { return kObservables ; }
    std::vector<unsigned int> GetNBins(void) const { return kNBins ; }
    std::vector<std::vector<double>> GetRange(void) const { return kRanges ; } 
    bool NormalizeHist(void) { return kNormalize ; }

    std::string GetOutputFile(void) const { return kOutputFile ; }
    std::string GetInputFile(void) const { return kInputFile ; }
    std::string GetXSecFile(void) const { return kXSecFile ; }

    virtual ~ConfigureI();
      
  protected: 
    bool InitializeFiducial(void) ;
    
    // Members
    bool kIsData = false ; // Is data (class?)
    bool kUseAllSectors = false ; // Are there any dead sectors?
    bool kApplyFiducial = true ; // Set to false to remove fiducial cuts
    bool kApplyEFiducial = true ; // Set to false to remove fiducial cuts
    bool kApplyHadFiducial = true ; // Set to false to remove fiducial cuts
    bool kApplyAccWeights = true ; // Set to false to ignore acceptance weights 
    bool kApplyMottWeight = true ; // Apply mott xsec convertion
    bool kApplyReso = true ; // Apply particle resolution
    bool kApplyPhiOpeningAngle = false ; // Apply PhiOpening Angle check
    bool kUsePhiThetaBand = false ; // Use a restricted Phi Theta band
    bool kApplyThetaSlice = false; // Use a theta slice
    bool kApplyGoodSectorPhiSlice = false ; 
    bool kApplyMomCut = true ; // Apply min momentum cuts on particles
    bool kQ2Cut = true ; // Apply Q2 cut
    bool kWCut = true ; // Apply W2 cut
    bool kOutEMomCut = true ; // Apply outgoing E cut
    bool kIsElectron = true ; // Is EM data  
    double koffset = 0 ;  // ofset for oscillation studies
    bool kSubtractBkg = false ; // Apply background correction
    bool kNoFSI = false ;

    double kEBeam = 1.161 ; 
    unsigned int kTargetPdg = 1000060120 ;
    unsigned int kNEvents = 0;
    unsigned int kFirstEvent = 0 ; 
    unsigned int kMult_signal = 0;

    Fiducial * kFiducialCut = nullptr ;
    TF1 * kElectronFit = nullptr ; 

    // Topology
    std::map<int,unsigned int> kTopology_map ; // Pdg, multiplicity
    unsigned int kMaxBkgMult = 2 ; 
    unsigned int kNRotations = 100; 

    // Histogram configurables
    std::vector< std::string > kObservables ;
    std::vector< unsigned int > kNBins ;
    std::vector<std::vector<double>> kRanges ; 
    std::string kInputFile ;
    std::string kOutputFile = "";
    std::string kXSecFile = "";
    bool kNormalize = true ; // Normalize histograms to cross section
    bool kApplyCorrWeights = true ; // Set to false to ignore correction weights to be applied to the histograms
    bool kDebugBkg = false ; 

    // Information for output file
    std::unique_ptr<TFile> kOutFile ;
    std::unique_ptr<TTree> kAnalysisTree ; 
    std::vector<TH1D*> kHistograms ; 

    // Configuration validity checks:
    bool kIsDataLoaded = false ;
    bool kIsConfigured = true ; 

  };
}

#endif
