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
    void Initialize(void);

    // Get generic information about your analysis
    bool IsData(void) const { return kIsData ; } 
    bool IsCLAS6Analysis(void) const { return kIsCLAS6Analysis ; } 
    bool IsCLAS12Analysis(void) const { return kIsCLAS12Analysis ; } 
    bool IsElectronData(void) const { return kIsElectron ; }
    bool IsConfigured(void) const { return kIsConfigured ; }
    unsigned int GetAnalysisTypeID(void) const{ return kAnalysisTypeID ; }
    unsigned int GetNEventsToRun(void) const { return kNEvents ; } 
    unsigned int GetFirstEventToRun(void) const { return kFirstEvent ; } 
    bool ComputeAccCorrection(void) const { return kComputeAccCorr; }
 
    // Get physics information about the analysis      
    double GetConfiguredEBeam(void) const { return kEBeam ; }
    unsigned int GetConfiguredTarget(void) const { return kTargetPdg ; }
    std::map<int,unsigned int> GetTopology(void) const{ return kTopology_map ; } 
    unsigned int GetNTopologyParticles(void) ;    
    
    // Get informtion about cuts:
    bool IsTrueSignal(void) const { return kTrueSignal ; }
    bool UseAllSectors(void) const { return kUseAllSectors ; } 
    bool ApplyFiducial(void) const { return kApplyFiducial ; }
    bool ApplyCorrWeights(void) const { return kApplyCorrWeights ; }
    bool ApplyAccWeights(void) const { return kApplyAccWeights ; }
    bool ApplyMottScaling(void) const { return kApplyMottWeight ; }
    bool ApplyReso(void) const { return kApplyReso ; }
    bool ApplyPhiOpeningAngle(void) const { return kApplyPhiOpeningAngle ; }
    bool UsePhiThetaBand(void) const { return kUsePhiThetaBand ; } 
    bool ApplyThetaSlice(void) const { return kApplyThetaSlice ; } 
    bool ApplyGoodSectorPhiSlice(void) const { return kApplyGoodSectorPhiSlice ; }
    bool GetOffSet(void) const { return koffset ; } 
    bool ApplyQ2Cut(void) const{ return kQ2Cut ; }
    bool ApplyWCut(void) const{ return kWCut ; }
    bool ApplyMomCut(void) const { return kApplyMomCut ; } 
    bool ApplyOutElectronCut(void) const { return kOutEMomCut ; }      
    bool IsNoFSI(void) const { return kNoFSI ; }
    Fiducial * GetFiducialCut(void) { return kFiducialCut ; } 

    // Define setter functions to be able to change the configuration 
    void SetTrueSignal( const bool b ) { kTrueSignal = b ; } 
    void SetApplyFiducial( const bool b ) { kApplyFiducial = b ; }
    void SetApplyAccWeights( const bool b ) { kApplyAccWeights = b ; }
    void SetApplyReso( const bool b ) { kApplyReso = b ; }
    void SetUseAllSectors( const bool b ) { kUseAllSectors = b ; }
    void SetSubtractBkg( const bool b ) { kSubtractBkg = b ; }

    // Get information about the background subtraction method
    unsigned int GetMaxBkgMult(void) const { return kMaxBkgMult ; }
    unsigned int GetMinBkgMult(void) const { return kMult_signal ; }
    unsigned int GetNRotations(void) const { return kNRotations ; } 
    bool GetSubtractBkg(void) const { return kSubtractBkg ; }
    bool GetDebugBkg(void) const { return kDebugBkg ; } 
    
    // Histogram Configurables
    std::vector<std::string> GetObservablesTag(void) const { return kObservables ; }
    std::vector<unsigned int> GetNBins(void) const { return kNBins ; }
    std::vector<std::vector<double>> GetRange(void) const { return kRanges ; } 
    bool NormalizeHist(void) { return kNormalize ; }

    // Output file information
    std::string GetOutputFile(void) const { return kOutputFile ; }
    std::string GetInputFile(void) const { return kInputFile ; }
    std::string GetXSecFile(void) const { return kXSecFile ; }

    // Others 
    void PrintConfiguration(void) const ; 
      
  protected: 
    virtual ~ConfigureI();
    bool InitializeFiducial(void) ;
    
    // Members
    bool kIsData = false ; // Is data
    bool kIsCLAS6Analysis = true ; 
    bool kIsCLAS12Analysis = false ; // Disabled for now
    bool kComputeAccCorr = false ; // Bool used to compute acc correction files
    bool kUseAllSectors = false ; // Are there any dead sectors? It removes sectors 2 and 4
    bool kApplyFiducial = true ; // Set to false to remove fiducial cuts
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
    bool kNoFSI = false ; // Can be use it to turn off FSI from GENIE files
    bool kTrueSignal = false ; // It removes background events before any cuts are applied. This is used for acceptance calculations

    double kEBeam = 1.161 ; 
    unsigned int kTargetPdg = 1000060120 ;
    unsigned int kNEvents = 0;
    unsigned int kFirstEvent = 0 ; 
    unsigned int kMult_signal = 0;

    Fiducial * kFiducialCut = nullptr ;

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

    // Analysis ID 
    unsigned int kAnalysisTypeID = 0 ;  // 0 -> Generic 

    // Analysis Record ID's
    const unsigned int kid_bcuts = 0 ;
    const unsigned int kid_acuts = 1 ;
    const unsigned int kid_fid = 2 ;
    const unsigned int kid_acc = 3 ;
    const unsigned int kid_bkgcorr = 4 ;
    
  };
}

#endif
