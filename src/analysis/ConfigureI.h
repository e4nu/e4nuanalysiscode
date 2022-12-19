/**
 * \info These parameters are configurable 
 * the default values are set here
 **/

#ifndef _CONFIGURABLES_I_H_
#define _CONFIGURABLES_I_H_

#include <vector>
#include <map>
#include "physics/EventI.h"

namespace e4nu { 

  namespace conf {

    class ConfigureI {

    public: 
      // Default constructor
      ConfigureI() ;

      // Configure with input file 
      ConfigureI( const std::string input_file ) ;
      ConfigureI( const double EBeam, const unsigned int TargetPdg ) ;
      
      bool IsConfigured(void) const { return kIsConfigured ; }
      void PrintConfiguration(void) const ; 

      std::map<int,unsigned int> GetTopology(void) const{ return kTopology_map ; } 

      double GetConfiguredEBeam(void) const { return kEBeam ; }
      unsigned int GetConfiguredTarget(void) const { return kTargetPdg ; }

      bool UseAllSectors(void) const { return kUseAllSectors ; } 
      bool ApplyFiducial(void) const { return kApplyFiducial ; }
      bool ApplyAccWeights(void) const { return kApplyAccWeights ; } 
      bool ApplyReso(void) const { return kApplyReso ; }
      bool ApplyPhiOpeningAngle(void) const { return kApplyPhiOpeningAngle ; }
      bool UsePhiThetaBand(void) const { return kUsePhiThetaBand ; } 
      bool ApplyThetaSlice(void) const { return kApplyThetaSlice ; } 
      bool ApplyGoodSectorPhiSlice(void) const { return kApplyGoodSectorPhiSlice ; }
      bool IsData(void) const { return kIsData ; } 
      bool GetOffSet(void) const { return koffset ; } 
      bool ApplyQ2Cut(void) const{ return kQ2Cut ; }
      bool ApplyWCut(void) const{ return kWCut ; }
      bool ApplyOutElectronCut(void) const { return fOutMomCut ; }
      
      bool IsElectronData(void) const { return kIsElectron ; }

      unsigned int GetNEventsToRun(void) const { return kNEvents ; } 
      unsigned int GetFirstEventToRun(void) const { return kFirstEvent ; } 
      unsigned int GetBkgMult(void) const { return kMaxBkgMult ; }

      // Histogram Configurables
      std::vector<std::string> GetObservablesTag(void) const { return kObservables ; }
      std::vector<unsigned int> GetNBins(void) const { return kNBins ; }
      std::vector<std::vector<double>> GetRange(void) const { return kRanges ; } 

      virtual ~ConfigureI();
      
    protected: 
      bool kUseAllSectors = false ; 
      bool kApplyFiducial = true ; 
      bool kApplyAccWeights = true ; 
      bool kApplyReso = true ; 
      bool kApplyPhiOpeningAngle = false ; 
      bool kUsePhiThetaBand = false ;
      bool kApplyThetaSlice = false; 
      bool kApplyGoodSectorPhiSlice = false ; 
      bool kIsData = false ; 
      bool kQ2Cut = true ; 
      bool kWCut = true ; 
      bool fOutMomCut = true ;
      bool kIsElectron = true ; 
      // ofset for oscillation studies
      double koffset = 0 ; 

      // Physics 
      double kEBeam = 1.161 ; 
      unsigned int kTargetPdg = 1000060120 ;

      // Number of events
      unsigned int kNEvents = 0;
      unsigned int kFirstEvent = 0 ; 

      // Histogram configurables
      std::vector< std::string > kObservables ;
      std::vector< unsigned int > kNBins ;
      std::vector<std::vector<double>> kRanges ; 

      // Topology
      std::map<int,unsigned int> kTopology_map ; // Pdg, multiplicity
      unsigned int kMaxBkgMult = 2 ; 

      bool kIsDataLoaded = false ;
      bool kIsConfigured = true ; 
    };
  }
}

#endif
