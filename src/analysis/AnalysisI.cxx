/** 
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/AnalysisI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/DetectorUtils.h"
#include "utils/Utils.h"

using namespace e4nu; 

AnalysisI::AnalysisI( ) { ; }

AnalysisI::AnalysisI( const std::string input_file ) : BackgroundI( input_file ) { if( kIsConfigured ) kIsConfigured = InitializeFiducial() ;  }

AnalysisI::AnalysisI( const double EBeam, const unsigned int TargetPdg ) : BackgroundI( EBeam, TargetPdg ) {;}    

AnalysisI::~AnalysisI() {;}

bool AnalysisI::Analyse( EventI * event ) {

  TLorentzVector in_mom = event -> GetInLepton4Mom() ;
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ;
  double EBeam = in_mom.E() ; 
  
  // Step 1 : Apply generic cuts
  if( !utils::IsValidSector( out_mom.Phi(), EBeam, UseAllSectors() ) ) return false ;   

  if( out_mom.Theta() * 180 / TMath::Pi() < GetElectronMinTheta( out_mom ) ) return false ;

  if( ApplyOutElectronCut() ){
    if( out_mom.P() < conf::GetMinMomentumCut( conf::kPdgElectron, EBeam ) ) return false ; 
  }

  if( ApplyThetaSlice() ) {
    if( out_mom.Theta() * 180./TMath::Pi() < conf::kMinEThetaSlice ) return false ; 
    if( out_mom.Theta() * 180./TMath::Pi() > conf::kMaxEThetaSlice ) return false ; 
  }

  if( ApplyPhiOpeningAngle() ) {
    double phi = out_mom.Phi() ; 
    if ( ! IsData() ) phi += TMath::Pi() ; 
    if ( ! conf::ValidPhiOpeningAngle( phi ) ) return false ;  
  }

  if( ApplyGoodSectorPhiSlice() ) {
    double phi = out_mom.Phi() ; 
    if ( ! IsData() ) phi += TMath::Pi() ; 
    if ( ! conf::GoodSectorPhiSlice( phi ) ) return false ; 
  }

  double reco_Q2 = utils::GetRecoQ2( out_mom, EBeam ) ; 
  double W_var = utils::GetRecoW( out_mom, EBeam ) ;

  if( ApplyQ2Cut() ) {
    double MaxQ2 = 0 ; 
    if( conf::GetQ2Cut( MaxQ2, EBeam ) ) {
      if( reco_Q2 < MaxQ2 ) return false ; 
    }
  }

  if( ApplyWCut() ) {
    double MinW = 0 ; 
    if( conf::GetWCut( MinW, EBeam ) ) {
      if( W_var > MinW ) return false ; 
    }
  }

  // Step 2 : Cook event
  // Remove particles not specified in topology maps
  // These are ignored in the analysis
  // No Cuts are applied on those
  this->CookEvent( event ) ; 
  
  return true ; 
}

void AnalysisI::CookEvent( EventI * event ) { 
  // Remove particles not specified in topology maps
  // These are ignored in the analysis
  // No Cuts are applied on those
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ; 
  std::map<int,std::vector<TLorentzVector>> part_map = event -> GetFinalParticlesUnCorr4Mom() ;
  std::map<int,unsigned int> Topology = GetTopology();
  std::map<int,std::vector<TLorentzVector>> cooked_part_map ; 
  for( auto it = part_map.begin() ; it != part_map.end() ; ++it ) {
    if( Topology.find(it->first) == Topology.end() ) continue ; 
    std::vector<TLorentzVector> topology_particles ;
    for( unsigned int i = 0 ; i < part_map[it->first].size() ; ++i ) {
      topology_particles.push_back( part_map[it->first][i] ) ;
    }
    cooked_part_map[it->first] = topology_particles ;
  }
  event -> SetFinalParticlesKinematics( cooked_part_map ) ;
  event -> SetFinalParticlesUnCorrKinematics( cooked_part_map ) ;
  return ; 
}

bool AnalysisI::Finalise(void) const{
  return true ; 
}
