/** 
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/AnalysisI.h"
#include "conf/FiducialCutI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/DetectorUtils.h"
#include "utils/Utils.h"

using namespace e4nu; 

AnalysisI::AnalysisI( ) { this->Initialize(); }

AnalysisI::AnalysisI( const std::string input_file ) : BackgroundI( input_file ) { this->Initialize(); }

AnalysisI::AnalysisI( const double EBeam, const unsigned int TargetPdg ) : BackgroundI( EBeam, TargetPdg ) { this->Initialize();}    

AnalysisI::~AnalysisI() { this->Initialize();}

bool AnalysisI::Analyse( Event & event ) {

  TLorentzVector in_mom = event.GetInLepton4Mom() ;
  TLorentzVector out_mom = event.GetOutLepton4Mom() ;

  // Step 1 : Apply generic cuts
  // Check run is correct
  double EBeam = GetConfiguredEBeam() ; 
  if ( in_mom.E() != EBeam ) {
    std::cout << " Electron energy is " << in_mom.E() << " instead of " << EBeam << "GeV. Configuration failed. Exit" << std::endl;
    exit(11); 
  }
  
  if ( (unsigned int) event.GetTargetPdg() != GetConfiguredTarget() ) {
    std::cout << "Target is " << event.GetTargetPdg() << " instead of " << GetConfiguredTarget() << ". Configuration failed. Exit" << std::endl;
    exit(11); 
  }

  // Check weight is physical
  double wght = event.GetEventWeight() ; 
  if ( wght < 0 || wght > 10 || wght == 0 ) {
    return false ; 
  }

  if( out_mom.Theta() * 180 / TMath::Pi() < GetElectronMinTheta( out_mom ) ) return false ;

  // GENIE coordinate system flipped with respect to CLAS
  if( !IsData() ) out_mom.SetPhi( out_mom.Phi() + TMath::Pi() );

  if( !utils::IsValidSector( out_mom.Phi(), EBeam, UseAllSectors() ) ) return false ; 
  if( !utils::IsValidSector( out_mom.Phi(), EnabledSectors() ) ) return false ; 

  if( ApplyOutElectronCut() ){
    if( out_mom.P() < conf::GetMinMomentumCut( conf::kPdgElectron, EBeam ) ) return false ; 
  }

  if( ApplyThetaSlice() ) {
    if( out_mom.Theta() * 180./TMath::Pi() < conf::kMinEThetaSlice ) return false ; 
    if( out_mom.Theta() * 180./TMath::Pi() > conf::kMaxEThetaSlice ) return false ; 
  }

  if( ApplyPhiOpeningAngle() ) {
    if ( ! conf::ValidPhiOpeningAngle( out_mom.Phi() ) ) return false ;  
  }

  if( ApplyGoodSectorPhiSlice() ) {
    if ( ! conf::GoodSectorPhiSlice( out_mom.Phi() ) ) return false ; 
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

  // Apply Mott Scaling to correct for different coupling
  if ( ApplyMottScaling() ) { 
    event.SetMottXSecWeight() ; 
  }

  // Store analysis record before momentum cuts (0) :
  event.StoreAnalysisRecord(kid_bcuts);

  // Step 3: Cook event
  // Remove particles not specified in topology maps
  // These are ignored in the analysis
  // No Cuts are applied on those
  this->CookEvent( event ) ; 
  
  // Step 4 : Apply momentum cut (detector specific) 
  if( ApplyMomCut() ) { 
    this->ApplyMomentumCut( event ) ; 
  }
  // Store analysis record after momentum cuts:
  event.StoreAnalysisRecord(kid_acuts);

  return true ; 
}


void AnalysisI::ApplyMomentumCut( Event & event ) {

  std::map<int,std::vector<TLorentzVector>> unsmeared_part_map = event.GetFinalParticlesUnCorr4Mom() ;
  TLorentzVector out_mom = event.GetOutLepton4Mom() ;
  // Remove particles below threshold
  for( auto it = unsmeared_part_map.begin() ; it != unsmeared_part_map.end() ; ++it ) {
    std::vector<TLorentzVector> above_th_particles ; 
    for( unsigned int i = 0 ; i < unsmeared_part_map[it->first].size() ; ++i ) {
      // Only store particles above threshold
      if( unsmeared_part_map[it->first][i].P() <= conf::GetMinMomentumCut( it->first, GetConfiguredEBeam() ) )  continue ; 
	 
      // Apply photon cuts for MC and data 
      if( it->first == conf::kPdgPhoton ) {
	if( ! conf::ApplyPhotRadCut( out_mom, unsmeared_part_map[it->first][i] ) ) continue ; 
      }
      above_th_particles.push_back( unsmeared_part_map[it->first][i] ) ;
    }
    unsmeared_part_map[it->first] = above_th_particles ;
  }
  event.SetFinalParticlesKinematics( unsmeared_part_map ) ;
  event.SetFinalParticlesUnCorrKinematics( unsmeared_part_map ) ; 
  
  return ;
}

void AnalysisI::CookEvent( Event & event ) { 
  // Remove particles not specified in topology maps
  // These are ignored in the analysis
  // No Cuts are applied on those
  TLorentzVector out_mom = event.GetOutLepton4Mom() ; 
  std::map<int,std::vector<TLorentzVector>> part_map = event.GetFinalParticlesUnCorr4Mom() ;
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
  event.SetFinalParticlesKinematics( cooked_part_map ) ;
  event.SetFinalParticlesUnCorrKinematics( cooked_part_map ) ;
  return ; 
}

void AnalysisI::PlotBkgInformation( Event event ) {
 
  if( ! GetDebugBkg() ) return ;
 
  // Store plots for Bakcground debugging
  std::map<unsigned int,std::pair<std::vector<int>,double>> AnalysisRecord = event.GetAnalysisRecord() ;
 
  // Signal multiplicity
  unsigned int min_mult = GetMinBkgMult() ;
  unsigned int max_mult = GetMaxBkgMult(); // Max multiplicity specified in conf file
  // Define status ids
  const std::pair<std::vector<int>,double> record_bmomcuts = AnalysisRecord[kid_bcuts] ; // Before mom cuts
  const std::pair<std::vector<int>,double> record_amomcuts = AnalysisRecord[kid_acuts] ; // After mom cuts
  const std::pair<std::vector<int>,double> record_afiducials = AnalysisRecord[kid_fid] ; // After fiducials
  const std::pair<std::vector<int>,double> record_acccorr = AnalysisRecord[kid_acc] ; // Acc Correction

  if( !kHistograms[kid_totestbkg] || !kHistograms[kid_signal] || !kHistograms[kid_tottruebkg] 
      || !kHistograms[kid_2p0pitruebkg] || !kHistograms[kid_1p1pitruebkg] || !kHistograms[kid_2p1pitruebkg] || !kHistograms[kid_1p2pitruebkg] 
      || !kHistograms[kid_2p0piestbkg] || !kHistograms[kid_1p1piestbkg] || !kHistograms[kid_2p1piestbkg] || !kHistograms[kid_1p2piestbkg] ) return ;

  if( (record_afiducials.first).size() > min_mult ) {
    // This is used to estimate the total background contribution 
    kHistograms[kid_totestbkg]->Fill( event.GetObservable("ECal"), - event.GetTotalWeight() ) ; 

    // Store contributions from different multiplicities 
    // Check only direct contribution : 
    unsigned int original_mult = min_mult + 1 ; 
    for( unsigned int j = 0 ; j < max_mult - min_mult ; ++j ) { 
      bool is_m_bkg = true ; 
      if( (AnalysisRecord[original_mult+kid_bkgcorr].first).size() != min_mult ) is_m_bkg = false ; 

      // Fill for direct contributions only
      unsigned int id2 = kid_totestbkg + original_mult - min_mult ; 	
      if( kHistograms[id2] && is_m_bkg ) {
	kHistograms[id2]->Fill( event.GetObservable("ECal"), -event.GetTotalWeight() ) ;
      }

      // Filling for specific topologies:
      // First, we need to get the correct mother pdg list.
      std::vector<int> pdgs ;
      if( (record_afiducials.first).size() == original_mult ) pdgs = record_afiducials.first ; // It comes directly from background event
      else if( (AnalysisRecord[original_mult+1+kid_bkgcorr].first).size() == original_mult ) pdgs= AnalysisRecord[original_mult+1+kid_bkgcorr].first ; // It comes from original_mult + 1 event
      // If none of the two cases above, the bkg event comes directly from a higher multiplicity event ... discard. Only considering direct contributions or corrections

      // Count number of protons and pions
      unsigned int nprotons = 0 ; 
      unsigned int npions = 0 ;
      for ( unsigned int k = 0 ; k < pdgs.size() ; ++k ) { 
	int particle_pdg = pdgs[k];
	if( particle_pdg == conf::kPdgProton ) ++nprotons ; 
	else if ( particle_pdg == conf::kPdgPiP || particle_pdg == conf::kPdgPiM || particle_pdg == conf::kPdgPi0 || particle_pdg == conf::kPdgPhoton ) ++npions ; 
      }

      // Add breakdown in topologies
      if( original_mult == 2 ) {
	// Store according to initial topology
	if( nprotons == 2 && npions == 0 ) kHistograms[kid_2p0piestbkg]->Fill( event.GetObservable("ECal"), -event.GetTotalWeight() ) ;
	if( nprotons == 1 && npions == 1 ) kHistograms[kid_1p1piestbkg]->Fill( event.GetObservable("ECal"), -event.GetTotalWeight() ) ;
      } else if (original_mult == 3 ) { 
	// Store according to initial topology
	if( nprotons == 2 && npions == 1 ) kHistograms[kid_2p1piestbkg]->Fill( event.GetObservable("ECal"), -event.GetTotalWeight() ) ;
	if( nprotons == 1 && npions == 2 ) kHistograms[kid_1p2piestbkg]->Fill( event.GetObservable("ECal"), -event.GetTotalWeight() ) ;
      }

      ++original_mult; 
    }
          
  } else { 
    // These are singal events. They are classified as either true signal or bkg events that contribute to signal after fiducial
    if( (record_afiducials.first).size() == (record_amomcuts.first).size() && (record_acccorr.first).size() == 0 ) {
      kHistograms[kid_signal]->Fill( event.GetObservable("ECal"), event.GetTotalWeight() ) ;
    } else if( (record_afiducials.first).size() == (record_amomcuts.first).size() && (record_acccorr.first).size() != 0 ) {
      kHistograms[kid_acccorr]->Fill( event.GetObservable("ECal"), event.GetTotalWeight() ) ;
    } else { 
      kHistograms[kid_tottruebkg]->Fill( event.GetObservable("ECal"), event.GetTotalWeight() ) ;
     
      // Fill each multiplicity contribution 
      unsigned int id = (record_amomcuts.first).size() - min_mult ; 
      if( kHistograms[kid_tottruebkg+id] ) kHistograms[kid_tottruebkg+id]->Fill( event.GetObservable("ECal"), event.GetTotalWeight() ) ;

      // Count number of protons and pions
      unsigned int nprotons = 0 ; 
      unsigned int npions = 0 ;
      for ( unsigned int k = 0 ; k < (record_amomcuts.first).size() ; ++k ) { 
	int particle_pdg = (record_amomcuts.first)[k];
	if( particle_pdg == conf::kPdgProton ) ++nprotons ; 
	else if ( particle_pdg == conf::kPdgPiP || particle_pdg == conf::kPdgPiM || particle_pdg == conf::kPdgPi0 || particle_pdg == conf::kPdgPhoton ) ++npions ; 
      }

      // Add breakdown in topologies
      if( (record_amomcuts.first).size() == 2 ) {
	// Store according to initial topology
	if( nprotons == 2 && npions == 0 ) kHistograms[kid_2p0pitruebkg]->Fill( event.GetObservable("ECal"), event.GetTotalWeight() ) ;
	if( nprotons == 1 && npions == 1 ) kHistograms[kid_1p1pitruebkg]->Fill( event.GetObservable("ECal"), event.GetTotalWeight() ) ;
      } else if ( (record_amomcuts.first).size() == 3 ) { 
	// Store according to initial topology
	if( nprotons == 2 && npions == 1 ) kHistograms[kid_2p1pitruebkg]->Fill( event.GetObservable("ECal"), event.GetTotalWeight() ) ;
	if( nprotons == 1 && npions == 2 ) kHistograms[kid_1p2pitruebkg]->Fill( event.GetObservable("ECal"), event.GetTotalWeight() ) ;
      }
    }
  }
}

double AnalysisI::GetElectronMinTheta( TLorentzVector emom ) {
  return kElectronFit ->Eval(emom.P()) ; 
}

void AnalysisI::Initialize(void) { 
  double Ebeam = GetConfiguredEBeam() ; 
  
  kElectronFit = new TF1( "myElectronFit", "[0]+[1]/x",0.,0.5);
  if( Ebeam == 1.161 ) { kElectronFit -> SetParameters(17,7) ; }
  if( Ebeam == 2.261 ) { kElectronFit -> SetParameters(16,10.5) ; }
  if( Ebeam == 4.461 ) { kElectronFit -> SetParameters(13.5,15) ; }
  if( !kElectronFit ) kIsConfigured = false ; 
  if( kIsConfigured ) kIsConfigured = InitializeFiducial() ; 
}

bool AnalysisI::Finalise(void) {
  if( kElectronFit ) delete kElectronFit ;
  return true ; 
}
