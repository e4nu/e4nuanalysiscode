/**
 * \info These parameters are configurable
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/AnalysisI.h"
#include "conf/ConstantsI.h"
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

AnalysisI::~AnalysisI() { this->Finalise();}

bool AnalysisI::Analyse( Event & event ) {

  TLorentzVector in_mom = event.GetInLepton4Mom() ;
  TLorentzVector out_mom = event.GetOutLepton4Mom() ;

  // Step 1 : Apply generic cuts
  // Check run is correct
  double EBeam = GetConfiguredEBeam() ;
  /*  if ( in_mom.E() != EBeam ) {
      std::cout << " Electron energy is " << in_mom.E() << " instead of " << EBeam << "GeV. Configuration failed. Exit" << std::endl;
      exit(11);
      }*/

  if ( (unsigned int) event.GetTargetPdg() != GetConfiguredTarget() ) {
    std::cout << "Target is " << event.GetTargetPdg() << " instead of " << GetConfiguredTarget() << ". Configuration failed. Exit" << std::endl;
    exit(11);
  }

  // Check weight is physical
  double wght = event.GetEventWeight() ;
  if ( wght < 0 || wght > 10 || wght == 0 ) {
    return false ;
  }

  if( kElectronFit && out_mom.Theta() * 180 / TMath::Pi() < GetElectronMinTheta( out_mom ) ) return false ;

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

  // Step 5: Apply fiducials
  // The detector has gaps where the particles cannot be detected
  // We need to account for these with fiducial cuts
  // It also takes into account angle cuts for particles
  if ( ! this->ApplyFiducialCut( event, ApplyFiducial() ) ) {
    return false ;
  }

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
  //  if( part_map.find(conf::kPdgPi0) != part_map.end() ) std::cout << "NOTICE: Pi0 present in generation file. They should be decayed" <<std::endl;
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

  // Store plots for Background debugging
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

  if( event.IsBkg() == true ) {
    // This is used to estimate the total background contribution
    kHistograms[kid_totestbkg]->Fill( event.GetObservable("ECal"), - event.GetTotalWeight() ) ;

    ///    std::cout << event.GetTotalWeight()  << std::endl;
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
  if( !kElectronFit ) return -999. ;
  return kElectronFit ->Eval(emom.P()) ;
}

void AnalysisI::Initialize(void) {
  double Ebeam = GetConfiguredEBeam() ;

  kElectronFit = new TF1( "myElectronFit", "[0]+[1]/x",0.,0.5);
  if( !kElectronFit ) kIsConfigured = false ;
  else {
    if( Ebeam == 1.161 ) { kElectronFit -> SetParameters(17,7) ; }
    else if( Ebeam == 2.261 ) { kElectronFit -> SetParameters(16,10.5) ; }
    else if( Ebeam == 4.461 ) { kElectronFit -> SetParameters(13.5,15) ; }
    else { kElectronFit -> SetParameters(13.5,15) ; }
    kIsConfigured = InitializeFiducial() ;
  }

}

bool AnalysisI::ApplyFiducialCut( Event & event, bool apply_fiducial ) {
  // First, we apply it to the electron
  // Apply fiducial cut to electron
  Fiducial * fiducial = GetFiducialCut() ;
  if( ! fiducial ) return true ;

  TLorentzVector out_mom = event.GetOutLepton4Mom() ;
  if (! fiducial -> FiducialCut(conf::kPdgElectron, GetConfiguredEBeam(), out_mom.Vect(), IsData(), apply_fiducial ) ) return false ;

  // Apply Fiducial cut for hadrons and photons
  std::map<int,std::vector<TLorentzVector>> part_map = event.GetFinalParticles4Mom() ;
  std::map<int,std::vector<TLorentzVector>> part_map_uncorr = event.GetFinalParticlesUnCorr4Mom() ;
  std::map<int,std::vector<TLorentzVector>> contained_part_map, contained_part_map_uncorr ;
  for( auto it = part_map.begin() ; it != part_map.end() ; ++it ) {
    std::vector<TLorentzVector> visible_part ;
    std::vector<TLorentzVector> visible_part_uncorr ;
    for( unsigned int i = 0 ; i < part_map[it->first].size() ; ++i ) {
      if( ! fiducial -> FiducialCut(it->first, GetConfiguredEBeam(), part_map[it->first][i].Vect(), IsData(), apply_fiducial ) ) continue ;

      // Consider case in which we shift the fiducial
      // Instead of changing the fiducial limits, we change the phi of each particle
      // To shrink the fiducial we must move the particle closer to the edge
      // Instead of chekcing which edge is closer, we compute both shifts... and check if it is still there
      TLorentzVector out_mom_part_shift = part_map[it->first][i] ;
      if( fFidAngleShift != 0 ) {
	out_mom_part_shift.SetPhi( part_map[it->first][i].Phi() + fFidAngleShift * TMath::Pi() / 180. ) ;
	if( ! fiducial -> FiducialCut(it->first, GetConfiguredEBeam(), out_mom_part_shift.Vect(), IsData(), apply_fiducial ) ) continue ;
	out_mom_part_shift.SetPhi( out_mom.Phi() - fFidAngleShift * TMath::Pi() / 180. ) ;
	if( ! fiducial -> FiducialCut(it->first, GetConfiguredEBeam(), out_mom_part_shift.Vect(), IsData(), apply_fiducial ) ) continue ;
      }

      visible_part.push_back( part_map[it->first][i] ) ;
      visible_part_uncorr.push_back( part_map_uncorr[it->first][i] ) ;
    }
    if( visible_part.size() == 0 ) continue ;
    contained_part_map[it->first] = visible_part ;
    contained_part_map_uncorr[it->first] = visible_part_uncorr ;
  }

  // Store changes in event after fiducial cut
  event.SetFinalParticlesKinematics( contained_part_map ) ;
  event.SetFinalParticlesUnCorrKinematics( contained_part_map_uncorr ) ;

  return true ;
}

bool AnalysisI::Finalise(void) {
  //  if( kElectronFit ) delete kElectronFit ;
  return true ;
}


void AnalysisI::LoadFinalObservables( Event event ){
  // Access observables
  fID = event.GetEventID() ;
  fTargetPdg = event.GetTargetPdg() ;
  fInLeptonPdg = event.GetInLeptPdg() ;
  fOutLeptonPdg = event.GetOutLeptPdg() ;
  fTotWeight = event.GetTotalWeight() ;
  fAccWght = event.GetAccWght() ;
  fEventWght = event.GetEventWeight() ;
  fBeamE = event.GetInLepton4Mom().E() ;

  fCC = event.IsCC();
  fNC = event.IsNC();
  fEM = event.IsEM();
  fQEL = event.IsQEL();
  fRES = event.IsRES();
  fMEC = event.IsMEC();
  fDIS = event.IsDIS();

  fTrueNProtons = event.GetTrueNProtons() ;
  fTrueNNeutrons = event.GetTrueNNeutrons();
  fTrueNPiP = event.GetTrueNPiP();
  fTrueNPiM = event.GetTrueNPiM();
  fTrueNPi0 = event.GetTrueNPi0();
  fTrueNKP = event.GetTrueNKP();
  fTrueNKM = event.GetTrueNKM();
  fTrueNK0 = event.GetTrueNK0();
  fTrueNEM = event.GetTrueNEM();
  fTrueNOther = event.GetTrueNOther();

  fRecoNProtons = event.GetRecoNProtons() ;
  fRecoNNeutrons = event.GetRecoNNeutrons();
  fRecoNPiP = event.GetRecoNPiP();
  fRecoNPiM = event.GetRecoNPiM();
  fRecoNPi0 = event.GetRecoNPi0();
  fRecoNKP = event.GetRecoNKP();
  fRecoNKM = event.GetRecoNKM();
  fRecoNK0 = event.GetRecoNK0();
  fRecoNEM = event.GetRecoNEM();

  fTrueQ2s = event.GetTrueQ2s();
  fTrueWs = event.GetTrueWs();
  fTruexs = event.GetTruexs();
  fTrueys = event.GetTrueys();
  fTrueQ2 = event.GetTrueQ2();
  fTrueW = event.GetTrueW();
  fTruex = event.GetTruex();
  fTruey = event.GetTruey();
  kresid = event.GetRESID();

  TLorentzVector out_mom = event.GetOutLepton4Mom();
  fEfl = out_mom.E();
  fpfl = out_mom.P();
  fpflx = out_mom.Px();
  fpfly = out_mom.Py();
  fpflz = out_mom.Pz();
  fpfl_theta = out_mom.Theta() * TMath::RadToDeg() ;
  out_mom.SetPhi( out_mom.Phi() + TMath::Pi() ) ;
  fElectronSector = utils::GetSector( out_mom.Phi() ) ;
  fpfl_phi = out_mom.Phi() * TMath::RadToDeg() ;

  fRecoQELEnu = utils::GetQELRecoEnu( out_mom, fTargetPdg ) ;
  fRecoEnergyTransfer = utils::GetEnergyTransfer( out_mom, fBeamE ) ;
  fRecoq3 = utils::GetRecoq3( out_mom, fBeamE ).Mag() ;
  fRecoQ2 = utils::GetRecoQ2( out_mom, fBeamE ) ;
  fRecoXBJK = utils::GetRecoXBJK( out_mom, fBeamE ) ;
  fRecoW = utils::GetRecoW(out_mom, fBeamE ) ;

  fMottXSecScale = event.GetMottXSecWeight();

  std::map<int,unsigned int> topology = GetTopology() ;
  static bool topology_has_protons = false ;
  if( topology.count(conf::kPdgProton) != 0 ) {
    if( topology[conf::kPdgProton] != 0 ) topology_has_protons = true ;
  }
  static bool topology_has_pip = false ;
  if( topology.count(conf::kPdgPiP) != 0 ) {
    if( topology[conf::kPdgPiP] != 0 ) topology_has_pip = true ;
  }

  static bool topology_has_pim = false ;
  if( topology.count(conf::kPdgPiM) != 0 ) {
    if( topology[conf::kPdgPiM] != 0 ) topology_has_pim = true ;
  }

  fTopMult = GetNTopologyParticles();

  // flip hadrons phi
  std::map<int,std::vector<TLorentzVector>> hadron_map = event.GetFinalParticles4Mom();
  for( auto it = hadron_map.begin() ; it!=hadron_map.end() ; ++it ) {
    for( unsigned int i = 0 ; i < (it->second).size() ; ++i ) {
      TLorentzVector had_4mom = (it->second)[i] ;
      had_4mom.SetPhi( (it->second)[i].Phi() + TMath::Pi() );
      (it->second)[i] = had_4mom;
    }
  }

  TLorentzVector p_max(0,0,0,0) ;
  if( topology_has_protons ) {
    double max_mom = 0 ;
    for( unsigned int i = 0 ; i < hadron_map[conf::kPdgProton].size() ; ++i ) {
      if( hadron_map[conf::kPdgProton][i].P() > max_mom ) {
	max_mom = hadron_map[conf::kPdgProton][i].P() ;
	p_max = hadron_map[conf::kPdgProton][i] ;
      }
    }
  }

  std::vector<TLorentzVector> particles;
  for( auto it = hadron_map.begin() ; it!=hadron_map.end() ; ++it ) {
    if( (it->second).size() != 1 ) continue ;
    for( unsigned int i = 0 ; i < (it->second).size() ; ++i ) {
      particles.push_back((it->second)[i]) ;
    }
  }
  if( particles.size() == 2 ) {
    kHadronsAngle = utils::Angle( particles[0].Vect(), particles[1].Vect() ) * TMath::RadToDeg() ;
  }
  
  kproton_E = p_max.E() ;
  kproton_mom = p_max.P() ;
  kproton_momx = p_max.Px() ;
  kproton_momy = p_max.Py() ;
  kproton_momz = p_max.Pz() ;
  kproton_theta = p_max.Theta() * TMath::RadToDeg() ;
  kproton_phi = p_max.Phi() * TMath::RadToDeg() ;
  kECal = utils::GetECal( out_mom.E(), event.GetFinalParticles4Mom(), fTargetPdg ) ;
  kDiffECal = utils::GetECal( out_mom.E(), event.GetFinalParticles4Mom(), fTargetPdg ) - fBeamE ;
  kAlphaT = utils::DeltaAlphaT( out_mom.Vect(), p_max.Vect() ) ;
  kDeltaPT = utils::DeltaPT( out_mom.Vect(), p_max.Vect() ).Mag() ;
  kDeltaPhiT = utils::DeltaPhiT( out_mom.Vect(), p_max.Vect() ) ;
  kHadAlphaT = utils::DeltaAlphaT( out_mom, hadron_map ) ;
  kHadDeltaPT = utils::DeltaPT( out_mom, hadron_map ).Mag() ;
  kHadDeltaPTx = utils::DeltaPTx( out_mom, hadron_map ) ;
  kHadDeltaPTy = utils::DeltaPTy( out_mom, hadron_map ) ;
  kHadDeltaPhiT = utils::DeltaPhiT( out_mom, hadron_map ) ;
  kInferedNucleonMom = utils::InferedNucleonMom( fBeamE, out_mom, hadron_map, fTargetPdg ) ;
  
  TLorentzVector pip_max(0,0,0,0) ;
  if( topology_has_pip ) {
    double max_mom = 0 ;
    for( unsigned int i = 0 ; i < hadron_map[conf::kPdgPiP].size() ; ++i ) {
      if( hadron_map[conf::kPdgPiP][i].P() > max_mom ) {
	max_mom = hadron_map[conf::kPdgPiP][i].P() ;
	pip_max = hadron_map[conf::kPdgPiP][i] ;
      }
    }
  }

  fpip_E = pip_max.E() ;
  fpip_mom = pip_max.P() ;
  fpip_momx = pip_max.Px() ;
  fpip_momy = pip_max.Py() ;
  fpip_momz = pip_max.Pz() ;
  fpip_theta = pip_max.Theta() * TMath::RadToDeg();
  fpip_phi = pip_max.Phi() * TMath::RadToDeg() ;

  TLorentzVector pim_max(0,0,0,0) ;
  if( topology_has_pim ) {
    double max_mom = 0 ;
    for( unsigned int i = 0 ; i < hadron_map[conf::kPdgPiM].size() ; ++i ) {
      if( hadron_map[conf::kPdgPiM][i].P() > max_mom ) {
	max_mom = hadron_map[conf::kPdgPiM][i].P() ;
	pim_max = hadron_map[conf::kPdgPiM][i] ;
      }
    }
  }

  //Adler angles
  fAdlerAngleThetaP = utils::GetAdlerAngleTheta( fBeamE, out_mom, hadron_map, conf::kPdgProton) * TMath::RadToDeg();
  fAdlerAnglePhiP = utils::GetAdlerAnglePhi( fBeamE, out_mom, hadron_map, conf::kPdgProton) * TMath::RadToDeg();
  fAdlerAngleThetaPi = utils::GetAdlerAngleTheta( fBeamE, out_mom, hadron_map, abs(conf::kPdgPiM)) * TMath::RadToDeg();
  fAdlerAnglePhiPi = utils::GetAdlerAnglePhi( fBeamE, out_mom, hadron_map, abs(conf::kPdgPiM)) * TMath::RadToDeg() ;

  // Angle between q and had system
  fAngleqvshad = utils::Angle( utils::GetRecoq3( out_mom, fBeamE), utils::TotHadron(hadron_map).Vect() ) * TMath::RadToDeg() ;

  fHadSystemMass = utils::HadSystemMass( hadron_map ) ;
  fpim_E = pim_max.E() ;
  fpim_mom = pim_max.P() ;
  fpim_momx = pim_max.Px() ;
  fpim_momy = pim_max.Py() ;
  fpim_momz = pim_max.Pz() ;
  fpim_theta = pim_max.Theta() * TMath::RadToDeg();
  fpim_phi = pim_max.Phi() * TMath::RadToDeg();

  fMissingEnergy = utils::Missing4Momenta( fBeamE, out_mom, hadron_map ).E();
  fMissingMomentum = utils::Missing4Momenta( fBeamE, out_mom, hadron_map ).P();
  fMissingAngle = utils::Missing4Momenta( fBeamE, out_mom, hadron_map ).Theta() * TMath::RadToDeg() ;

  fIsBkg = event.IsBkg() ;
  fInitialNEvents = GetNEventsToRun() ;
  fConversionFactor = conf::kConversionFactorCm2ToMicroBarn  * TMath::Power(10.,-38.) ;
  fTotalXSec = kXSec ;
  fMCNormalization = fTotalXSec * fConversionFactor / fInitialNEvents ;

  return;
}

double AnalysisI::GetObservable( std::string obs ){
  /*
  double fTotWeight, fAccWght, fEventWght, fBeamE, fTrueQ2s, fTrueWs,fTruexs, fTrueys,fTrueQ2,fTrueW,fTruex,fTruey ;
  double fEfl, fpfl, fpflx, fpfly, fpflz, fpfl_theta, fpfl_phi, fRecoQELEnu, fRecoEnergyTransfer, fRecoq3, fRecoQ2, fRecoXBJK, fRecoW, fMottXSecScale;
  double kHadronsAngle,kproton_E, kproton_mom, kproton_momx, kproton_momy, kproton_momz, kproton_theta, kproton_phi, kECal, kDiffECal, kAlphaT, kDeltaPT;
  double kDeltaPhiT, kHadAlphaT, kHadDeltaPT, kHadDeltaPTx, kHadDeltaPTy, kHadDeltaPhiT, kInferedNucleonMom ;
  double fpip_E, fpip_mom, fpip_momx, fpip_momy, fpip_momz, fpip_theta, fpip_phi;
  double fAdlerAngleThetaP, fAdlerAnglePhiP, fAdlerAngleThetaPi, fAdlerAnglePhiPi, fAngleqvshad, fHadSystemMass ;
  double fpim_E, fpim_mom, fpim_momx, fpim_momy, fpim_momz, fpim_theta, fpim_phi;
  double fMissingEnergy, fMissingMomentum, fMissingAngle, fConversionFactor, fTotalXSec, fMCNormalization;
*/
  if( obs == "ECal" ) {
    return kECal ;
  } 


  return 0 ; 
} 
