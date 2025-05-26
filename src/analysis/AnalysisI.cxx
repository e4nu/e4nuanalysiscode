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
#include "conf/ConstantsI.h"
#include "utils/ParticleUtils.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"
#include "conf/ConstantsI.h"

using namespace e4nu;

AnalysisI::AnalysisI( ) { this->Initialize(); }

AnalysisI::AnalysisI( const std::string input_file ) : BackgroundI( input_file ) { this->Initialize(); }

AnalysisI::AnalysisI( const double EBeam, const unsigned int TargetPdg ) : BackgroundI( EBeam, TargetPdg ) { this->Initialize();}

AnalysisI::~AnalysisI() { this->Finalise();}

bool AnalysisI::Analyse( Event & event ) {

  TLorentzVector in_mom = event.GetInLepton4Mom() ;
  TLorentzVector out_mom = event.GetOutLepton4Mom() ;
  
  // Step 1 : Apply generic cuts
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

  if( !utils::IsValidSector( out_mom.Phi(), kEBeam, UseAllSectors() ) ) return false ;
  if( !utils::IsValidSector( out_mom.Phi(), EnabledSectors() ) ) return false ;

  if( ApplyOutElectronCut() ){
    if( out_mom.P() < conf::GetMinMomentumCut( conf::kPdgElectron, kEBeam ) ) return false ;
  }

  if( ApplyThetaSlice() ) {
    if( out_mom.Theta() * 180./TMath::Pi() < GetEThetaSliceMin() ) return false ;
    if( out_mom.Theta() * 180./TMath::Pi() > GetEThetaSliceMax() ) return false ;
  }

  if( ApplyPhiOpeningAngle() ) {
    if ( ! conf::ValidPhiOpeningAngle( out_mom.Phi() ) ) return false ;
  }

  if( ApplyGoodSectorPhiSlice() ) {
    if ( ! conf::GoodSectorPhiSlice( out_mom.Phi() ) ) return false ;
  }

  double reco_Q2 = utils::GetRecoQ2( out_mom, kEBeam ) ;
  double W_var = utils::GetRecoW( out_mom, kEBeam ) ;

  if( ApplyQ2Cut() ) {
    double MaxQ2 = 0 ;
    if( conf::GetQ2Cut( MaxQ2, kEBeam ) ) {
      if( reco_Q2 < MaxQ2 ) return false ;
    }
  }

  if( ApplyWCut() ) {
    double MinW = 0 ;
    if( conf::GetWCut( MinW, kEBeam ) ) {
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

  // (*) Step 4 : Smear before appling mom cuts 
  if( ApplyReso() && !IsData() ) {
    this -> SmearParticles( event ) ; 
  }

  // Step 5 : Apply momentum cut (detector specific)
  if( ApplyMomCut() ) {
    this->ApplyMomentumCut( event ) ;
  }

  // Step 6 : Apply angle cuts, theta for electron protons and pions
  // these are applied to both data and MC
  if ( ! this->ApplyFiducialCutExtra( event ) ) {
    return false ;
  }

  // Store analysis record after momentum cuts and general angle cuts:
  event.StoreAnalysisRecord(kid_acuts);
  
  // Step 7 : Remove true Bkg events if requested before applying fiducial cuts:
  if (IsTrueSignal() && !IsData() )
    {
      std::map<int, unsigned int> Topology = GetTopology();
      std::map<int, std::vector<TLorentzVector>> hadrons = event.GetFinalParticles4Mom();
      for (auto it = Topology.begin(); it != Topology.end(); ++it)
	{
	  if (it->first == conf::kPdgElectron)
	    continue;
	  if (hadrons[it->first].size() != it->second)
	    {
	      return false;
	    }
	}
    }  

  // (*) Step 8 : apply fiducial if requested
  // The detector has gaps where the particles cannot be detected
  // We need to account for these with fiducial cuts
  // It also takes into account angle cuts for particles
  if ( ! this->ApplyFiducialCut( event, ApplyFiducial() ) ) {
    return false;
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
  auto tags = GetObservablesTag() ;

  // Signal multiplicity
  unsigned int min_mult = GetMinBkgMult() ;
  unsigned int max_mult = GetMaxBkgMult(); // Max multiplicity specified in conf file

  // Define status ids
  const std::pair<std::vector<int>,double> record_bmomcuts = AnalysisRecord[kid_bcuts] ; // Before mom cuts
  const std::pair<std::vector<int>,double> record_amomcuts = AnalysisRecord[kid_acuts] ; // After mom cuts
  const std::pair<std::vector<int>,double> record_afiducials = AnalysisRecord[kid_fid] ; // After fiducials
  const std::pair<std::vector<int>,double> record_acccorr = AnalysisRecord[kid_acc] ; // Acc Correction

  for ( unsigned obs_id = 0 ; obs_id < tags.size() ; ++obs_id ) {
    if( !kHistograms[tags[obs_id]][kid_totestbkg] || !kHistograms[tags[obs_id]][kid_signal] || !kHistograms[tags[obs_id]][kid_tottruebkg]
	|| !kHistograms[tags[obs_id]][kid_2p0pitruebkg] || !kHistograms[tags[obs_id]][kid_1p1pitruebkg] || !kHistograms[tags[obs_id]][kid_2p1pitruebkg] || !kHistograms[tags[obs_id]][kid_1p2pitruebkg]
	|| !kHistograms[tags[obs_id]][kid_2p0piestbkg] || !kHistograms[tags[obs_id]][kid_1p1piestbkg] || !kHistograms[tags[obs_id]][kid_2p1piestbkg] || !kHistograms[tags[obs_id]][kid_1p2piestbkg] ) return ;

    if( event.IsBkg() == true ) {
      // This is used to estimate the total background contribution
      kHistograms[tags[obs_id]][kid_totestbkg]->Fill( event.GetObservable(tags[obs_id]), - event.GetTotalWeight() ) ;

      // Store contributions from different multiplicities
      // Check only direct contribution :
      unsigned int original_mult = min_mult + 1 ;
      for( unsigned int j = 0 ; j < max_mult - min_mult ; ++j ) {
	bool is_m_bkg = true ;
	if( (AnalysisRecord[original_mult+kid_bkgcorr].first).size() != min_mult ) is_m_bkg = false ;

	// Fill for direct contributions only
	unsigned int id2 = kid_totestbkg + original_mult - min_mult ;
	if( kHistograms[tags[obs_id]][id2] && is_m_bkg ) {
	  kHistograms[tags[obs_id]][id2]->Fill( event.GetObservable(tags[obs_id]), -event.GetTotalWeight() ) ;
	}

	// Filling for specific topologies:
	// First, we need to get the correct mother pdg list.
	std::vector<int> pdgs ;
	if( (record_afiducials.first).size() == original_mult ) pdgs = record_afiducials.first ; // It comes directly from background event
	else if( (AnalysisRecord[original_mult+1+kid_bkgcorr].first).size() == original_mult ) pdgs= AnalysisRecord[original_mult+1+kid_bkgcorr].first ;
	// It comes from original_mult + 1 event
	// If none of the two cases above, the bkg event comes directly from a higher multiplicity event ... discard.
	// Only considering direct contributions or corrections

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
	  if( nprotons == 2 && npions == 0 ) kHistograms[tags[obs_id]][kid_2p0piestbkg]->Fill( event.GetObservable(tags[obs_id]), -event.GetTotalWeight() ) ;
	  if( nprotons == 1 && npions == 1 ) kHistograms[tags[obs_id]][kid_1p1piestbkg]->Fill( event.GetObservable(tags[obs_id]), -event.GetTotalWeight() ) ;
	} else if (original_mult == 3 ) {
	  // Store according to initial topology
	  if( nprotons == 2 && npions == 1 ) kHistograms[tags[obs_id]][kid_2p1piestbkg]->Fill( event.GetObservable(tags[obs_id]), -event.GetTotalWeight() ) ;
	  if( nprotons == 1 && npions == 2 ) kHistograms[tags[obs_id]][kid_1p2piestbkg]->Fill( event.GetObservable(tags[obs_id]), -event.GetTotalWeight() ) ;
	}

	++original_mult;
      }
    } else {
      // These are singal events. They are classified as either true signal or bkg events that contribute to signal after fiducial
      if( (record_afiducials.first).size() == (record_amomcuts.first).size() && (record_acccorr.first).size() == 0 ) {
	kHistograms[tags[obs_id]][kid_signal]->Fill( event.GetObservable(tags[obs_id]), event.GetTotalWeight() ) ;
      } else {
	kHistograms[tags[obs_id]][kid_tottruebkg]->Fill( event.GetObservable(tags[obs_id]), event.GetTotalWeight() ) ;

	// Fill each multiplicity contribution
	unsigned int id = (record_amomcuts.first).size() - min_mult ;
	if( kHistograms[tags[obs_id]][kid_tottruebkg+id] ) kHistograms[tags[obs_id]][kid_tottruebkg+id]->Fill( event.GetObservable(tags[obs_id]), event.GetTotalWeight() ) ;

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
	  if( nprotons == 2 && npions == 0 ) kHistograms[tags[obs_id]][kid_2p0pitruebkg]->Fill( event.GetObservable(tags[obs_id]), event.GetTotalWeight() ) ;
	  if( nprotons == 1 && npions == 1 ) kHistograms[tags[obs_id]][kid_1p1pitruebkg]->Fill( event.GetObservable(tags[obs_id]), event.GetTotalWeight() ) ;
	} else if ( (record_amomcuts.first).size() == 3 ) {
	  // Store according to initial topology
	  if( nprotons == 2 && npions == 1 ) kHistograms[tags[obs_id]][kid_2p1pitruebkg]->Fill( event.GetObservable(tags[obs_id]), event.GetTotalWeight() ) ;
	  if( nprotons == 1 && npions == 2 ) kHistograms[tags[obs_id]][kid_1p2pitruebkg]->Fill( event.GetObservable(tags[obs_id]), event.GetTotalWeight() ) ;
	}
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
	out_mom_part_shift.SetPhi( part_map[it->first][i].Phi() - fFidAngleShift * TMath::Pi() / 180. ) ;
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

bool AnalysisI::ApplyFiducialCutExtra( Event & event ) {
  // First, we apply it to the electron
  // Apply fiducial cut to electron
  Fiducial * fiducial = GetFiducialCut() ;
  if( ! fiducial ) return true ;

  // apply theta cuts
  TLorentzVector out_mom = event.GetOutLepton4Mom() ;
  if (! fiducial -> FiducialCutExtra(conf::kPdgElectron, GetConfiguredEBeam(), out_mom.Vect(), IsData() ) ) return false ;

  // Apply Fiducial cut for hadrons and photons
  std::map<int,std::vector<TLorentzVector>> part_map = event.GetFinalParticles4Mom() ;
  std::map<int,std::vector<TLorentzVector>> part_map_uncorr = event.GetFinalParticlesUnCorr4Mom() ;
  std::map<int,std::vector<TLorentzVector>> contained_part_map, contained_part_map_uncorr ;
  for( auto it = part_map.begin() ; it != part_map.end() ; ++it ) {
    std::vector<TLorentzVector> visible_part ;
    std::vector<TLorentzVector> visible_part_uncorr ;
    for( unsigned int i = 0 ; i < part_map[it->first].size() ; ++i ) {
      if( ! fiducial -> FiducialCutExtra(it->first, GetConfiguredEBeam(), part_map[it->first][i].Vect(), IsData() ) ) continue ;
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

void AnalysisI::SmearParticles(Event &event)
{
  double EBeam = GetConfiguredEBeam();
  TLorentzVector out_mom = event.GetOutLepton4Mom();

  utils::ApplyResolution(conf::kPdgElectron, out_mom, EBeam);
  event.SetOutLeptonKinematics(out_mom);

  // Apply for other particles
  std::map<int, std::vector<TLorentzVector>> part_map = event.GetFinalParticles4Mom();
  for (std::map<int, std::vector<TLorentzVector>>::iterator it = part_map.begin(); it != part_map.end(); ++it)
    {
      std::vector<TLorentzVector> vtemp;
      for (unsigned int i = 0; i < (it->second).size(); ++i)
	{
	  TLorentzVector temp = (it->second)[i];
	  utils::ApplyResolution(it->first, temp, EBeam);
	  vtemp.push_back(temp);
	}
      part_map[it->first] = vtemp;
    }
  event.SetFinalParticlesKinematics(part_map);
}

bool AnalysisI::StoreTree(Event event)
{
  static bool isFirstCall = true;
  unsigned int InitialNEvents = GetNEventsToRun();
  double ConversionFactor = conf::kConversionFactorCm2ToMicroBarn * TMath::Power(10., -38.);
  int ID = event.GetEventID();
  double BeamE = event.GetInLepton4Mom().E();
  int TargetPdg = event.GetTargetPdg();
  int InLeptonPdg = event.GetInLeptPdg();
  int OutLeptonPdg = event.GetOutLeptPdg();
  double TotWeight = event.GetTotalWeight();
  double EventWght = event.GetEventWeight();

  unsigned int RecoNProtons = event.GetRecoNProtons();
  unsigned int RecoNNeutrons = event.GetRecoNNeutrons();
  unsigned int RecoNPiP = event.GetRecoNPiP();
  unsigned int RecoNPiM = event.GetRecoNPiM();
  unsigned int RecoNPi0 = event.GetRecoNPi0();
  unsigned int RecoNKP = event.GetRecoNKP();
  unsigned int RecoNKM = event.GetRecoNKM();
  unsigned int RecoNK0 = event.GetRecoNK0();
  unsigned int RecoNEM = event.GetRecoNEM();

  // Angle rotated if is MC, otherwise not rotated.
  TLorentzVector out_mom = event.GetOutLepton4Mom();
  double Efl = out_mom.E();
  double pfl = out_mom.P();
  double pflx = out_mom.Px();
  double pfly = out_mom.Py();
  double pflz = out_mom.Pz();
  double pfl_theta = out_mom.Theta() * TMath::RadToDeg();
  unsigned int ElectronSector = utils::GetSector(out_mom.Phi());
  double pfl_phi = out_mom.Phi() * TMath::RadToDeg();
  double pfl_T = utils::GetPT(out_mom.Vect()).Mag();

  double RecoQELEnu = utils::GetQELRecoEnu(out_mom, TargetPdg);
  double RecoEnergyTransfer = utils::GetEnergyTransfer(out_mom, BeamE);
  double Recoq3 = utils::GetRecoq3(out_mom, BeamE).Mag();
  double RecoQ2 = utils::GetRecoQ2(out_mom, BeamE);
  double RecoXBJK = utils::GetRecoXBJK(out_mom, BeamE);
  double RecoW = utils::GetRecoW(out_mom, BeamE);

  double MottXSecScale = event.GetMottXSecWeight();

  std::map<int, unsigned int> topology = GetTopology();
  static bool topology_has_protons = false;
  if (topology.count(conf::kPdgProton) != 0)
    {
      if (topology[conf::kPdgProton] != 0)
	topology_has_protons = true;
    }
  static bool topology_has_pip = false;
  if (topology.count(conf::kPdgPiP) != 0)
    {
      if (topology[conf::kPdgPiP] != 0)
	topology_has_pip = true;
    }

  static bool topology_has_pim = false;
  if (topology.count(conf::kPdgPiM) != 0)
    {
      if (topology[conf::kPdgPiM] != 0)
	topology_has_pim = true;
    }

  unsigned int TopMult = GetNTopologyParticles();

  // flip hadrons phi
  std::map<int, std::vector<TLorentzVector>> hadron_map = event.GetFinalParticles4Mom();
  TLorentzVector p_max(0, 0, 0, 0);
  if (topology_has_protons)
    {
      double max_mom = 0;
      for (unsigned int i = 0; i < hadron_map[conf::kPdgProton].size(); ++i)
	{
	  if (hadron_map[conf::kPdgProton][i].P() > max_mom)
	    {
	      max_mom = hadron_map[conf::kPdgProton][i].P();
	      p_max = hadron_map[conf::kPdgProton][i];
	    }
	}
    }

  double HadronsAngle = 0;
  std::vector<TLorentzVector> particles;
  for (auto it = hadron_map.begin(); it != hadron_map.end(); ++it)
    {
      if ((it->second).size() != 1)
	continue;
      for (unsigned int i = 0; i < (it->second).size(); ++i)
	{
	  particles.push_back((it->second)[i]);
	}
    }
  if (particles.size() == 2)
    {
      HadronsAngle = utils::Angle(particles[0].Vect(), particles[1].Vect()) * TMath::RadToDeg();
    }

  double proton_E = p_max.E();
  double proton_mom = p_max.P();
  double proton_momx = p_max.Px();
  double proton_momy = p_max.Py();
  double proton_momz = p_max.Pz();
  double proton_theta = p_max.Theta() * TMath::RadToDeg();
  double proton_phi = p_max.Phi() * TMath::RadToDeg();
  double ECal = utils::GetECal(out_mom.E(), event.GetFinalParticles4Mom(), TargetPdg);
  double DiffECal = utils::GetECal(out_mom.E(), event.GetFinalParticles4Mom(), TargetPdg) - BeamE;
  double AlphaT = utils::DeltaAlphaT(out_mom.Vect(), p_max.Vect());
  double DeltaPT = utils::DeltaPT(out_mom.Vect(), p_max.Vect()).Mag();
  double DeltaPhiT = utils::DeltaPhiT(out_mom.Vect(), p_max.Vect());
  double HadAlphaT = utils::DeltaAlphaT(out_mom, hadron_map);
  double HadDeltaPT = utils::DeltaPT(out_mom, hadron_map).Mag();
  double HadDeltaPTx = utils::DeltaPTx(out_mom, hadron_map);
  double HadDeltaPTy = utils::DeltaPTy(out_mom, hadron_map);
  double HadDeltaPhiT = utils::DeltaPhiT(out_mom, hadron_map);
  double InferedNucleonMom = utils::InferedNucleonMom(BeamE, out_mom, hadron_map, TargetPdg);

  TLorentzVector pip_max(0, 0, 0, 0);
  if (topology_has_pip)
    {
      double max_mom = 0;
      for (unsigned int i = 0; i < hadron_map[conf::kPdgPiP].size(); ++i)
	{
	  if (hadron_map[conf::kPdgPiP][i].P() > max_mom)
	    {
	      max_mom = hadron_map[conf::kPdgPiP][i].P();
	      pip_max = hadron_map[conf::kPdgPiP][i];
	    }
	}
    }

  double pip_E = pip_max.E();
  double pip_mom = pip_max.P();
  double pip_momx = pip_max.Px();
  double pip_momy = pip_max.Py();
  double pip_momz = pip_max.Pz();
  double pip_theta = pip_max.Theta() * TMath::RadToDeg();
  double pip_phi = pip_max.Phi() * TMath::RadToDeg();

  TLorentzVector pim_max(0, 0, 0, 0);
  if (topology_has_pim)
    {
      double max_mom = 0;
      for (unsigned int i = 0; i < hadron_map[conf::kPdgPiM].size(); ++i)
	{
	  if (hadron_map[conf::kPdgPiM][i].P() > max_mom)
	    {
	      max_mom = hadron_map[conf::kPdgPiM][i].P();
	      pim_max = hadron_map[conf::kPdgPiM][i];
	    }
	}
    }

  // Adler angles
  double AdlerAngleThetaP = utils::GetAdlerAngleTheta(BeamE, out_mom, hadron_map, conf::kPdgProton) * TMath::RadToDeg();
  double AdlerAnglePhiP = utils::GetAdlerAnglePhi(BeamE, out_mom, hadron_map, conf::kPdgProton) * TMath::RadToDeg();
  double AdlerAngleThetaPi = utils::GetAdlerAngleTheta(BeamE, out_mom, hadron_map, abs(conf::kPdgPiM)) * TMath::RadToDeg();
  double AdlerAnglePhiPi = utils::GetAdlerAnglePhi(BeamE, out_mom, hadron_map, abs(conf::kPdgPiM)) * TMath::RadToDeg();

  // Angle between q and had system
  double Angleqvshad = utils::Angle(utils::GetRecoq3(out_mom, BeamE), utils::TotHadron(hadron_map).Vect()) * TMath::RadToDeg();

  double HadSystemMass = utils::HadSystemMass(hadron_map);
  double pim_E = pim_max.E();
  double pim_mom = pim_max.P();
  double pim_momx = pim_max.Px();
  double pim_momy = pim_max.Py();
  double pim_momz = pim_max.Pz();
  double pim_theta = pim_max.Theta() * TMath::RadToDeg();
  double pim_phi = pim_max.Phi() * TMath::RadToDeg();

  double MissingEnergy = utils::Missing4Momenta(BeamE, out_mom, hadron_map, TargetPdg).E();
  double MissingMomentum = utils::Missing4Momenta(BeamE, out_mom, hadron_map, TargetPdg).P();
  double MissingAngle = utils::Missing4Momenta(BeamE, out_mom, hadron_map, TargetPdg).Theta() * TMath::RadToDeg();
  double MissingTransMomentum = utils::GetPT(utils::Missing4Momenta(BeamE, out_mom, hadron_map, TargetPdg).Vect()).Mag();

  TLorentzVector pi_mom(0, 0, 0, 0);
  if (topology_has_pip)
    {
      double max_mom = 0;
      for (unsigned int i = 0; i < hadron_map[conf::kPdgPiP].size(); ++i)
	{
	  if (hadron_map[conf::kPdgPiP][i].P() > max_mom)
	    {
	      max_mom = hadron_map[conf::kPdgPiP][i].P();
	      pi_mom = hadron_map[conf::kPdgPiP][i];
	    }
	}
    }
  else if (topology_has_pim)
    {
      double max_mom = 0;
      for (unsigned int i = 0; i < hadron_map[conf::kPdgPiM].size(); ++i)
	{
	  if (hadron_map[conf::kPdgPiM][i].P() > max_mom)
	    {
	      max_mom = hadron_map[conf::kPdgPiM][i].P();
	      pi_mom = hadron_map[conf::kPdgPiM][i];
	    }
	}
    }

  double RecoEvPion = utils::GetRecoEvPionProduction(out_mom, pi_mom);
  double RecoWPion = utils::GetRecoWPionProduction(out_mom, pi_mom);
  double ElectronPT = utils::GetPT(out_mom.Vect()).Mag();
  double PionPT = utils::GetPT(pi_mom.Vect()).Mag();

  bool IsBkg = event.IsBkg();

  // Store name of files used
  const char *InputROOTFile = kInputFile.c_str();
  const char *OutputROOTFile = kOutputFile.c_str();

  if (isFirstCall == true)
    {
      kAnalysisTree->Branch("InputROOTFile", &InputROOTFile, "InputROOTFile/C");
      kAnalysisTree->Branch("OutputROOTFile", &OutputROOTFile, "OutputROOTFile/C");
      kAnalysisTree->Branch("InitialNEvents", &InitialNEvents, "InitialNEvents/I");
      kAnalysisTree->Branch("ConversionFactor", &ConversionFactor, "ConversionFactor/D");
      kAnalysisTree->Branch("ID", &ID, "ID/I");
      kAnalysisTree->Branch("TargetPdg", &TargetPdg, "TargetPdg/I");
      kAnalysisTree->Branch("InLeptonPdg", &InLeptonPdg, "InLeptonPdg/I");
      kAnalysisTree->Branch("OutLeptonPdg", &OutLeptonPdg, "OutLeptonPdg/I");
      kAnalysisTree->Branch("BeamE", &BeamE, "BeamE/D");
      kAnalysisTree->Branch("TotWeight", &TotWeight, "TotWeight/D");
      kAnalysisTree->Branch("EventWght", &EventWght, "EventWght/D");
      kAnalysisTree->Branch("MottXSecScale", &MottXSecScale, "MottXSecScale/D");
      kAnalysisTree->Branch("IsBkg", &IsBkg, "IsBkg/O");
      kAnalysisTree->Branch("RecoNProtons", &RecoNProtons, "RecoNProtons/I");
      kAnalysisTree->Branch("RecoNNeutrons", &RecoNNeutrons, "RecoNNeutrons/I");
      kAnalysisTree->Branch("RecoNPiP", &RecoNPiP, "RecoNPiP/I");
      kAnalysisTree->Branch("RecoNPiM", &RecoNPiM, "RecoNPiM/I");
      kAnalysisTree->Branch("RecoNPi0", &RecoNPi0, "RecoNPi0/I");
      kAnalysisTree->Branch("RecoNKP", &RecoNKP, "RecoNKP/I");
      kAnalysisTree->Branch("RecoNKM", &RecoNKM, "RecoNKM/I");
      kAnalysisTree->Branch("RecoNK0", &RecoNK0, "RecoNK0/I");
      kAnalysisTree->Branch("RecoNEM", &RecoNEM, "RecoNEM/I");
      kAnalysisTree->Branch("TopMult", &TopMult, "TopMult/I");
      kAnalysisTree->Branch("Efl", &Efl, "Efl/D");
      kAnalysisTree->Branch("pfl", &pfl, "pfl/D");
      kAnalysisTree->Branch("pflx", &pflx, "pflx/D");
      kAnalysisTree->Branch("pfly", &pfly, "pfly/D");
      kAnalysisTree->Branch("pflz", &pflz, "pflz/D");
      kAnalysisTree->Branch("pfl_theta", &pfl_theta, "pfl_theta/D");
      kAnalysisTree->Branch("pfl_phi", &pfl_phi, "pfl_phi/D");
      kAnalysisTree->Branch("pfl_T", &pfl_T, "pfl_T/D");
      kAnalysisTree->Branch("RecoQELEnu", &RecoQELEnu, "RecoQELEnu/D");
      kAnalysisTree->Branch("RecoEnergyTransfer", &RecoEnergyTransfer, "RecoEnergyTransfer/D");
      kAnalysisTree->Branch("Recoq3", &Recoq3, "Recoq3/D");
      kAnalysisTree->Branch("RecoQ2", &RecoQ2, "RecoQ2/D");
      kAnalysisTree->Branch("RecoW", &RecoW, "RecoW/D");
      kAnalysisTree->Branch("RecoXBJK", &RecoXBJK, "RecoXBJK/D");
      kAnalysisTree->Branch("HadSystemMass", &HadSystemMass, "HadSystemMass/D");
      kAnalysisTree->Branch("ElectronSector", &ElectronSector, "ElectronSector/I");
      kAnalysisTree->Branch("MissingEnergy", &MissingEnergy, "MissingEnergy/D");
      kAnalysisTree->Branch("MissingMomentum", &MissingMomentum, "MissingMomentum/D");
      kAnalysisTree->Branch("MissingTransMomentum", &MissingTransMomentum, "MissingTransMomentum/D");
      kAnalysisTree->Branch("MissingAngle", &MissingAngle, "MissingAngle/D");
      kAnalysisTree->Branch("ECal", &ECal, "ECal/D");
      kAnalysisTree->Branch("DiffECal", &DiffECal, "DiffECal/D");
      kAnalysisTree->Branch("HadronsAngle", &HadronsAngle, "HadronsAngle/D");
      kAnalysisTree->Branch("Angleqvshad", &Angleqvshad, "Angleqvshad/D");
      kAnalysisTree->Branch("RecoEvPion", &RecoEvPion, "RecoEvPion/D");
      kAnalysisTree->Branch("RecoWPion", &RecoWPion, "RecoWPion/D");
      kAnalysisTree->Branch("ElectronPT", &ElectronPT, "ElectronPT/D");
      kAnalysisTree->Branch("PionPT", &PionPT, "PionPT/D");

      if (topology_has_protons)
	{
	  kAnalysisTree->Branch("proton_E", &proton_E, "proton_E/D");
	  kAnalysisTree->Branch("proton_mom", &proton_mom, "proton_mom/D");
	  kAnalysisTree->Branch("proton_momx", &proton_momx, "proton_momx/D");
	  kAnalysisTree->Branch("proton_momy", &proton_momy, "proton_momy/D");
	  kAnalysisTree->Branch("proton_momz", &proton_momz, "proton_momz/D");
	  kAnalysisTree->Branch("proton_theta", &proton_theta, "proton_theta/D");
	  kAnalysisTree->Branch("proton_phi", &proton_phi, "proton_phi/D");
	  kAnalysisTree->Branch("AlphaT", &AlphaT, "AlphaT/D");
	  kAnalysisTree->Branch("DeltaPT", &DeltaPT, "DeltaPT/D");
	  kAnalysisTree->Branch("DeltaPhiT", &DeltaPhiT, "DeltaPhiT/D");
	  kAnalysisTree->Branch("AdlerAngleThetaP", &AdlerAngleThetaP, "AdlerAngleThetaP/D");
	  kAnalysisTree->Branch("AdlerAnglePhiP", &AdlerAnglePhiP, "AdlerAnglePhiP/D");
	}

      if (topology_has_pip)
	{
	  kAnalysisTree->Branch("pip_E", &pip_E, "pip_E/D");
	  kAnalysisTree->Branch("pip_mom", &pip_mom, "pip_mom/D");
	  kAnalysisTree->Branch("pip_momx", &pip_momx, "pip_momx/D");
	  kAnalysisTree->Branch("pip_momy", &pip_momy, "pip_momy/D");
	  kAnalysisTree->Branch("pip_momz", &pip_momz, "pip_momz/D");
	  kAnalysisTree->Branch("pip_theta", &pip_theta, "pip_theta/D");
	  kAnalysisTree->Branch("pip_phi", &pip_phi, "pip_phi/D");
	}

      if (topology_has_pim)
	{
	  kAnalysisTree->Branch("pim_E", &pim_E, "pim_E/D");
	  kAnalysisTree->Branch("pim_mom", &pim_mom, "pim_mom/D");
	  kAnalysisTree->Branch("pim_momx", &pim_momx, "pim_momx/D");
	  kAnalysisTree->Branch("pim_momy", &pim_momy, "pim_momy/D");
	  kAnalysisTree->Branch("pim_momz", &pim_momz, "pim_momz/D");
	  kAnalysisTree->Branch("pim_theta", &pim_theta, "pim_theta/D");
	  kAnalysisTree->Branch("pim_phi", &pim_phi, "pim_phi/D");
	}

      if (topology_has_pip || topology_has_pim)
	{
	  kAnalysisTree->Branch("AdlerAngleThetaPi", &AdlerAngleThetaPi, "AdlerAngleThetaPi/D");
	  kAnalysisTree->Branch("AdlerAnglePhiPi", &AdlerAnglePhiPi, "AdlerAnglePhiPi/D");
	}

      kAnalysisTree->Branch("HadAlphaT", &HadAlphaT, "HadAlphaT/D");
      kAnalysisTree->Branch("HadDeltaPT", &HadDeltaPT, "HadDeltaPT/D");
      kAnalysisTree->Branch("HadDeltaPTx", &HadDeltaPTx, "HadDeltaPTx/D");
      kAnalysisTree->Branch("HadDeltaPTy", &HadDeltaPTy, "HadDeltaPTy/D");
      kAnalysisTree->Branch("HadDeltaPhiT", &HadDeltaPhiT, "HadDeltaPhiT/D");
      kAnalysisTree->Branch("InferedNucleonMom", &InferedNucleonMom, "InferedNucleonMom/D");

      isFirstCall = false;
    }

  // Fill
  kAnalysisTree->Fill();

  return true;
}

bool AnalysisI::Finalise(void) {
  //  if( kElectronFit ) delete kElectronFit ;
  return true ;
}
