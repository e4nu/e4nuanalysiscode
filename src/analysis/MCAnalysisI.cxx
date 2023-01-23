// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "TFile.h"
#include "TDirectoryFile.h"
#include <TRandom3.h>
#include "analysis/MCAnalysisI.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"
#include "utils/DetectorUtils.h"
#include "conf/FiducialCutI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "conf/AnalysisConstantsI.h"

using namespace e4nu ; 

MCAnalysisI::MCAnalysisI() {
  kAcceptanceMap.clear();
  kAccMap.clear();
  kGenMap.clear();
  kAnalysisTree = std::unique_ptr<TTree>( new TTree("MCTree","GENIE Tree") ) ; 
  kMult_signal = GetNTopologyParticles() ; 
  this->Initialize() ;
}

MCAnalysisI::~MCAnalysisI() {
  delete fData;

  kAcceptanceMap.clear();
  kAccMap.clear();
  kGenMap.clear();
}

bool MCAnalysisI::LoadData( void ) {
  if( ! IsConfigured() ) return false ; 

  std::string file = GetInputFile() ; 
  double nevents = GetNEventsToRun() ; 
  double first_event = GetFirstEventToRun() ; 

  if( ! kIsDataLoaded ) { 
    fData = new MCEventHolder( file, first_event, nevents ) ;
    kNEvents = fData->GetNEvents() ; 
    kIsDataLoaded = true ;
  }
  return kIsDataLoaded ; 
}

EventI * MCAnalysisI::GetEvent( const unsigned int event_id ) {
  return fData -> GetEvent(event_id) ; 
}

EventI * MCAnalysisI::GetValidEvent( const unsigned int event_id ) {
  
  MCEvent * event = (MCEvent*) fData -> GetEvent(event_id) ; 
  if( !event ) {
    delete event ; 
    return nullptr ; 
  }

  ++fEventsBeforeCuts ;

  TLorentzVector in_mom = event -> GetInLepton4Mom() ; 
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ; 

  // Check run is correct
  double EBeam = GetConfiguredEBeam() ; 
  if ( in_mom.E() != EBeam ) {
    std::cout << " Electron energy is " << in_mom.E() << " instead of " << EBeam << "GeV. Configuration failed. Exit" << std::endl;
    delete event ;
    exit(11); 
  }

  if ( (unsigned int) event -> GetTargetPdg() != GetConfiguredTarget() ) {
    std::cout << "Target is " << event -> GetTargetPdg() << " instead of " << GetConfiguredTarget() << ". Configuration failed. Exit" << std::endl;
    delete event ;
    exit(11); 
  }

  // Check weight is physical
  double wght = event->GetEventWeight() ; 
  if ( wght < 0 || wght > 10 || wght == 0 ) {
    delete event ;
    return nullptr ; 
  }

  // Apply momentum cut before smearing
  // Get Topology Definition
  std::map<int,std::vector<TLorentzVector>> unsmeared_part_map = event -> GetFinalParticles4Mom() ;

  // Remove particles below threshold
  for( auto it = unsmeared_part_map.begin() ; it != unsmeared_part_map.end() ; ++it ) {
    std::vector<TLorentzVector> above_th_particles ; 
    for( unsigned int i = 0 ; i < unsmeared_part_map[it->first].size() ; ++i ) {
      // Only store particles above threshold
      if( unsmeared_part_map[it->first][i].P() < conf::GetMinMomentumCut( it->first, EBeam ) ) continue ; 
      
      // Apply photon cuts for MC and data 
      if( it->first == conf::kPdgPhoton ) {
	if( !conf::ApplyPhotRadCut( out_mom, unsmeared_part_map[it->first][i] ) ) continue ; 
      }
      above_th_particles.push_back( unsmeared_part_map[it->first][i] ) ;
    }
    unsmeared_part_map[it->first] = above_th_particles ;
  }
  event -> SetFinalParticlesKinematics( unsmeared_part_map ) ;

  // Apply smaring to particles
  if( ApplyReso() ) {
    this -> SmearParticles( event ) ; 
  }

  // Apply fiducial cut to electron
  if( ApplyFiducial() ) {
    if (! kFiducialCut -> EFiducialCut(EBeam, out_mom.Vect() ) ) { delete event ; return nullptr ; } 
  }
  ++fNEventsAfterFiducial;

  // Get Topology Definition
  std::map<int,unsigned int> Topology = GetTopology(); 
  std::map<int,std::vector<TLorentzVector>> part_map = event -> GetFinalParticles4Mom() ;

  // Apply Fiducial cut for hadrons and photons
  if( ApplyFiducial() ) {  
    for( auto it = part_map.begin() ; it != part_map.end() ; ++it ) {
      std::vector<TLorentzVector> visible_part ; 
      for( unsigned int i = 0 ; i < part_map[it->first].size() ; ++i ) {
	if( it->first == conf::kPdgElectron ) {
	  if (! kFiducialCut -> EFiducialCut(EBeam, part_map[it->first][i].Vect() ) ) continue ; 
        } else if ( it->first == conf::kPdgProton ) {
	  if( ! kFiducialCut -> PFiducialCut( EBeam, part_map[it->first][i].Vect() ) ) continue ; 
        } else if ( it->first == conf::kPdgPiP ) {
	  if( ! kFiducialCut -> Pi_phot_fid_united( EBeam, part_map[it->first][i].Vect(), 1 ) ) continue ;
	} else if ( it->first == conf::kPdgPiM ) {
	  if( ! kFiducialCut -> Pi_phot_fid_united( EBeam, part_map[it->first][i].Vect(), -1 ) ) continue ;
	} else if ( it->first == conf::kPdgPhoton ) {
	  if( ! kFiducialCut -> Pi_phot_fid_united( EBeam, part_map[it->first][i].Vect(), 0 ) ) continue ; 
	}
	visible_part.push_back( part_map[it->first][i] ) ; 
      }
      part_map[it->first] = visible_part ; 
    }
  }
  // Store changes in event after fiducial and momentum cuts
  event -> SetFinalParticlesKinematics( part_map ) ; 
  
  // Apply acceptance to all particles 
  double acc_wght = 1 ;
  if( ApplyAccWeights() ) {
    // Electron acceptance
    if( kAccMap[conf::kPdgElectron] && kGenMap[conf::kPdgElectron] ) acc_wght *= utils::GetAcceptanceMapWeight( *kAccMap[conf::kPdgElectron], *kGenMap[conf::kPdgElectron], out_mom ) ; 
    // Others
    for( auto it = Topology.begin() ; it != Topology.end() ; ++it ) {
      if ( part_map.find(it->first) == part_map.end()) continue ;
      if ( it->first == conf::kPdgElectron ) continue ; 
      else { 
	for( unsigned int i = 0 ; i < part_map[it->first].size() ; ++i ) {
	  if( kAccMap[it->first] && kGenMap[it->first] ) acc_wght *= utils::GetAcceptanceMapWeight( *kAccMap[it->first], *kGenMap[it->first], part_map[it->first][i] ) ;
	}
      }
    }
  }
  event->SetAccWght(acc_wght);

  // Signal event
  bool is_signal = false ; 
  //Topology ID
  for( auto it = Topology.begin() ; it != Topology.end() ; ++it ) {
    if( it->first == conf::kPdgElectron ) {
      is_signal = true ; 
    } else if( part_map.count( it->first ) && part_map[it->first].size() == it->second ) {
      is_signal = true ; 
    } else {
      is_signal = false ; 
    }
  }

  if( is_signal ) ++fNEventsAfterTopologyCut ;
  else { // BACKGROUND 
    event->SetIsBkg(true); 
    ++fNBkgEvents ;

    // Get Number of signal particles, "multiplicity"
    unsigned int mult_bkg = event->GetNSignalParticles( part_map, Topology ) ; 

    // Only store background events with multiplicity > mult_signal
    // Also ignore background events above the maximum multiplicity
    if( mult_bkg > kMult_signal && mult_bkg <= GetMaxBkgMult() ) {
      if( fBkg.find(mult_bkg) == fBkg.end() ) {
	std::vector<EventI*> temp ( 1, event ) ;
	fBkg[mult_bkg] = temp ; 
      } else { 
	fBkg[mult_bkg].push_back( event ) ; 
      }
    }
    delete event ; 
    return nullptr; 
  }
  
  StoreTree(event) ;
  
  return event ; 
    
}

void MCAnalysisI::SmearParticles( MCEvent * event ) {
  
  double EBeam = GetConfiguredEBeam() ; 
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ; 

  utils::ApplyResolution( conf::kPdgElectron, out_mom, EBeam ) ; 
  event -> EventI::SetOutLeptonKinematics( out_mom ) ; 
  
  // Apply for other particles
  std::map<int,std::vector<TLorentzVector>> part_map = event -> GetFinalParticles4Mom() ;
  for( std::map<int,std::vector<TLorentzVector>>::iterator it = part_map.begin() ; it != part_map.end() ; ++it ) {
    std::vector<TLorentzVector> vtemp ; 
    for( unsigned int i = 0 ; i < (it->second).size() ; ++i ) { 
      TLorentzVector temp = (it->second)[i] ; 
      utils::ApplyResolution( it->first, temp, EBeam ) ;
      vtemp.push_back(temp) ; 
    }
    part_map[it->first] = vtemp ; 
  }
  event -> EventI::SetFinalParticlesKinematics( part_map ) ; 
  
} 

unsigned int MCAnalysisI::GetNEvents( void ) const {
  return (unsigned int) fData ->GetNEvents() ; 
}

void MCAnalysisI::Initialize() { 

  fData = nullptr ; 

  // Get run configurables
  double EBeam = GetConfiguredEBeam() ; 
  unsigned int Target = GetConfiguredTarget() ;

  if( ApplyFiducial() ) {
    // Initialize fiducial for this run
    kFiducialCut = new Fiducial() ; 
    kFiducialCut -> InitPiMinusFit( EBeam ) ; 
    kFiducialCut -> InitEClimits(); 
    kFiducialCut -> up_lim1_ec -> Eval(60) ;
    kFiducialCut -> SetConstants( conf::GetTorusCurrent( EBeam ), Target , EBeam ) ;
    kFiducialCut -> SetFiducialCutParameters( EBeam ) ;
  }

  // Initialize acceptance map histograms from file
  if( ApplyAccWeights() ) { 
    kAcceptanceMap = conf::GetAcceptanceFileMap2( Target, EBeam ) ; 

    // THESE ONE BELOW CAUSE A SMALL MEMORY LEAK - INVESTIGATE
    kAccMap[conf::kPdgElectron] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgElectron] -> Get("Accepted Particles") ) ) ;
    kAccMap[conf::kPdgProton] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgProton] -> Get("Accepted Particles") ) );
    kAccMap[conf::kPdgPiP] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiP] -> Get("Accepted Particles") ) );
    kAccMap[conf::kPdgPiM] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiM] -> Get("Accepted Particles") ) );

    kGenMap[conf::kPdgElectron] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgElectron] -> Get("Generated Particles") ) );
    kGenMap[conf::kPdgProton] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgProton] -> Get("Generated Particles") ) ) ;
    kGenMap[conf::kPdgPiP] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiP] -> Get("Generated Particles") ) ) ;
    kGenMap[conf::kPdgPiM] = std::unique_ptr<TH3D>( dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiM] -> Get("Generated Particles") ) ) ;
    
    for( auto it = kAccMap.begin() ; it != kAccMap.end(); ++it ) {
      kAccMap[it->first]->SetDirectory(nullptr);
      kGenMap[it->first]->SetDirectory(nullptr);
    }
  }  

  // Get xsec from cross section spline
  std::string target_tag = "" ;
  if( GetConfiguredTarget() == conf::kPdgC12 ) target_tag = "e-_C12" ; 
  if( GetConfiguredTarget() == conf::kPdgHe3 ) target_tag = "e-_He3" ;
  if( GetConfiguredTarget() == conf::kPdgHe4 ) target_tag = "e-_He4" ;
  if( GetConfiguredTarget() == conf::kPdgFe56 ) target_tag = "e-_Fe56" ;
  if( GetConfiguredTarget() == conf::kPdgO16 ) target_tag = "e-_O16" ;

  std::unique_ptr<TFile> xsec_file = std::unique_ptr<TFile>( new TFile( ( GetXSecFile() ).c_str(),"READ") );
  if( !xsec_file ) {
    std::cout << " ERROR: Xsec file does not exist: " << GetXSecFile() << std::endl;
    kIsConfigured = false ; 
    return ; 
  }
  
  TDirectoryFile * xsec_dir = (TDirectoryFile *) xsec_file -> Get(target_tag.c_str());
  if( !xsec_dir ) {
    std::cout << " ERROR: Xsec dir does not exist: " << target_tag << std::endl;
    kIsConfigured = false ; 
    return ; 
  }
    
  TGraph * gxsec  = (TGraph*) xsec_dir  -> Get("tot_em");
  if( !gxsec ) {
    std::cout << " ERROR: Cannot create graph for " << "tot_em" << std::endl;
    kIsConfigured = false; 
    return ; 
  }

  fXSec = gxsec->Eval( GetConfiguredEBeam() ) ; 

  xsec_file->Close();

  gRandom = new TRandom3() ; 
  gRandom->SetSeed(10);

}

bool MCAnalysisI::Finalise( void ) {

  // Normalize
  double domega = 0.01; // sr
  double ConversionFactorCm2ToMicroBarn = TMath::Power(10.,30.); // cm^2 to μbarn

  if ( NormalizeHist() ) {
    for( unsigned int j = 0 ; j < kHistograms.size() ; ++j ) {
      double NBins = kHistograms[j]->GetNbinsX(); 
    
      for (int k = 1; k <= NBins; k++) { 
	double content = kHistograms[j]->GetBinContent(k);
	double error = kHistograms[j]->GetBinError(k);
	double width = kHistograms[j]->GetBinWidth(k);
	double newcontent = content / width;
	double newerror = error / width;
	kHistograms[j]->SetBinContent(k,newcontent);
	kHistograms[j]->SetBinError(k,newerror);
      }

      kHistograms[j]->Scale( fXSec * ConversionFactorCm2ToMicroBarn  * TMath::Power(10.,-38.) / ( GetNEventsToRun() * domega ) );
    }
  }
  std::cout << " Total Number of Events Processed = " << fEventsBeforeCuts << std::endl;
  std::cout << " Total number of true signal events = " << fNEventsAfterTopologyCut << std::endl;
  std::cout << " Events after electron fiducial cut = " << fNEventsAfterFiducial << std::endl;

  return true ; 
}

bool MCAnalysisI::StoreTree(MCEvent * event){
  static int n = true ;

  int ID = event->GetEventID() ; 
  int TargetPdg = event->GetTargetPdg() ;
  int InLeptonPdg = event->GetInLeptPdg() ; 
  int OutLeptonPdg = event->GetOutLeptPdg() ; 
  double TotWeight = event->GetTotalWeight() ; 
  double AccWght = event->GetAccWght() ; 
  double EventWght = event->GetEventWeight() ; 
  double BeamE = event->GetInLepton4Mom().E() ; 

  bool CC = event->IsCC();
  bool NC = event->IsNC();
  bool EM = event->IsEM();
  bool QEL = event->IsQEL();
  bool RES = event->IsRES();
  bool MEC = event->IsMEC();
  bool DIS = event->IsDIS();

  unsigned int TrueNProtons = event->GetTrueNProtons() ; 
  unsigned int TrueNNeutrons = event->GetTrueNNeutrons();
  unsigned int TrueNPiP = event->GetTrueNPiP();
  unsigned int TrueNPiM = event->GetTrueNPiM();
  unsigned int TrueNPi0 = event->GetTrueNPi0();
  unsigned int TrueNKP = event->GetTrueNKP();
  unsigned int TrueNKM = event->GetTrueNKM(); 
  unsigned int TrueNK0 = event->GetTrueNK0();
  unsigned int TrueNEM = event->GetTrueNEM();
  unsigned int TrueNOther = event->GetTrueNOther();

  unsigned int RecoNProtons = event->GetRecoNProtons() ; 
  unsigned int RecoNNeutrons = event->GetRecoNNeutrons();
  unsigned int RecoNPiP = event->GetRecoNPiP();
  unsigned int RecoNPiM = event->GetRecoNPiM();
  unsigned int RecoNPi0 = event->GetRecoNPi0();
  unsigned int RecoNKP = event->GetRecoNKP();
  unsigned int RecoNKM = event->GetRecoNKM(); 
  unsigned int RecoNK0 = event->GetRecoNK0();
  unsigned int RecoNEM = event->GetRecoNEM();

  double TrueQ2s = event->GetTrueQ2s();
  double TrueWs = event->GetTrueWs();
  double Truexs = event->GetTruexs();
  double Trueys = event->GetTrueys();
  double TrueQ2 = event->GetTrueQ2();
  double TrueW = event->GetTrueW();
  double Truex = event->GetTruex();
  double Truey = event->GetTruey();

  TLorentzVector out_mom = event->GetOutLepton4Mom();
  double Efl = out_mom.E();
  double pfl = out_mom.P();
  double pflx = out_mom.Px();
  double pfly = out_mom.Py();
  double pflz = out_mom.Pz();
  double pfl_theta = out_mom.Theta();
  double pfl_phi = out_mom.Phi() + TMath::Pi();
  unsigned int ElectronSector = utils::GetSector( pfl_phi ) ; 

  double RecoQELEnu = utils::GetQELRecoEnu( out_mom, TargetPdg ) ; 
  double RecoEnergyTransfer = utils::GetEnergyTransfer( out_mom, TargetPdg ) ; 
  double Recoq3 = utils::GetRecoq3( out_mom, BeamE ).Mag() ; 
  double RecoQ2 = utils::GetRecoQ2( out_mom, BeamE ) ; 
  double RecoXBJK = utils::GetRecoXBJK( out_mom, BeamE ) ; 
  double RecoW = utils::GetRecoW(out_mom, BeamE ) ;
  double MottXSecScale = event->GetMottXSecWeight();

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

  unsigned int TopMult = GetNTopologyParticles();
  std::map<int,std::vector<TLorentzVector>> hadron_map = event->GetFinalParticles4Mom();
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
  double proton_mom = p_max.P() ; 
  double proton_momx = p_max.Px() ; 
  double proton_momy = p_max.Py() ; 
  double proton_momz = p_max.Pz() ; 
  double proton_theta = p_max.Theta() ; 
  double proton_phi = p_max.Phi() + TMath::Pi() ; 
  double ECal = utils::GetECal( out_mom, p_max, TargetPdg ) ; 
  double AlphaT = utils::DeltaAlphaT( out_mom.Vect(), p_max.Vect() ) ; 
  double DeltaPT = utils::DeltaPT( out_mom.Vect(), p_max.Vect() ).Mag() ; 
  double DeltaPhiT = utils::DeltaPhiT( out_mom.Vect(), p_max.Vect() ) ; 

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
  double pip_mom = pip_max.P() ;
  double pip_momx = pip_max.Px() ;
  double pip_momy = pip_max.Py() ;
  double pip_momz = pip_max.Pz() ;
  double pip_theta = pip_max.Theta() ;
  double pip_phi = pip_max.Phi() + TMath::Pi();

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
  double pim_mom = pim_max.P() ;
  double pim_momx = pim_max.Px() ;
  double pim_momy = pim_max.Py() ;
  double pim_momz = pim_max.Pz() ;
  double pim_theta = pim_max.Theta() ;
  double pim_phi = pim_max.Phi() + TMath::Pi() ;

  bool IsBkg = event->IsBkg() ; 
   
  if( n == true ) {
    kAnalysisTree -> Branch( "ID", &ID, "ID/I"); 
    kAnalysisTree -> Branch( "TargetPdg", &TargetPdg, "TargetPdg/I");
    kAnalysisTree -> Branch( "InLeptonPdg", &InLeptonPdg, "InLeptonPdg/I");
    kAnalysisTree -> Branch( "OutLeptonPdg", &OutLeptonPdg, "OutLeptonPdg/I");
    kAnalysisTree -> Branch( "BeamE", &BeamE, "BeamE/D");
    kAnalysisTree -> Branch( "CC", &CC, "CC/B");
    kAnalysisTree -> Branch( "NC", &NC, "NC/B");
    kAnalysisTree -> Branch( "EM", &EM, "EC/B");
    kAnalysisTree -> Branch( "QEL", &QEL, "QEL/B");
    kAnalysisTree -> Branch( "RES", &RES, "RES/B");
    kAnalysisTree -> Branch( "MEC", &MEC, "MEC/B");
    kAnalysisTree -> Branch( "DIS", &DIS, "DIS/B");
    kAnalysisTree -> Branch( "TrueQ2s", &TrueQ2s, "TrueQ2s/D");
    kAnalysisTree -> Branch( "TrueWs", &TrueWs, "TrueWs/D");
    kAnalysisTree -> Branch( "Truexs", &Truexs, "Truexs/D");
    kAnalysisTree -> Branch( "Trueys", &Trueys, "Trueys/D");
    kAnalysisTree -> Branch( "TrueQ2", &TrueQ2, "TrueQ2/D");
    kAnalysisTree -> Branch( "TrueW", &TrueW, "TrueW/D");
    kAnalysisTree -> Branch( "Truex", &Truex, "Truex/D");
    kAnalysisTree -> Branch( "Truey", &Truey, "Truey/D");
    kAnalysisTree -> Branch( "TotWeight", &TotWeight, "TotWeight/D");
    kAnalysisTree -> Branch( "EventWght", &EventWght, "EventWght/D");
    kAnalysisTree -> Branch( "AccWght", &AccWght, "AccWght/D");
    kAnalysisTree -> Branch( "MottXSecScale", &MottXSecScale, "MottXSecScale/D");
    kAnalysisTree -> Branch( "IsBkg", &IsBkg, "IsBkg/B");
    kAnalysisTree -> Branch( "TrueNProtons", &TrueNProtons, "TrueNProtons/I");
    kAnalysisTree -> Branch( "TrueNNeutrons", &TrueNNeutrons, "TrueNNeutrons/I");
    kAnalysisTree -> Branch( "TrueNPiP", &TrueNPiP, "TrueNPiP/I");
    kAnalysisTree -> Branch( "TrueNPiM", &TrueNPiM, "TrueNPiM/I");
    kAnalysisTree -> Branch( "TrueNPi0", &TrueNPi0, "TrueNPi0/I");
    kAnalysisTree -> Branch( "TrueNKP", &TrueNKP, "TrueNKP/I");
    kAnalysisTree -> Branch( "TrueNKM", &TrueNKM, "TrueNKM/I");
    kAnalysisTree -> Branch( "TrueNK0", &TrueNK0, "TrueNK0/I");
    kAnalysisTree -> Branch( "TrueNEM", &TrueNEM, "TrueNEM/I");
    kAnalysisTree -> Branch( "TrueNOther", &TrueNOther, "TrueNOther/I"); 
    kAnalysisTree -> Branch( "RecoNProtons", &RecoNProtons, "RecoNProtons/I");
    kAnalysisTree -> Branch( "RecoNNeutrons", &RecoNNeutrons, "RecoNNeutrons/I");
    kAnalysisTree -> Branch( "RecoNPiP", &RecoNPiP, "RecoNPiP/I");
    kAnalysisTree -> Branch( "RecoNPiM", &RecoNPiM, "RecoNPiM/I");
    kAnalysisTree -> Branch( "RecoNPi0", &RecoNPi0, "RecoNPi0/I");
    kAnalysisTree -> Branch( "RecoNKP", &RecoNKP, "RecoNKP/I");
    kAnalysisTree -> Branch( "RecoNKM", &RecoNKM, "RecoNKM/I");
    kAnalysisTree -> Branch( "RecoNK0", &RecoNK0, "RecoNK0/I");
    kAnalysisTree -> Branch( "RecoNEM", &RecoNEM, "RecoNEM/I");
    kAnalysisTree -> Branch( "TopMult", &TopMult, "TopMult/I");
    kAnalysisTree -> Branch( "Efl", &Efl, "Efl/D");
    kAnalysisTree -> Branch( "pfl", &pfl, "pfl/D");
    kAnalysisTree -> Branch( "pflx", &pflx, "pflx/D");
    kAnalysisTree -> Branch( "pfly", &pfly, "pfly/D");
    kAnalysisTree -> Branch( "pflz", &pflz, "pflz/D");
    kAnalysisTree -> Branch( "pfl_theta", &pfl_theta, "pfl_theta/D");
    kAnalysisTree -> Branch( "pfl_phi", &pfl_phi, "pfl_phi/D");
    kAnalysisTree -> Branch( "RecoQELEnu", &RecoQELEnu, "RecoQELEnu/D");
    kAnalysisTree -> Branch( "RecoEnergyTransfer", &RecoEnergyTransfer, "RecoEnergyTransfer/D");
    kAnalysisTree -> Branch( "Recoq3", &Recoq3, "Recoq3/D");
    kAnalysisTree -> Branch( "RecoQ2", &RecoQ2, "RecoQ2/D");
    kAnalysisTree -> Branch( "RecoW", &RecoW, "RecoW/D");
    kAnalysisTree -> Branch( "RecoXBJK", &RecoXBJK, "RecoXBJK/D");
    kAnalysisTree -> Branch( "ElectronSector", &ElectronSector, "ElectronSector/I");
    if( topology_has_protons ) {
      kAnalysisTree -> Branch( "proton_mom", &proton_mom, "proton_mom/D");
      kAnalysisTree -> Branch( "proton_momx", &proton_momx, "proton_momx/D");
      kAnalysisTree -> Branch( "proton_momy", &proton_momy, "proton_momy/D");
      kAnalysisTree -> Branch( "proton_momz", &proton_momz, "proton_momz/D");
      kAnalysisTree -> Branch( "proton_theta", &proton_theta, "proton_theta/D");
      kAnalysisTree -> Branch( "proton_phi", &proton_phi, "proton_phi/D");
      kAnalysisTree -> Branch( "ECal", &ECal, "ECal/D");
      kAnalysisTree -> Branch( "AlphaT", &AlphaT, "AlphaT/D");
      kAnalysisTree -> Branch( "DeltaPT", &DeltaPT, "DeltaPT/D");
      kAnalysisTree -> Branch( "DeltaPhiT", &DeltaPhiT, "DeltaPhiT/D");
    }

    if( topology_has_pip ) {
      kAnalysisTree -> Branch( "pip_mom", &pip_mom, "pip_mom/D");
      kAnalysisTree -> Branch( "pip_momx", &pip_momx, "pip_momx/D");
      kAnalysisTree -> Branch( "pip_momy", &pip_momy, "pip_momy/D");
      kAnalysisTree -> Branch( "pip_momz", &pip_momz, "pip_momz/D");
      kAnalysisTree -> Branch( "pip_theta", &pip_theta, "pip_theta/D");
      kAnalysisTree -> Branch( "pip_phi", &pip_phi, "pip_phi/D");
    }

    if( topology_has_pim ) {
      kAnalysisTree -> Branch( "pim_mom", &pim_mom, "pim_mom/D");
      kAnalysisTree -> Branch( "pim_momx", &pim_momx, "pim_momx/D");
      kAnalysisTree -> Branch( "pim_momy", &pim_momy, "pim_momy/D");
      kAnalysisTree -> Branch( "pim_momz", &pim_momz, "pim_momz/D");
      kAnalysisTree -> Branch( "pim_theta", &pim_theta, "pim_theta/D");
      kAnalysisTree -> Branch( "pim_phi", &pim_phi, "pim_phi/D");
    }

    n = false ; 
  }
  
  kAnalysisTree -> Fill();
  return true ; 
}
