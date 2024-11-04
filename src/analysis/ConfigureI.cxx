/** 
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/ConfigureI.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/ParticleUtils.h"
#include "utils/DetectorUtils.h"
#include "utils/Utils.h"

using namespace e4nu; 

ConfigureI::ConfigureI( ) {
  this->Initialize();
}

ConfigureI::ConfigureI( const double EBeam, const unsigned int TargetPdg ) { 
  kEBeam = EBeam ; 
  kTargetPdg = TargetPdg ;
  kIsDataLoaded = false ;
  kIsConfigured = true ; 

  this->Initialize();
}
    
ConfigureI::~ConfigureI() {
  for ( auto it = kAcceptanceMap.begin(); it != kAcceptanceMap.end(); it++){
    delete it->second ;
  }

  for ( auto it = kAccMap.begin(); it != kAccMap.end(); it++){
    delete it->second ;
  }

  for ( auto it = kGenMap.begin(); it != kGenMap.end(); it++){
    delete it->second ;
  }

  kAcceptanceMap.clear();
  kAccMap.clear();
  kGenMap.clear();

  kTopology_map.clear() ; 
  kObservables.clear();
  if( kFiducialCut ) delete kFiducialCut ;
}

ConfigureI::ConfigureI( const std::string input_file ) {
  std::cout << " Configuring analysis from " << input_file << " ..." << std::endl;
  std::ifstream file (input_file.c_str()) ;
      
  std::vector<std::string> param, value ; 
  std::string parameter_name, line, parameter_value ; 
  if( file.is_open() ) { 
    while ( getline (file,line) ) {
      if(line.find("#", 0) != std::string::npos ) continue;
      if(line == "" ) continue;
      std::istringstream line_as_stream( line ) ; 
      line_as_stream >> parameter_name >> parameter_value ; 
      param.push_back(parameter_name) ; 
      value.push_back(parameter_value) ; 
    }
  }
      
  for ( unsigned int i = 0 ; i < param.size(); ++i ) {
    if ( param[i] == "UseAllSectors" ) {
      if( value[i] == "true" ) kUseAllSectors = true ; 
      else kUseAllSectors = false ; 
    } else if ( param[i] == "DisableSector" ) { 
      std::istringstream iss(value[i]);
      std::string token ; 
      while (std::getline(iss, token, ',')) {
        unsigned int sector = stoi(token);
        kEnabledSectors[sector] = false ; 
      }
    } else if ( param[i] == "ApplyFiducial" ){
      if( value[i] == "false" ) kApplyFiducial = false ; 
      else {
	kApplyFiducial = true ; 
      }
    } else if ( param[i] == "ApplyAccWeights" ){
      if ( value[i] == "false" ) kApplyAccWeights = false ; 
      else kApplyAccWeights = true ; 
    } else if ( param[i] == "ApplyCorrWeights" ){
      if ( value[i] == "false" ) kApplyCorrWeights = false ; 
      else kApplyCorrWeights = true ; 
    } else if ( param[i] == "ApplyMottWeight" ){
      if ( value[i] == "false" ) kApplyMottWeight = false ; 
      else kApplyMottWeight = true ; 
    } else if ( param[i] == "ApplyReso" ) {
      if ( value[i] == "false" ) kApplyReso = false ;
      else kApplyReso = true ; 
    } else if ( param[i] == "ApplyMomCut" ) {
      if ( value[i] == "false" ) kApplyMomCut = false ;
      else kApplyMomCut = true ; 
    } else if ( param[i] == "ApplyPhiOpeningAngle" ) {
      if( value[i] == "true" ) kApplyPhiOpeningAngle = true ; 
      else kApplyPhiOpeningAngle = false ; 
    } else if ( param[i] == "UsePhiThetaBand" ) {
      if( value[i] == "true" ) kUsePhiThetaBand = true ;
      else kUsePhiThetaBand = false ; 
    } else if ( param[i] == "ApplyThetaSlice" ) {
      if ( value[i] == "true" ) kApplyThetaSlice = true ; 
      else kApplyThetaSlice = false ; 
    } else if ( param[i] == "ApplyGoodSectorPhiSlice" ) { 
      if ( value[i] == "true" ) kApplyGoodSectorPhiSlice = true ; 
      else kApplyGoodSectorPhiSlice = false ;
    } else if ( param[i] == "IsData" ) {
      if( value[i] == "true" ) kIsData = true ; 
      else kIsData = false ; 
    } else if ( param[i] == "IsCLAS6Analysis") { 
      if( value[i] == "true" ) { kIsCLAS6Analysis = true ; kIsCLAS12Analysis = false ;} 
      else kIsCLAS6Analysis = false ; 
    } else if ( param[i] == "IsCLAS12Analysis") { 
      if( value[i] == "true" ) { kIsCLAS12Analysis = true ; kIsCLAS6Analysis = false; } 
      else kIsCLAS12Analysis = false ; 
    } else if ( param[i] == "IsRadiated" ) { 
      if( value[i] == "true" ) { kIsRadiated = true ; }
      else kIsRadiated = false ; 
    } else if( param[i] == "ApplyQ2Cut" ) { 
      if( value[i] == "true" ) kQ2Cut = true ; 
      else kQ2Cut = false ; 
    } else if( param[i] == "ApplyWCut" ) { 
      if( value[i] == "true" ) kWCut = true ; 
      else kWCut = false ; 	
    } else if ( param[i] == "ApplyOutEMomCut" ) {
      if( value[i] == "true" ) kOutEMomCut = true ; 
      else kOutEMomCut = false ; 
    } else if( param[i] == "IsElectronData" ) { 
      if( value[i] == "true" ) kIsElectron = true ; 
      else kIsElectron = false ; 
    } else if ( param[i] == "AnalysisTypeID" ) { 
      kAnalysisTypeID = std::stoi(value[i]) ; 
    } else if ( param[i] == "offset" ) koffset = std::stod( value[i] ) ; 
    else if ( param[i] == "NoFSI") {
      if( value[i] == "true" ) kNoFSI = true ; 
      else kNoFSI = false ; 
    } else if ( param[i] == "EBeam" ) kEBeam = std::stod( value[i] ) ; 
    else if ( param[i] == "TargetPdg" ) kTargetPdg = (unsigned int) std::stoi( value[i] ) ; 
    else if ( param[i] == "NEvents" ) kNEvents = (unsigned int) std::stoi( value[i] ) ;
    else if ( param[i] == "FirstEvent" ) kFirstEvent = (unsigned int) std::stoi( value[i] ) ;
    else if ( param[i] == "Toplogy") {
      std::string element, m_element ;
      std::istringstream particle_list( value[i] ) ;
      while( getline( particle_list, element, ',' ) ) {
	std::istringstream part_list_mult ( element ) ;
	std::vector<std::string> temp ; 
	while ( getline( part_list_mult, m_element, ':' ) ) {
	  temp.push_back( m_element ) ; 
	}
	if( temp.size() == 2 ) {
	  std::pair<int, unsigned int> pair ( std::stoi( temp[0]), (unsigned int) std::stoi( temp[1] ) ) ; // Pdg, multiplicity
	  kTopology_map.insert( pair ) ;
	}
      }
    } else if ( param[i] == "TrueSignal" ) { 
      if ( value[i] == "true" ) kTrueSignal = true ; 
      else kTrueSignal = false ; 
    } else if ( param[i] == "SubtractBkg" ) { 
      if( value[i] == "true" ) kSubtractBkg = true ; 
      else kSubtractBkg = false ; 
    } else if ( param[i] == "MaxBackgroundMultiplicity" ) { kMaxBkgMult = (unsigned int) std::stoi( value[i] ) ;
    } else if ( param[i] == "NRotations" ) { kNRotations = (unsigned int) std::stoi( value[i] ) ;
    } else if ( param[i] == "AnalysisKey" ) { kAnalysisKey = value[i] ;
    } else if ( param[i] == "ObservableList" ) {
      std::string obs ; 
      std::istringstream obs_list( value[i] ) ; 
      while( getline( obs_list, obs, ',' ) ) { 
	kObservables.push_back(obs) ; 
      }
    } else if ( param[i] == "NormalizeHists" ) {
      if ( value[i] == "true" ) kNormalize = true ; 
      else kNormalize = false ; 
      if( kNormalize && ! kApplyCorrWeights ) {
	std::cout << " WARNING: You are trying to Normalize the cross section without correction weights. Aborting..." << std::endl;
	kIsConfigured = false ; 
	break ;
      } 
    } else if ( param[i] == "ComputeTrueAccCorr" ) { 
      if( value[i] == "true" ) { kComputeTrueAccCorr = true ; }
      else kComputeTrueAccCorr = false ; 
    } else if ( param[i] == "ComputeTrueRecoAccCorr" ) { 
      if( value[i] == "true" ) { kComputeTrueRecoAccCorr = true ; }
      else kComputeTrueRecoAccCorr = false ; 
    } else if ( param[i] == "DebugBkg") {
      if( value[i] == "false" ) {
	kDebugBkg = false ; 
      } else { kDebugBkg = true ; }  
    } else if ( param[i] == "OutputFile" ) {
      kOutputFile = value[i] ;
    } else if ( param[i] == "InputFile" ) {
      kInputFile = value[i] ;
    } else if ( param[i] == "XSecFile" ) {
      kXSecFile = value[i] ;
    } else if ( param[i] == "SystFidAngleShift" ) { 
      fFidAngleShift = stod(value[i]) ;
    }
  }

  if( kInputFile == "" ) {
    std::cout << " ERROR : Input file not specified " << std::endl;
    kIsConfigured = false ; 
  } else { std::cout << " Input file: " << kInputFile << std::endl;}

  if( kOutputFile == "" ) {
    std::cout << " ERROR : Output file not specified " << std::endl;
    kIsConfigured = false ; 
  }
  if( kXSecFile == "" && !IsData() ) {
    std::cout << " ERROR : XSec file not specified " << std::endl;
    kIsConfigured = false ; 
  }
  
  if( !kIsCLAS6Analysis && !kIsCLAS6Analysis ) {
    std::cout << " WARN : Analysis type not configured. Using CLAS6... " << std::endl;
    kIsCLAS6Analysis = true ;
  }

  if( kIsRadiated ) std::cout << " Using MC file with radiative effects " << std::endl;

  if( !kIsCLAS6Analysis && !kIsCLAS6Analysis ) {
    std::cout << " ERROR : CLAS12 analysis not available yet..."<<std::endl;
    kIsConfigured = false ;
  }
    
  this->Initialize();
}

void ConfigureI::Initialize(void){

  gRandom = new TRandom3() ; 
  gRandom->SetSeed(10);

  if( !kIsConfigured ) std::cout << " CONFIGURATION FAILED..." << std::endl;
  if( !IsData() ) kAnalysisTree = std::unique_ptr<TTree>( new TTree("MCCLAS6Tree","GENIE CLAS6 Tree") ) ; 
  else kAnalysisTree = std::unique_ptr<TTree>( new TTree("CLAS6Tree","CLAS6 Tree") ) ; 

  if( ApplyFiducial() &&  kIsConfigured ) kIsConfigured = InitializeFiducial() ; 
  
  // Initialize acceptance map histograms from file
  if( ApplyAccWeights() ) {
    double EBeam = GetConfiguredEBeam() ; 
    unsigned int Target = GetConfiguredTarget() ;
    kAcceptanceMap = conf::GetAcceptanceFileMap2( Target, EBeam ) ;

    kAccMap[conf::kPdgElectron] = dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgElectron] -> Get("Accepted Particles") ) ;
    kAccMap[conf::kPdgProton] = dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgProton] -> Get("Accepted Particles") );
    kAccMap[conf::kPdgPiP] = dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiP] -> Get("Accepted Particles") );
    kAccMap[conf::kPdgPiM] = dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiM] -> Get("Accepted Particles") );

    kGenMap[conf::kPdgElectron] = dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgElectron] -> Get("Generated Particles") );
    kGenMap[conf::kPdgProton] = dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgProton] -> Get("Generated Particles") ) ;
    kGenMap[conf::kPdgPiP] = dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiP] -> Get("Generated Particles") ) ;
    kGenMap[conf::kPdgPiM] = dynamic_cast<TH3D*>( kAcceptanceMap[conf::kPdgPiM] -> Get("Generated Particles") ) ;

    for( auto it = kAccMap.begin() ; it != kAccMap.end(); ++it ) {
      kAccMap[it->first]->SetDirectory(nullptr);
      kGenMap[it->first]->SetDirectory(nullptr);
    }
  }
}
      
void ConfigureI::PrintConfiguration(void) const { 
  std::cout << "*********************************************************************" << std::endl;
  std::cout << "*                         E4NU ANALYSIS CONF                       **" << std::endl;
  std::cout << "*********************************************************************" << std::endl;
  if( !UseAllSectors() ) std::cout << " UseAllSectors: " << UseAllSectors() << std::endl;
  else { 
    std::cout << "Enabled Sectors: " << std::endl;
    for( unsigned int i = 0 ; i < kEnabledSectors.size() ; ++i ) { 
      std::cout << "Sector " << i << " : " << kEnabledSectors[i] << std::endl; 
    }
  }
  std::cout << "ApplyOutMomCut: " <<kOutEMomCut << std::endl;
  if( kQ2Cut ) std::cout << "ApplyQ2Cut: "<<kQ2Cut<< std::endl;
  if( kWCut ) std::cout << "ApplyWCut: "<<kWCut<< std::endl;
  if( kApplyFiducial ) std::cout << "ApplyFiducial:" << kApplyFiducial << std::endl;
  if( kApplyAccWeights ) std::cout << "ApplyAccWeights: " << kApplyAccWeights << std::endl;
  if( kApplyMottWeight )std::cout << "ApplyMottWeight: " << kApplyMottWeight << std::endl;
  if( kApplyReso ) std::cout << "ApplyReso:" << kApplyReso << std::endl;
  if( kApplyMomCut ) std::cout << "ApplyMomCut:" << kApplyMomCut << std::endl;
  if( kApplyPhiOpeningAngle) std::cout << "ApplyPhiOpeningAngle:" << kApplyPhiOpeningAngle << std::endl;
  if( kUsePhiThetaBand ) std::cout << "UsePhiThetaBand:"<< kUsePhiThetaBand << std::endl;
  if( kApplyThetaSlice ) std::cout << "ApplyThetaSlice:"<< kApplyThetaSlice << std::endl;
  if( kApplyGoodSectorPhiSlice ) std::cout << "ApplyGoodSectorPhiSlice:"<<kApplyGoodSectorPhiSlice << "\n"<<std::endl;
  if( fFidAngleShift!=0 ) std::cout << "Phi shifted by " << fFidAngleShift << std::endl;
  std::cout << "EBeam: " << kEBeam << " GeV " << std::endl;
  std::cout << "Target Pdg : " << kTargetPdg << "\n" <<std::endl;
  if ( kIsData ) std::cout << "\nIsData" << std::endl;
  else {
    std::cout << "Is MC Data " << std::endl;
    if ( kNoFSI ) std::cout << " No FSI " << std::endl ; 
  }
  if( kIsCLAS6Analysis ) std::cout << " Analysing CLAS6 "<<std::endl;
  if( kIsCLAS12Analysis ) std::cout << " Analysing CLAS12 "<<std::endl;
  if( kIsElectron ) std::cout << " Electron scattering \n" <<std::endl;
  else std::cout << " Neutrino scattering \n" <<std::endl;
  std::cout << " Analysis Type ID: " << GetAnalysisTypeID() << std::endl;
  std::cout << "Topology: " << std::endl;
  for ( auto it = kTopology_map.begin() ; it != kTopology_map.end() ; ++it ) { 
    std::cout << "    " << utils::PdgToString(it->first) << ", multiplicity " << it->second << std::endl;
  }

  if( kSubtractBkg ) {
    std::cout << "\nBackground Subtraction enabled : " << std::endl;
    std::cout << "Maximum Background Multiplicity: "<< kMaxBkgMult << std::endl;
    std::cout << "Number of rotations: "<< kNRotations << "\n" << std::endl;
  }

  if( kTrueSignal ) std::cout << "Only analysing true signal events : " << std::endl;

  for( unsigned int i = 0 ; i < kObservables.size(); ++i ) {
    std::cout << "Observable " << kObservables[i] << std::endl;
  }
  if( kDebugBkg ) std::cout << " Storing debugging plots for background " << std::endl;

  std::cout << "\nXSecFile " << kXSecFile << std::endl;
  std::cout << "\nStoring output in " << kOutputFile << std::endl;
  std::cout << "Analizing " << kNEvents << " ... " <<std::endl;
  if( kFirstEvent != 0 ) std::cout << " startint from event " << kFirstEvent << std::endl;
  std::cout << "*********************************************************************" << std::endl;

}

unsigned int ConfigureI::GetNTopologyParticles(void) {
  unsigned int N_signal = 0 ;
  for( auto it = kTopology_map.begin() ; it != kTopology_map.end() ; ++it ) {
    if( it->first == conf::kPdgElectron ) continue ; 
    N_signal += it->second ; 
  }
  return N_signal ;
}


bool ConfigureI::InitializeFiducial(void) {
  double EBeam = GetConfiguredEBeam() ; 
  unsigned int Target = GetConfiguredTarget() ;
  
  // Initialize fiducial for this run
  kFiducialCut = new Fiducial() ; 
  kFiducialCut -> InitPiMinusFit( EBeam ) ; 
  kFiducialCut -> InitEClimits(); 
  kFiducialCut -> up_lim1_ec -> Eval(60) ;
  kFiducialCut -> SetConstants( conf::GetTorusCurrent( EBeam ), Target , EBeam ) ;

  kFiducialCut -> SetFiducialCutParameters( EBeam ) ;
  if( !kFiducialCut ) return false ; 

  return true ; 
}


void ConfigureI::ApplyAcceptanceCorrection( Event & event, bool invert ) {
  double acc_wght = 1 ;
  if( ApplyAccWeights() || kIsData ) {
    // We only apply the acceptance if configured
    // For the data, acceptance is used to correct background events for this effect
    TLorentzVector out_mom = event.GetOutLepton4Mom() ;
    std::map<int,std::vector<TLorentzVector>> part_map = event.GetFinalParticles4Mom() ;
    std::map<int,unsigned int> Topology = GetTopology();
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
    double initial_accwght = event.GetAccWght(); 
    
    if( invert && acc_wght != 0 ) acc_wght = initial_accwght / acc_wght ; 
    event.SetAccWght(acc_wght);
  }
  return ;
}


