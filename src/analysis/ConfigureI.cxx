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
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/ParticleUtils.h"
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
  kTopology_map.clear() ; 
  kObservables.clear();
  kNBins.clear();
  kRanges.clear();
  delete kFiducialCut ;
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
    } else if ( param[i] == "ApplyFiducial" ){
      if( value[i] == "false" ) kApplyFiducial = false ; 
      else {
	kApplyFiducial = true ; 
	kApplyHadFiducial = true ; 
	kApplyEFiducial = true ; 
      }
    } else if ( param[i] == "ApplyHadFiducial" ){
      if( value[i] == "false" ) kApplyHadFiducial = false ; 
      else kApplyHadFiducial = true ; 
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
    } else if ( param[i] == "SubtractBkg" ) { 
      if( value[i] == "true" ) kSubtractBkg = true ; 
      else kSubtractBkg = false ; 
    } else if ( param[i] == "MaxBackgroundMultiplicity" ) { kMaxBkgMult = (unsigned int) std::stoi( value[i] ) ;
    } else if ( param[i] == "NRotations" ) { kNRotations = (unsigned int) std::stoi( value[i] ) ;
    } else if ( param[i] == "ObservableList" ) {
      std::string obs ; 
      std::istringstream obs_list( value[i] ) ; 
      while( getline( obs_list, obs, ',' ) ) { 
	kObservables.push_back(obs) ; 
      }
    } else if ( param[i] == "NBinsList" ) {
      std::string nb ;
      std::istringstream nbs_list( value[i] ) ;
      while( getline( nbs_list, nb, ',' ) ) {
	kNBins.push_back((unsigned int) std::stoi(nb)) ;
      }
    } else if ( param[i] == "RangeList" ) {
      std::istringstream ranges_list( value[i] ) ;
      std::string m_element,range ;
      while ( getline( ranges_list, range, ',' ) ) {
	std::istringstream ele(range) ; 
	std::vector<double> myrange ; 
	while ( getline( ele, m_element, ':' ) ) {
	  myrange.push_back(std::stod(m_element)); 
	}
	if( myrange.size() != 2 ) {
	  std::cout<< " Range size is not 2 !! " << std::endl;
	  kIsConfigured = false ; 
	  break ; 
	}
	kRanges.push_back( myrange ) ; 
      }
    } else if ( param[i] == "NormalizeHists" ) {
      if ( value[i] == "true" ) kNormalize = true ; 
      else kNormalize = false ; 
      if( kNormalize && ! kApplyCorrWeights ) {
	std::cout << " WARNING: You are trying to Normalize the cross section without correction weights. Aborting..." << std::endl;
	kIsConfigured = false ; 
	break ;
      } 
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
    }
  }

  if( kObservables.size() != kRanges.size() || kObservables.size() != kNBins.size() || kRanges.size()!= kNBins.size() ){
    std::cout << " ERROR : Ranges don't match !! " << std::endl;
    kIsConfigured = false ; 
  }

  if( kInputFile == "" ) {
    std::cout << " ERROR : Input file not specified " << std::endl;
    kIsConfigured = false ; 
  }
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

  if( !kIsCLAS6Analysis && !kIsCLAS6Analysis ) {
    std::cout << " ERROR : CLAS12 analysis not available yet..."<<std::endl;
    kIsConfigured = false ;
  }

  this->Initialize();
}

void ConfigureI::Initialize(void){

  gRandom = new TRandom3() ; 
  gRandom->SetSeed(10);

  if( ApplyFiducial() &&  kIsConfigured ) kIsConfigured = InitializeFiducial() ; 

  if( kIsConfigured ) PrintConfiguration() ;
  else std::cout << " CONFIGURATION FAILED..." << std::endl;

}
      
void ConfigureI::PrintConfiguration(void) const { 
  std::cout << "*********************************************************************" << std::endl;
  std::cout << "*                         E4NU ANALYSIS CONF                       **" << std::endl;
  std::cout << "*********************************************************************" << std::endl;
  std::cout << "UseAllSectors: " << kUseAllSectors << std::endl;
  std::cout << "ApplyOutMomCut: " <<kOutEMomCut << std::endl;
  std::cout << "ApplyQ2Cut: "<<kQ2Cut<< std::endl;
  std::cout << "ApplyWCut: "<<kWCut<< std::endl;
  std::cout << "ApplyFiducial:" << kApplyFiducial << std::endl;
  std::cout << "ApplyEFiducial:" << kApplyEFiducial << std::endl;
  std::cout << "ApplyHadFiducial:" << kApplyHadFiducial << std::endl;
  std::cout << "ApplyAccWeights: " << kApplyAccWeights << std::endl;
  std::cout << "ApplyMottWeight: " << kApplyMottWeight << std::endl;
  std::cout << "ApplyReso:" << kApplyReso << std::endl;
  std::cout << "ApplyMomCut:" << kApplyMomCut << std::endl;
  std::cout << "ApplyQ2Cut: " << kQ2Cut << std::endl;
  std::cout << "ApplyWCut: " << kWCut << std::endl;
  std::cout << "ApplyPhiOpeningAngle:" << kApplyPhiOpeningAngle << std::endl;
  std::cout << "UsePhiThetaBand:"<< kUsePhiThetaBand << std::endl;
  std::cout << "ApplyThetaSlice:"<< kApplyThetaSlice << std::endl;
  std::cout << "ApplyGoodSectorPhiSlice:"<<kApplyGoodSectorPhiSlice << "\n"<<std::endl;
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
  
  std::cout << "Topology: " << std::endl;
  for ( auto it = kTopology_map.begin() ; it != kTopology_map.end() ; ++it ) { 
    std::cout << "    " << utils::PdgToString(it->first) << ", multiplicity " << it->second << std::endl;
  }

  if( kSubtractBkg ) {
    std::cout << "\nBackground Subtraction enabled : " << std::endl;
    std::cout << "Maximum Background Multiplicity: "<< kMaxBkgMult << std::endl;
    std::cout << "Number of rotations: "<< kNRotations << "\n" << std::endl;
  }

  for( unsigned int i = 0 ; i < kObservables.size(); ++i ) {
    std::cout << "Observable " << kObservables[i] << std::endl;
    std::cout << "Number of bins = " << kNBins[i] << std::endl;
    std::cout << "Range = {"<<kRanges[i][0]<<","<<kRanges[i][1]<<"}\n"<<std::endl;
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

  if( ApplyFiducial() ) {
    // Initialize fiducial for this run
    kFiducialCut = new Fiducial() ; 
    kFiducialCut -> InitPiMinusFit( EBeam ) ; 
    kFiducialCut -> InitEClimits(); 
    kFiducialCut -> up_lim1_ec -> Eval(60) ;
    kFiducialCut -> SetConstants( conf::GetTorusCurrent( EBeam ), Target , EBeam ) ;
    kFiducialCut -> SetFiducialCutParameters( EBeam ) ;
    if( !kFiducialCut ) return false ; 
  } else { return true ; }

  return true ; 
}


