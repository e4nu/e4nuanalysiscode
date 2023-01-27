/** 
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/AnalysisI.h"

using namespace e4nu; 

AnalysisI::AnalysisI( ) {;}

AnalysisI::AnalysisI( const double EBeam, const unsigned int TargetPdg ) { 
  kEBeam = EBeam ; 
  kTargetPdg = TargetPdg ;
  kIsDataLoaded = false ;
  kIsConfigured = true ; 
  InitializeFiducial();
  PrintConfiguration();
}
    
AnalysisI::~AnalysisI() {
  kTopology_map.clear() ; 
  kObservables.clear();
  kNBins.clear();
  kRanges.clear();
  delete kFiducialCut ; 
}


AnalysisI::AnalysisI( const std::string input_file ) {
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
      else kApplyFiducial = true ; 
    } else if ( param[i] == "ApplyAccWeights" ){
      if ( value[i] == "false" ) kApplyAccWeights = false ; 
      else kApplyAccWeights = true ; 
    } else if ( param[i] == "ApplyReso" ) {
      if ( value[i] == "false" ) kApplyReso = false ;
      else kApplyReso = true ; 
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
    } else if( param[i] == "ApplyQ2Cut" ) { 
      if( value[i] == "true" ) kQ2Cut = true ; 
      else kQ2Cut = false ; 
    } else if( param[i] == "ApplyWCut" ) { 
      if( value[i] == "true" ) kWCut = true ; 
      else kWCut = false ; 	
    } else if ( param[i] == "ApplyOutMomCut" ) {
      if( value[i] == "true" ) fOutMomCut = true ; 
      else fOutMomCut = false ; 
    }else if( param[i] == "IsElectronData" ) { 
      if( value[i] == "true" ) kIsElectron = true ; 
      else kIsElectron = false ; 
    }else if ( param[i] == "offset" ) koffset = std::stod( value[i] ) ; 
    else if ( param[i] == "EBeam" ) kEBeam = std::stod( value[i] ) ; 
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
      
  if( kIsConfigured ) kIsConfigured = InitializeFiducial() ; 
  if( kIsConfigured ) PrintConfiguration() ;
}

bool AnalysisI::InitializeFiducial(void) {
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
  } else { return true ; }
  if( !kFiducialCut ) return false ; 
  return true ; 
}
      
void AnalysisI::PrintConfiguration(void) const { 
  std::cout << "*********************************************************************" << std::endl;
  std::cout << "*                         E4NU ANALYSIS CONF                       **" << std::endl;
  std::cout << "*********************************************************************" << std::endl;
  std::cout << "UseAllSectors: " << kUseAllSectors << std::endl;
  std::cout << "ApplyOutMomCut: " <<fOutMomCut << std::endl;
  std::cout << "ApplyQ2Cut: "<<kQ2Cut<< std::endl;
  std::cout << "ApplyWCut: "<<kWCut<< std::endl;
  std::cout << "ApplyFiducial:" << kApplyFiducial << std::endl;
  std::cout << "ApplyAccWeights: " << kApplyAccWeights << std::endl;
  std::cout << "ApplyReso:" << kApplyReso << std::endl;
  std::cout << "ApplyQ2Cut: " << kQ2Cut << std::endl;
  std::cout << "ApplyWCut: " << kWCut << std::endl;
  std::cout << "ApplyPhiOpeningAngle:" << kApplyPhiOpeningAngle << std::endl;
  std::cout << "UsePhiThetaBand:"<< kUsePhiThetaBand << std::endl;
  std::cout << "ApplyThetaSlice:"<< kApplyThetaSlice << std::endl;
  std::cout << "ApplyGoodSectorPhiSlice:"<<kApplyGoodSectorPhiSlice << "\n"<<std::endl;
  std::cout << "EBeam: " << kEBeam << " GeV " << std::endl;
  std::cout << "Target Pdg : " << kTargetPdg << "\n" <<std::endl;
  if ( kIsData ) std::cout << "\nIsData" << std::endl;
  else std::cout << "Is MC Data " << std::endl;
  if( kIsElectron) std::cout << " Electron scattering data \n" <<std::endl;
  std::cout << "Topology: " << std::endl;
  for ( auto it = kTopology_map.begin() ; it != kTopology_map.end() ; ++it ) { 
    std::cout << "    " << it->first << ", multiplicity " << it->second << std::endl;
  }
  std::cout << "Maximum Background Multiplicity: "<< kMaxBkgMult << "\n" << std::endl;
  std::cout << "Number of rotations: "<< kNRotations << "\n" << std::endl;
  for( unsigned int i = 0 ; i < kObservables.size(); ++i ) {
    std::cout << "Observable " << kObservables[i] << std::endl;
    std::cout << "Number of bins = " << kNBins[i] << std::endl;
    std::cout << "Range = {"<<kRanges[i][0]<<","<<kRanges[i][1]<<"}\n"<<std::endl;
  }
  std::cout << "\nXSecFile " << kXSecFile << std::endl;
  std::cout << "\nStoring output in " << kOutputFile << std::endl;
  std::cout << "Analizing " << kNEvents << " ... " <<std::endl;
  if( kFirstEvent != 0 ) std::cout << " startint from event " << kFirstEvent << std::endl;
  std::cout << "*********************************************************************" << std::endl;
}

unsigned int AnalysisI::GetNTopologyParticles(void) {
  unsigned int N_signal = 0 ;
  for( auto it = kTopology_map.begin() ; it != kTopology_map.end() ; ++it ) {
    if( it->first == conf::kPdgElectron ) continue ; 
    if ( it -> second != 0 ) ++N_signal ; 
  }
  return N_signal ;
}
