/*
 * This class reads the run information of CLAS6 e2a data
 *
 */

#include <fstream> 
#include <sstream>
#include "analysis/CLAS6RunInfo.h"

using namespace e4nu;

CLAS6RunInfo::CLAS6RunInfo()
{
  this->LoadData();
}

CLAS6RunInfo::~CLAS6RunInfo()
{

}

void CLAS6RunInfo::LoadData(void)
{
  static const  std::string local_path = std::getenv("E4NUANALYSIS");
  static const std::string run_info = local_path + "/data/Larry_e2GoodRunList.csv";
  std::cout << " Opening file " << run_info << std::endl;
  ifstream inputFile ;
  inputFile.open(run_info.c_str());
  if (inputFile.is_open()){
    std::string line;
    
    while (std::getline(inputFile, line)) {
      std::stringstream ss(line);
      std::string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15;

      std::getline(ss, col1, ';');
      if( atoi(col1.c_str()) == 0 ) continue ; // Skip line if it does not contain run number
      else kRunID = atoi(col1.c_str()) ; 
      
      std::getline(ss, col2, ';');
      std::getline(ss, col3, ';');
      std::getline(ss, col4, ';');
      std::getline(ss, col5, ';');
      std::getline(ss, col6, ';');
      std::getline(ss, col7, ';');
      std::getline(ss, col8, ';');
      std::getline(ss, col9, ';');
      std::getline(ss, col10, ';');
      std::getline(ss, col11, ';');
      std::getline(ss, col12, ';');
      std::getline(ss, col13, ';');
      std::getline(ss, col14, ';');
      std::getline(ss, col15);

      kDate[kRunID] = col2 ; 
      kTime[kRunID] = col3 ; 
      if( atoi(col4.c_str()) != 0 ) kNFiles[kRunID] = atoi(col4.c_str()) ;
      if( atoi(col5.c_str()) != 0 ) kNRunEvents[kRunID] = atoi(col5.c_str()) ;
      kTarget[kRunID] = col6 ;
      if( atof(col7.c_str()) != 0 ) kLiveC[kRunID] = atof(col7.c_str());
      if( atof(col8.c_str()) !=0 ) kBeam[kRunID] = atof(col8.c_str())/1000.;
      if( atof(col9.c_str()) !=0 ) kI[kRunID] = atof(col9.c_str());
      if( atoi(col10.c_str()) != 0 ) kTorus[kRunID] = atoi(col10.c_str());
      if( atoi(col11.c_str()) != 0 ) kMini[kRunID] = atoi(col11.c_str());
      if( atoi(col12.c_str()) != 0 ) kechi[kRunID] = atoi(col12.c_str());
      if( atof(col13.c_str()) != 0 ) kFcupA[kRunID] = atof(col13.c_str())*0.0000001;
      kComments[kRunID] = col15;

      // By default enable all 
      kEnabledRun[kRunID] = true ; 
      
      // In case we missed one, we will check on run time if the run is used or not. If it is not used it will stay in false 
      kIsRunInFile[kRunID] = false ;
    }
    //while (inputFile.good()) cout << (char) inputFile.get();
    inputFile.close();
  }
  else {
    cout << "Error opening file";
  }
  
}

void CLAS6RunInfo::RetrieveAllRuns( const double beam, const string target ) {

  for ( const auto & runset : kEnabledRun ) {
    if( kBeam[runset.first] != beam ) kEnabledRun[runset.first] = false;
    if( kTarget[runset.first] != target ) kEnabledRun[runset.first] = false;
  }

}

void CLAS6RunInfo::NeglectRuns( const std::vector<unsigned int> run_ids ) {
  // Remove runs if asked. 
  for( unsigned int i = 0 ; i < run_ids.size() ; ++i ) { 
    kEnabledRun[run_ids[i]] = false;
  }
}

double CLAS6RunInfo::GetTotalIntegratedCharge(void) {
  double total_ch = 0 ;
  for ( const auto & runset : kEnabledRun ) {
    if( runset.second == true && kIsRunInFile[runset.first] == true ) { 
      total_ch += GetRunIntegratedCharge(runset.first) ; 
    }
  }
  return total_ch ; 
} 

bool CLAS6RunInfo::IsRunValid( const unsigned int run_id) {
  bool valid = kEnabledRun[run_id];
  if( !valid ) {
    //std::cout << " Trying to use " << run_id << " run. Target: " << kTarget[run_id] << " Beam: " << kBeam[run_id] << " GeV " << std::endl;
  }
  return valid;
}
void CLAS6RunInfo::SetIsRunInFile( const unsigned int run_id ) { 
  kIsRunInFile[run_id] = true; 
  if( kEnabledRun[run_id] == false ) { 
    std::cout << " WARNING: Run in File " << run_id << " run. Target: " << kTarget[run_id] << " Beam: " << kBeam[run_id] << " GeV " << std::endl;
  }
}
