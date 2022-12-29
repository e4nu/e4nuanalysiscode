// _____________________________________________________________
/* This app is used to test aspects of the e4nu software only */
// _____________________________________________________________

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TH1D.h"
#include "analysis/E4NuAnalysis.h"

using namespace std; 
using namespace e4nu;

int main( void ) {
  std::cout << "Test ongoing..." << std::endl;
  std::string file_name = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_1161GeV/apapadop_SuSav2_C12_1161GeV_master.root";

  E4NuAnalysis * analysis = new E4NuAnalysis("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vfork/ConfFiles/example_configuration.txt") ;
  if( ! analysis -> LoadData( file_name.c_str() ) ) return 0 ;  

  analysis -> Analyse() ; 
  analysis -> SubstractBackground() ; 
  analysis -> Finalise("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vfork/output.txt") ; 

  delete analysis ;
  return 0 ; 
}
