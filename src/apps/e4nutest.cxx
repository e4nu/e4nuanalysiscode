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
  std::string file_name = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/apapadop_SuSav2_C12_2261GeV_master.root";

  E4NuAnalysis * analysis = new E4NuAnalysis("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster/data/ConfFiles/example_configuration.txt") ;
  analysis -> LoadData( file_name.c_str(), 10 ) ; 
  for ( unsigned i = 0 ; i < 10 ; ++i ) {
    analysis -> GetEvent(i) ;
  }

  return 0 ; 
}
