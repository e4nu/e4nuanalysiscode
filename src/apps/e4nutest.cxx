// _____________________________________________________________
/* This app is used to test aspects of the e4nu software only */
// _____________________________________________________________

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TH1D.h"
#include "../physics/MCEventHolder.h"
#include "../conf/ConfigurablesI.h"

using namespace std; 
using namespace e4nu;

int main( void ) {
  std::cout << "Test ongoing..." << std::endl;
  std::string file_name = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/apapadop_SuSav2_C12_2261GeV_master.root";

  ConfigurablesI conf() ; 
  MCEventHolder * mc_events = new MCEventHolder(file_name.c_str(), 10);
   //EventHolderI * events = new EventHolderI( file_name.c_str() );

  return 0 ; 
}
