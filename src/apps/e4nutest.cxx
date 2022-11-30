// _____________________________________________________________
/* This app is used to test aspects of the e4nu software only */
// _____________________________________________________________

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TH1D.h"
#include "../physics/EventHolderI.h"

using namespace std; 
using namespace e4nu;

int main( void ) {
  std::cout << "Test ongoing..." << std::endl;
  std::string file_name = "/pnfs/genie/scratch/users/jtenavid/TestScripts/G18_02a_EM/master-routine_validation_01-eScattering/e_on_1000060120_2000MeV_0.gst.root";
  
  EventHolderI * events = new EventHolderI(file_name.c_str());
   //EventHolderI * events = new EventHolderI( file_name.c_str() );

  return 0 ; 
}
