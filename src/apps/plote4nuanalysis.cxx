// _____________________________________________________________
/* This app is used to plot the output of an e4nu analysis    */
// _____________________________________________________________

#include "plotting/AcceptanceUtils.h"
#include "plotting/XSecUtils.h"

using namespace std; 
using namespace e4nu;
using namespace e4nu::plotting;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////
// 1) MC files location                                        //
// 2) Data files location                                      //
// 3) OutputLocation                                           //
// 4) Input MC Files                                           //
// 5) Input data file                                          // 
// 6) Observables List                                         //
// 7) Input Model Names (*)                                    //
// 8) Title (*)                                                //
// 9) Data Name (*)                                            //
// 10) Compute systematics (* def false )                      //
// 11) AddSystematics                                          //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting e4nu analysis..." << std::endl;

  if( argc <= 1 ) { 
    std::cout << " Missing arguments " << std::endl;
    std::cout << " plote4nuanalysis --mc_location <mc_location> \n --data_location <data_location> \n --output_location <output_loc> \n ";
    std::cout << " --input_mc_files <file1,file2,...,fileN> \n --input_data_file <data> --observable_list <obs1,obs2,...,obsM> " << std::endl;
    std::cout << " optional arguments are : --model_names <name1,name2,...,nameN> --title <title> --data_name <data> --systematics --add-systematics name,value:name2,value2:...:nameK,valueK" << std::endl;
    return 0 ; 
  } else { 
    //Loading configuration
    std::cout << " OPTIONS EXIST " << std::endl;
    /*while(argc>1) {
    if(argv[1][0] == '-') {
      if (argv[1][1] == op) set = true;
    }
    argc--;
    argv++;
    */
  }
    
  //  if ( argv[1] ) config_file = argv[1] ;

  return 0 ; 
}
