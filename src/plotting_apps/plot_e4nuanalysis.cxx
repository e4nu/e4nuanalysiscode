// _____________________________________________________________
/* This app is used to plot the output of an e4nu analysis    */
// _____________________________________________________________

#include "plotting/AcceptanceUtils.h"
#include "plotting/XSecUtils.h"
#include "plotting/Systematics.h"
#include <iomanip>
#include <filesystem>
#include <map>
#include <vector>

using namespace std;
using namespace e4nu;
using namespace e4nu::plotting;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////
// 1) MC files location                                        //
// 2) Data files location (*)                                  //
// 3) OutputLocation                                           //
// 4) Input MC Files (without _true.root or _truereco.root)    //
// 5) Input data file (without _clas6.root extension) (*)      //
// 6) Observables List (X axis)                                //
// 7) Observables List (Y axis) - Optional                     //
// 8) Slices for Y axis (mandatory if Y axis used)             //
// 9) Normalize to cross-section                               //
// 10) Plot MC, default true                                   //
// 11) Input Model Names (*)                                   //
// 12) Title (*)                                               //
// 13) Data Name (*)                                           //
// 14) Compute systematics (* def false )                      //
// 15) AddSystematics                                          //
// 16) NoFSI ROOT file (without .root extension)               //
// 17) Analysis id - used to set the ranges. default none      //
//     i.e. 1p1pim                                             //
// 18) Plot root file output. Stores individual histograms     //
// 19) Apply cut on observable                                 //
//     --cut-observables Obs1,min1,max1:...:ObsN,minN,maxN     //
// 20) log-scale : option to use log scale                     //
// 21) mott-scalet: scale by mott (Q4)                         //
// 22) units : mb or nb                                        //
// 23) scale : scaling factor to multiply data and MC          //
/////////////////////////////////////////////////////////////////

string mc_location="", data_location="", output_location ="", output_name ="", analysis_id="default";
vector<string> mc_files, rad_files,observables_x, observables_y, model_names ;
string data_file ="", nofsi_file="", title="", data_name="" ;
bool compute_systematics ;
string bkg_syst;
map<string,double> systematic_map ;
bool plot_data = true ;
bool store_root = false ;
bool log_scale = false ;
bool mott_scale = false ;
double scaling = 1;
std::string units = "mb";
void PrintFormat(string s);
std::vector<double> y_cuts ; // for 2D slicing.

// We want to add a map which contains the observable to cut on and its range
std::map<string,vector<double>> cuts;

int main( int argc, char* argv[] ) {

  std::cout << "Plotting e4nu analysis..." << std::endl;

  if( argc <= 1 ) {
    PrintFormat("");
    return 0 ;
  } else {
    //Loading configuration
    if( ExistArg("data_location",argc,argv)) {
      data_location = GetArg("data_location",argc,argv) ;
    } else plot_data = false ;

    if( ExistArg("output_location",argc,argv)) {
      output_location = GetArg("output_location",argc,argv) ;

      if( ! std::filesystem::exists(output_location) ) std::filesystem::create_directory(output_location);
    } else PrintFormat("output_location") ;

    if( ExistArg("output_name",argc,argv)) {
      output_name = GetArg("output_name",argc,argv) ;
    } else PrintFormat("output_name") ;

    if( ExistArg("mc_location",argc,argv)) {
      mc_location = GetArg("mc_location",argc,argv) ;
      cout << "Reading MC files from " << mc_location << std::endl;
    } else PrintFormat("mc_location") ;

    if( ExistArg("input_mc_files",argc,argv)) {
      string files ;
      files = GetArg("input_mc_files",argc,argv) ;
      mc_files = SplitString(files);
    } else PrintFormat("input_mc_files") ;

    if( ExistArg("input_rad_files",argc,argv)) {
      string files ;
      files = GetArg("input_rad_files",argc,argv) ;
      rad_files = SplitString(files);
    } else PrintFormat("input_rad_files") ;

    if( ExistArg("input_data_file",argc,argv)) {
      data_file = GetArg("input_data_file",argc,argv) ;
      cout << "Reading data file from " << data_location << data_file << std::endl;
    } else plot_data = false ;

    if( ExistArg("analysis_id",argc,argv)){
      analysis_id = GetArg("analysis_id",argc,argv);
    }

    string obs_x ;
    if( ExistArg("observable_list",argc,argv)) {
      obs_x = GetArg("observable_list",argc,argv) ;
      observables_x = SplitString(obs_x);
      cout << "Plotting xsec for the following observables: \n- ";
      for( unsigned s = 0 ; s<observables_x.size();++s ) {
        cout << observables_x[s] << std::endl;
      }
    } else PrintFormat("observable_list") ;

    string obs_y;
    if( ExistArg("observable_list_y",argc,argv)) {
      obs_y = GetArg("observable_list_y",argc,argv) ;
      observables_y = SplitString(obs_y);
      cout << "Plotting xsec for the following observables: \n- ";
      for( unsigned s = 0 ; s<observables_y.size();++s ) {
        cout << observables_y[s] << std::endl;
      }
      string obs_cuts_y;
      if( ExistArg("observable_y_cuts",argc,argv)) {
        obs_cuts_y = GetArg("observable_y_cuts",argc,argv) ;
        auto string_cuts_y = SplitString(obs_cuts_y);
        for( unsigned s = 0 ; s<string_cuts_y.size();++s ) {
          y_cuts.push_back(stod(string_cuts_y[s]));
        }
      } else {
        std::cout << "Using 2D plotting, needs slice information."<< std::endl;
        return 0;
      }

    } else PrintFormat("observable_list") ;

    string mdl ;
    if( ExistArg("model_names",argc,argv)) {
      mdl = GetArg("model_names",argc,argv) ;
      model_names = SplitString(mdl);
      if( model_names.size() != mc_files.size() ){
        std::cout << "Number of mc files does not match the number of models"<< std::endl;
        return 0;
      }
    }
    if( ExistArg("store_root",argc,argv)) store_root = true ;
    if( ExistArg("log-scale",argc,argv)) {
      std::cout << " Enabling log scale..." << std::endl;
      log_scale = true ;
    }
    if( ExistArg("mott-scale",argc,argv)) {
      std::cout << " Enabling mott scaling..." << std::endl;
      mott_scale = true ;
    }

    cout<<"Loading MC Files:"<<endl;
    for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
      cout << " -> " ;
      if( model_names.size() != 0 ) cout << model_names[i] << ": " ;
      cout << mc_files[i] << endl;
    }

    if( ExistArg("title",argc,argv)) {
      title = GetArg("title",argc,argv) ;
    }

    if( ExistArg("data_name",argc,argv)) {
      data_name = GetArg("data_name",argc,argv) ;
    }

    if( ExistArg("nofsi_file",argc,argv)) {
      nofsi_file = GetArg("nofsi_file",argc,argv) ;
      cout << "Adding No FSI file: " << nofsi_file << endl;
    }

    string sys ;
    if( ExistArg("add-systematics",argc,argv)) {
      sys = GetArg("add-systematics",argc,argv) ;
      vector<string> sys_names = SplitString(sys,',');
      for( unsigned s = 0 ; s < sys_names.size() ; ++s ) {
        std::cout << sys_names[s]<<std::endl;
        vector<string> tmpsys = SplitString(sys_names[s],':') ;
        if( tmpsys.size() != 2 ) continue ;
        systematic_map[tmpsys[0]] = stod(tmpsys[1]);
      }
    }

    if( ExistArg("bkg-systematics",argc,argv)) {
      bkg_syst = GetArg("bkg-systematics",argc,argv);
    }

    string cut_obs ;
    if( ExistArg("cut-observables",argc,argv)) {
      cut_obs = GetArg("cut-observables",argc,argv) ;
      vector<string> obs_names = SplitString(cut_obs,':');
      for( unsigned s = 0 ; s < obs_names.size() ; ++s ) {
        vector<string> obs_details = SplitString(obs_names[s],',');
        if( obs_details.size() != 3 ) {
          std::cout << " Need three entries for observable " << obs_names[s] << std::endl;
          continue ;
        }
        string observable_s = obs_details[0] ;
        double min_obs = stod( obs_details[1] );
        double max_obs = stod( obs_details[2] );
        // And add in map
        cuts[observable_s] = {min_obs, max_obs};
        // print out information:
        std::cout << " Adding cuts on " << observable_s << " : " << std::endl;
        std::cout << " --> Min: " << min_obs << std::endl;
        std::cout << " --> Max: " << max_obs << std::endl;
      }
    }

    if( ExistArg("units",argc,argv)) {
      units = GetArg("units",argc,argv) ;
      if (units != "mb" && units != "nb") units = "mb";
    }

    if( ExistArg("scaling",argc,argv)) {
      scaling = stod(GetArg("scaling",argc,argv)) ;
    }

  }

  // Fill 2D Graph for oscillation study:
  // We only fill it with the first MC sample
  plotting::GetMissingEnergyGraph( (mc_location+"/"+mc_files[0]+"_true.root").c_str() );

  // Loop over observables_x
  for( unsigned int i = 0 ; i < observables_x.size(); ++i ){
    vector<string> root_files = mc_files;
    vector<string> root_rad_files = rad_files;
    vector<string> names = model_names ;
    string acceptance_file_1D = Compute1DAcceptance( root_files, observables_x[i], title, mc_location, output_location, output_name, cuts, analysis_id, store_root ) ;

    // compute 2D acceptance if requested
    string acceptance_file_2D = "", acceptance_file_slices = "";
    if( observables_y.size() > 0 ){
      // For now it only works for one single alternative observable
      std::cout << " Computing 2D acceptance: "<<std::endl;
      acceptance_file_2D = Compute2DAcceptance( root_files, observables_x[i], observables_y[0], title, mc_location, output_location, output_name, cuts, analysis_id, store_root ) ;
      acceptance_file_slices = Compute2DAccCorrSlice(root_files, observables_x[i], observables_y[i], y_cuts, title, mc_location, output_location, output_name, cuts, analysis_id, store_root );

    }

    string radcorr_file = "";
    if( rad_files.size() != 0 ) radcorr_file = Compute1DRadCorr( root_files, root_rad_files, observables_x[i], title, mc_location, output_location, output_name, cuts, analysis_id, store_root ) ;
    if( nofsi_file != "" ) { root_files.push_back(nofsi_file); names.push_back("No FSI");}

    // compute 2D rad corr if requested
    string radcorr_file_2D = "", radcorr_file_slices = "";
    if( observables_y.size() > 0 && rad_files.size() != 0 ){
      // For now it only works for one single alternative observable
      std::cout << " Computing 2D rad corr: "<<std::endl;
      radcorr_file_2D = Compute2DRadCorr( root_files, root_rad_files, observables_x[i], observables_y[i], title, mc_location, output_location, output_name, cuts, analysis_id, store_root ) ;
      radcorr_file_slices = Compute2DRadCorrSlice(root_files, root_rad_files, observables_x[i], observables_y[i], y_cuts, title, mc_location, output_location, output_name, cuts, analysis_id, store_root );

    }

    Plot1DXSec( root_files, data_file, acceptance_file_1D, radcorr_file, observables_x[i], title, data_name, names, mc_location, data_location, output_location, output_name, plot_data, systematic_map, bkg_syst, cuts, analysis_id, store_root, log_scale, mott_scale, units, scaling ) ;

    if( observables_y.size() > 0 ){
      Plot2DXSec( root_files, data_file, acceptance_file_2D, radcorr_file_2D, observables_x[i], observables_y[i], y_cuts, title, data_name, names, mc_location, data_location, output_location, output_name, plot_data, systematic_map, cuts, analysis_id, units, store_root, log_scale, mott_scale, scaling ) ;
      Plot2DSlicesXSec(root_files, data_file, acceptance_file_slices, radcorr_file_slices, observables_x[i], observables_y[i], y_cuts, title, data_name, names, mc_location, data_location, output_location, output_name, plot_data, systematic_map, bkg_syst, cuts, analysis_id, units, store_root, log_scale, mott_scale, scaling);
    }
  }

  return 0 ;
}

void PrintFormat(string s){
  if( s!="") cout << " Missing " << s << endl;
  cout << "RequiredArguments:\n";
  cout << " plote4nuanalysis --mc_location <mc_location> \n --data_location <data_location> \n --output_location <output_loc> \n --output_name <outputname> \n ";
  cout << "--input_mc_files <file1,file2,...,fileN> \n --input_data_file <data> --observable_list <obs1,obs2,...,obsM> " << endl;
  cout << " optional arguments are : \n --model_names <name1,name2,...,nameN> \n --title <title> \n --data_name <data> \n --systematics \n";
  cout << " --add-systematics name,value:name2,value2:...:nameK,valueK \n --nofsi_file <rootfile>" << endl;
}
