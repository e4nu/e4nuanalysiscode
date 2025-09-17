#ifndef _ACCEPTANCE_UTILS_H_
#define _ACCEPTANCE_UTILS_H_

#include "plotting/PlottingUtils.h"

namespace e4nu {
  namespace plotting {
    // This function computes the acceptance for a given observable
    // It uses the mc file location to find the _true and _truereco.root files
    std::string Compute1DAcceptance(std::vector<std::string> mc_files, std::string observable, std::string title, std::string input_MC_location, std::string output_location, std::string output_file_name, std::map<std::string,std::vector<double>> cuts, std::string analysis_id = "default", bool store_root=false) ;

    std::string Compute2DAcceptance(std::vector<std::string> mc_files, std::string x_observable, std::string y_observable, std::string title, std::string input_MC_location, std::string output_location, std::string output_file_name, std::map<std::string,std::vector<double>> cuts, std::string analysis_id = "default", bool store_root=false) ;

    std::string Compute1DRadCorr(std::vector<std::string> mc_files, std::vector<std::string> rad_files, std::string observable, std::string title, std::string input_MC_location, std::string output_location,  std::string output_file_name, std::map<std::string,std::vector<double>> cuts, std::string analysis_id="default", bool store_root=false) ;

    std::string Compute2DRadCorr(std::vector<std::string> mc_files, std::vector<std::string> rad_files, std::string x_observable, std::string y_observable, std::string title, std::string input_MC_location, std::string output_location, std::string output_file_name, std::map<std::string,std::vector<double>> cuts, std::string analysis_id = "default", bool store_root=false) ;
    std::string Compute2DAccCorrSlice(std::vector<std::string> mc_files, std::string x_observable, std::string y_observable, std::vector<double> y_cuts, std::string title, std::string input_MC_location, std::string output_location, std::string output_file_name, std::map<std::string,std::vector<double>> alt_cuts, std::string analysis_id, bool store_root ) ;
    std::string Compute2DRadCorrSlice(std::vector<std::string> mc_files, std::vector<std::string> rad_files, std::string x_observable, std::string y_observable, std::vector<double> y_cuts, std::string title, std::string input_MC_location, std::string output_location, std::string output_file_name, std::map<std::string,std::vector<double>> alt_cuts, std::string analysis_id, bool store_root ) ;
  }
}

#endif
