#ifndef _XSEC_UTILS_H_
#define _XSEC_UTILS_H_

#include "plotting/PlottingUtils.h"
#include "plotting/AcceptanceUtils.h"

namespace e4nu {
  namespace plotting {
    // Input paramters:
    // MC_file_name : true MC file, without detector effects, after e4nu analysis
    // data_file_name: data file, after e4nu Analysis
    // acceptance_file_name: acceptance file obtained with compute_acceptance.C
    // target target_pdg
    // beam energy
    // Number of events in original MC file (pre-analysis)
    void Plot1DXSec(std::vector<std::string> MC_files_name, std::string data_file_name,
		    std::string acceptance_file_name, std::string observable,
		    std::string title, std::string data_name, std::vector<std::string> model,
		    std::string input_MC_location, std::string input_data_location, std::string output_location,
		    std::string output_file_name, bool plot_data ) ; 
  }
}

#endif



