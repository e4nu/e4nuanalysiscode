#ifndef _ACCEPTANCE_UTILS_H_
#define _ACCEPTANCE_UTILS_H_

#include "plotting/PlottingUtils.h"

namespace e4nu {
  namespace plotting {

    // This function computes the acceptance for a given observable
    // It uses the mc file location to find the _true and _truereco.root files
    std::string ComputeAcceptance(std::vector<std::string> mc_files, std::string observable, std::string title,
				  std::string input_MC_location, std::string output_location,  std::string output_file_name ) ; 

  }
}

#endif


