#ifndef _SYST_UTILS_H_
#define _SYST_UTILS_H_

#include "plotting/PlottingUtils.h"

namespace e4nu {
  namespace systematics {
    // This function calculates the corresponding xsec distributions for a given observable
    // No acceptance correction is used
    void ComputeHistSyst( std::vector<std::string> input_files, std::vector<std::string> tags, std::string observable, bool is_data, 
			  std::string input_location, std::string output_location, std::string analysis_id );
  }
}

#endif


