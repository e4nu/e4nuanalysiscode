#ifndef _SYST_UTILS_H_
#define _SYST_UTILS_H_

#include "plotting/PlottingUtils.h"

namespace e4nu {
  namespace systematics {
    // This function calculates the corresponding xsec distributions for a given observable
    // No acceptance correction is used
    void ComputeHistSyst( std::vector<std::string> input_files, std::vector<std::string> tags, std::string observable, bool is_data,
			                    std::string input_location, std::string output_location, std::string analysis_id );
    void AddSystematic( TH1D & hist, const double rel_error, const std::string name ) ;
    void AddSystematic( TH1D & hist, const TH1D & hist_w_error ) ;
    TH1D * SectorVariationError( TH1D & hist, const std::vector<TH1D*> h_per_sector ) ; 
  }
}

#endif
