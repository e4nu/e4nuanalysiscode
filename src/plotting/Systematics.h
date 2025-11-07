#ifndef _SYST_UTILS_H_
#define _SYST_UTILS_H_

#include "plotting/PlottingUtils.h"

namespace e4nu {
  namespace systematics {
    // This function calculates the corresponding xsec distributions for a given observable
    // No acceptance correction is used
    void AddSystematic( TH1D & hist, const double rel_error, const std::string name ) ;
    void AddSystematic( TH2D & hist, const double rel_error, const std::string name ) ;
    void AddRadSystematic( TH1D & hist, const double rel_error, const std::string name );
    TH1D * AddSystematic( TH1D & hist, const TH1D & hist_w_error ) ;
    TH2D * AddSystematic( TH2D & hist, const TH2D & hist_w_error ) ;
    TH1D * SectorVariationError( TH1D & hist, const std::vector<TH1D*> h_per_sector ) ;
    TH2D * SectorVariationError( TH2D & hist, const std::vector<TH2D*> h_per_sector ) ;
  }
}

#endif
