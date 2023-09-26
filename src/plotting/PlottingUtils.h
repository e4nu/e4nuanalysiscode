#ifndef _PLOTTING_UTILS_H_
#define _PLOTTING_UTILS_H_

#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"
#include <iomanip>

namespace e4nu {
  namespace plotting {
    void NormalizeHist( TH1D * h, double normalization_factor );
    void CorrectData(TH1D* h, TH1D* acc);
    std::string GetAxisLabel( std::string observable, unsigned int id_axis );
    std::vector<double> GetUniformBinning( unsigned int nbins, double min, double max);
    std::vector<double> GetECalBinning( unsigned int nbins_tail, unsigned int nbins_peak, double min, double max, double EBeam);
    std::vector<double> GetBinning( std::string observable, double EBeam );
    std::vector<double> GetAdditionalBinning( std::string second_observable, double EBeam ) ;
    std::string GetAlternativeObs( std::string observable );
    std::string GetObsName( std::string observable );
    std::string GetUnit( std::string observable );
    double GetMaximum( std::vector<TH1D*> predictions);
    void StandardFormat( TH1D * prediction, std::string title, int color, int style, std::string observable, double y_max = 0 );
    vector<string> SplitString(string s, char d ) ;
    string GetArg(string op, int argc, char ** argv );
    bool ExistArg(string op, int argc, char ** argv );
  }
}

#endif



