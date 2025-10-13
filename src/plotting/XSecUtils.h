
#ifndef _XSEC_UTILS_H_
#define _XSEC_UTILS_H_

#include <string>
#include <iostream>
#include "plotting/PlottingUtils.h"
#include "plotting/AcceptanceUtils.h"

using namespace std;
namespace e4nu {
  namespace plotting {

    void Plot1DXSec(std::vector<std::string> MC_files_name, std::string data_file_name,
		    std::string acceptance_file_name, std::string radcorr_file, std::string observable,
		    std::string title, std::string data_name, std::vector<std::string> model,
		    std::string input_MC_location, std::string input_data_location, std::string output_location,
		    std::string output_file_name, bool plot_data, std::map<string,double> systematic_map, string bkg_syst,
        std::map<std::string,std::vector<double>> cuts,
		    std::string analysis_id = "default", bool store_root = false, bool log_scale = false, bool scale_mott = false ) ;

    void Plot2DXSec(std::vector<std::string> MC_files_name, std::string data_file_name,
        std::string acceptance_file_name, std::string radcorr_file, std::string x_observable, std::string y_observable,
        std::vector<double> & y_cuts, std::string title, std::string data_name, std::vector<std::string> model,
        std::string input_MC_location, std::string input_data_location, std::string output_location,
        std::string output_file_name, bool plot_data, std::map<string,double> systematic_map,
        std::map<std::string,std::vector<double>> cuts,
        std::string analysis_id = "default", bool store_root = false, bool log_scale = false, bool scale_mott = false ) ;

    void PlotXsecDataTotal( TH1D * data, std::string observable, std::string title, std::string data_name,
			    std::string input_data_location, std::string output_location,
			    std::string output_file_name, std::map<string,double> systematic_map,
			    std::string analysis_id, bool store_root = false, bool log_scale = false ) ;

    void PlotComparisonDataNormalized( std::vector<TH1D> mc_hists, std::vector<TH1D> breakdown,
        	TH1D * data, std::string observable, std::string title, std::string data_name, std::vector<std::string> model,
          std::string input_MC_location, std::string input_data_location, std::string output_location,
          std::string output_file_name, std::map<string,double> systematic_map, bool show_breakdown = true,
          std::string analysis_id= "default", bool store_root = false, bool log_scale = false) ;

    void PlotTotalXSec( std::vector<TH1D*> mc_hists, std::vector<TH1D*> breakdown, TH1D * data,
			std::string observable, std::string title, std::string data_name, std::vector<std::string> model,
			std::string input_MC_location, std::string input_data_location, std::string output_location,
			std::string output_file_name, std::map<string,double> systematic_map, bool show_breakdown = true,
			std::string analysis_id = "default", bool store_root = false, bool log_scale = false, std::string slice_title = "") ;

    void PlotTotalXSec( std::vector<TH1D*> mc_hists, TH1D * data, std::string observable, std::string title, std::string data_name, std::vector<std::string> model,
      std::string input_MC_location, std::string input_data_location, std::string output_location,
      std::string output_file_name, std::map<string,double> systematic_map, bool show_breakdown = true,
      std::string analysis_id = "default", bool store_root = false, bool log_scale = false, std::string slice_title = "") ;

    void PlotTotal2DXSec( std::vector<TH2D*> mc_hists, std::vector<TH2D*> breakdown, TH2D * data,
  		std::string x_observable, std::string y_observable, std::vector<double> & y_cuts, std::string title, std::string data_name,
      std::vector<std::string> model, std::string input_MC_location,
  		std::string input_data_location, std::string output_location,
  		std::string output_file_name, std::map<string,double> systematic_map, bool show_breakdown = true,
  		std::string analysis_id = "default", bool store_root = false, bool log_scale = false) ;

    void PlotEventRate( TH1D * data, std::string observable, std::string title, std::string data_name, std::string input_data_location,
			std::string output_location, std::string output_file_name, std::string analysis_id, bool store_root, bool log_scale = false ) ;

    void PlotEventRatePerSector( std::vector<TH1D*> data_per_sector, std::string observable, std::string title, std::string data_name, std::string input_data_location,
      std::string output_location, std::string output_file_name, std::string analysis_id, bool store_root, bool log_scale = false ) ;

    void PlotLegend( std::vector<TH1D*> mc_hists, std::vector<TH1D*> breakdown, TH1D * data, std::string observable,
		     std::string data_name, std::vector<std::string> model,std::string output_location,
		     std::string output_file_name, bool store_root = false, bool log_scale = false) ;


    void PlotPerSector( std::vector<TH1D*> mc_per_sector,std::vector<TH1D*> data_per_sector, std::string observable,
			std::string title, std::string data_name, std::vector<std::string> model,
			std::string input_MC_location, std::string input_data_location, std::string output_location,
			std::string output_file_name, std::map<string,double> systematic_map,
			std::string analysis_id = "default", bool store_root = false, bool log_scale = false) ;

    void CreateCanvasWithPads(TPad *& pad, std::vector<TPad*>& topPad, std::vector<TPad*>& bottomPad, const std::string& canvasName);
    void PlotProjection( const std::vector<TH2D*>& mcHists, const std::vector<TH2D*>& data, const std::string& xobservable, const std::string& yobservable, TCanvas* canvas, bool logScale, const std::string& axis, double y_cut_min, double y_cut_max );
    void PlotProjectionWithRatio( const std::vector<TH2D*>& mcHists, const std::vector<TH2D*>& breakdown, TH2D* data, TH2D* acceptance, TH2D* radcorr, const std::string& xobservable, const std::string& yobservable, TPad* topPad, TPad* bottomPad, bool logScale, const std::string& axis, double y_cut_min, double y_cut_max) ;
    void PlotProjectionsStack( const std::vector<TH2D*>& mcHists, const std::vector<TH2D*>& breakdown, const TH2D* data, const TH2D* acceptance, const TH2D* radcorr, const std::string& xobservable, const std::string& yobservable, bool logScale, const std::string& axis, std::vector<double> y_cuts, const std::string& outputLocation, const std::string& outputName, bool store_root );
    void PlotSlicesGeneralized(const std::vector<TH2D*>& mcHists, const std::vector<TH2D*>& breakdown, TH2D* data, TH2D* acceptance, TH2D* radcorr,
        const std::string& xObservable, const std::string& yObservable, std::vector<double> & y_cuts, const std::string& title, const std::string& outputLocation,
        const std::string& outputName, bool storeRoot, bool logScale );
    void Plot2DSlicesXSec(std::vector<std::string> MC_files_name, std::string data_file_name, std::string acceptance_file_name, std::string radcorr_file, std::string x_observable, std::string y_observable, std::vector<double> & y_cuts, std::string title, std::string data_name, std::vector<std::string> model, std::string input_MC_location, std::string input_data_location, std::string output_location, std::string output_file_name, bool plot_data, std::map<string, double> systematic_map, string bkg_syst, std::map<std::string,std::vector<double>> cuts, std::string analysis_id, bool store_root, bool log_scale, bool scale_mott);

  }
}

#endif
