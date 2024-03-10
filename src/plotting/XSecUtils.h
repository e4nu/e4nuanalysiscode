#ifndef _XSEC_UTILS_H_
#define _XSEC_UTILS_H_

#include <string>
#include <iostream>
#include "plotting/PlottingUtils.h"
#include "plotting/AcceptanceUtils.h"

using namespace std;
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
		    std::string acceptance_file_name, std::string radcorr_file, std::string observable,
		    std::string title, std::string data_name, std::vector<std::string> model,
		    std::string input_MC_location, std::string input_data_location, std::string output_location,
		    std::string output_file_name, bool plot_data, std::map<string,double> systematic_map, 
		    std::string analysis_id = "default", bool store_root = false) ; 

    void PlotXsecDataTotal( TH1D * data, std::string observable, std::string title, std::string data_name,
			    std::string input_data_location, std::string output_location,
			    std::string output_file_name, std::map<string,double> systematic_map,
			    std::string analysis_id, bool store_root = false ) ;
 
    void PlotTotalXSec( std::vector<TH1D*> mc_hists, std::vector<TH1D*> breakdown, TH1D * data, 
			std::string observable, std::string title, std::string data_name, std::vector<std::string> model,
			std::string input_MC_location, std::string input_data_location, std::string output_location,
			std::string output_file_name, std::map<string,double> systematic_map, 
			std::string analysis_id = "default", bool store_root = false) ;
    
    void PlotEventRate( TH1D * data, std::string observable, std::string title, std::string data_name, std::string input_data_location, 
			std::string output_location, std::string output_file_name, std::string analysis_id, bool store_root ) ;

    void PlotLegend( std::vector<TH1D*> mc_hists, std::vector<TH1D*> breakdown, TH1D * data, std::string observable,
		     std::string data_name, std::vector<std::string> model,std::string output_location,
		     std::string output_file_name, bool store_root = false) ;


    void PlotPerSector( std::vector<TH1D*> mc_per_sector,std::vector<TH1D*> data_per_sector, std::string observable,
			std::string title, std::string data_name, std::vector<std::string> model,
			std::string input_MC_location, std::string input_data_location, std::string output_location,
			std::string output_file_name, std::map<string,double> systematic_map,
			std::string analysis_id = "default", bool store_root = false) ;

    void PlotSlices( std::vector<std::vector<TH1D*>> all_slices, std::vector<double> addbinning, std::string observable,
		     std::string title, std::string data_name, std::vector<std::string> model,
		     std::string input_MC_location, std::string input_data_location, std::string output_location,
		     std::string output_file_name, std::map<string,double> systematic_map,
		     std::string analysis_id = "default", bool store_root = false ) ;
  }
}

#endif



