#include "plotting/Systematics.h"
#include "plotting/PlottingUtils.h"
#include "plotting/XSecUtils.h"
#include "TLegend.h"
#include <iomanip>
#include "TMath.h"
#include <filesystem>
#include <sstream>
#include <iostream>
#include <string>

using namespace e4nu ;
using namespace e4nu::systematics ;
using namespace e4nu::plotting;

void systematics::ComputeHistSyst( std::vector<std::string> input_files, std::vector<std::string> tags, std::string observable, bool is_data,
				   std::string input_location, std::string output_location, std::string analysis_id ){
  if( input_files.size() ==0 ) return ;

  std::string name_fullpath = output_location+"/"+input_files[0]+"_"+observable+"_systematics.root";
  std::unique_ptr<TFile> myFile( TFile::Open(name_fullpath.c_str(), "RECREATE") );

  std::vector<TFile*> ifiles ;
  for( unsigned int id = 0 ; id < input_files.size(); ++id ){
    ifiles.push_back(new TFile((input_location+"/"+input_files[id]+".root").c_str(),"ROOT"));
    if( !ifiles[id] ) { std::cout << "ERROR: the "<< input_location<<input_files[id]<<".root does not exist." << std::endl; return ;}
  }

  // Get Tree for main model
  std::string treename = "MCCLAS6Tree";
  if( is_data ) treename = "CLAS6Tree";

  std::vector<TH1D*> hists;
  for( unsigned int i = 0 ; i < input_files.size(); ++i ) {
    double BeamE ;
    TTree * tree = (TTree*)ifiles[i]->Get(treename.c_str());
    if( !tree ) return ;
    tree ->SetBranchAddress("BeamE",&BeamE);
    tree->GetEntry(0);
    std::vector<double> binning = plotting::GetBinning(observable,BeamE,analysis_id);
    if( binning.size()==0 ) return ;

    // Create histogram for total and total xsec per sector
    TH1D * hist = new TH1D( tags[i].c_str(), tags[i].c_str(), binning.size()-1, &binning[0] ) ;
    if( !hist ) return ;

    // OBSERVABLE DEFINITION:
    std::vector<double> mc_norm ;

    if( !tree ) continue ;
    plotting::SetAnalysisBranch( tree ) ;

    double norm =0;
    for( int j = 0 ; j < NEntries ; ++j ) {
      tree->GetEntry(j) ;
      double content = 0 ;
      double w = TotWeight ;
      content = GetObservable(observable);

      norm = DataNormalization ;
      if( !is_data ) norm = MCNormalization ;

      hist -> Fill( content, w ) ;
      hist -> SetLineWidth(3);
    }// entries loop

    plotting::NormalizeHist(hist,norm);
    myFile->WriteObject(hist,tags[i].c_str());
    hists.push_back(hist);

  } // input files

  if( hists.size() > 1 ) {
    for( unsigned int i = 1 ; i < hists.size() ; ++i ) {
      TH1D * diff = (TH1D*) hists[0]->Clone();
      diff -> Add( hists[i], -1);
      diff -> SetName( ("diff_def_vs_"+std::to_string(i)).c_str()) ;
      myFile->WriteObject(diff, ("diff_def_vs_"+std::to_string(i)).c_str() );
    }
  }

}


void systematics::AddSystematic( TH1D & hist, const double rel_error, const std::string name ) {
  double NBins = hist.GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
    double error = hist.GetBinError(i);
    double content = hist.GetBinContent(i);
    double newerror = TMath::Sqrt( TMath::Power(error,2.) + TMath::Power(rel_error*content/100.,2.));
    hist.SetBinError(i,newerror);
  }
}

void systematics::AddSystematic( TH2D & hist, const double rel_error, const std::string name ) {
  double NBins = hist.GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
    double error = hist.GetBinError(i);
    double content = hist.GetBinContent(i);
    double newerror = TMath::Sqrt( TMath::Power(error,2.) + TMath::Power(rel_error*content/100.,2.));
    hist.SetBinError(i,newerror);
  }
}


TH1D * systematics::AddSystematic( TH1D & hist, const TH1D & hist_w_error ) {
  TH1D * hist_syst = (TH1D*)hist_w_error.Clone();
  hist_syst->Reset();
  hist_syst->SetName("h_syst");
  hist_syst->GetYaxis()->SetTitle("#sigma/#hat{x} [%]");

  double NBins = hist.GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
    double stat_error = hist.GetBinError(i);
    double syst_error = hist_w_error.GetBinError(i);

		// The error calculation assumes that it is a multiplicative factor
		// The error is added as a relative
		double newerror = TMath::Sqrt( TMath::Power(stat_error,2.) + TMath::Power(syst_error*hist.GetBinContent(i)/hist_w_error.GetBinContent(i),2.));
		if( hist.GetBinContent(i) > 0 ) hist.SetBinError(i,newerror);
    else hist.SetBinError(i,0); // remove odd bins

    hist_syst->SetBinContent(i,hist_w_error.GetBinError(i)/hist_w_error.GetBinContent(i)*100); // Store sector to sector uncertainty
  }
  return hist_syst;
}


TH2D * systematics::AddSystematic( TH2D & hist, const TH2D & hist_w_error ) {
  TH2D * hist_syst = (TH2D*)hist_w_error.Clone();
  hist_syst->Reset();
  hist_syst->SetName("h_syst");
  hist_syst->GetYaxis()->SetTitle("#sigma/#hat{x} [%]");

  double NBins = hist.GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
    double stat_error = hist.GetBinError(i);
    double syst_error = hist_w_error.GetBinError(i);

		// The error calculation assumes that it is a multiplicative factor
		// The error is added as a relative
		double newerror = TMath::Sqrt( TMath::Power(stat_error,2.) + TMath::Power(syst_error*hist.GetBinContent(i)/hist_w_error.GetBinContent(i),2.));
		if( hist.GetBinContent(i) > 0 ) hist.SetBinError(i,newerror);
    else hist.SetBinError(i,0); // remove odd bins

    hist_syst->SetBinContent(i,hist_w_error.GetBinError(i)/hist_w_error.GetBinContent(i)*100); // Store sector to sector uncertainty
  }
  return hist_syst;
}

TH1D * systematics::SectorVariationError( TH1D & hist, const std::vector<TH1D*> h_per_sector ) {

  TH1D * hist_syst_sector = (TH1D*)h_per_sector[0]->Clone();
  hist_syst_sector->Reset();
  hist_syst_sector->SetName("h_syst_sectors");
  hist_syst_sector->GetYaxis()->SetTitle("#sigma_{Sector}/#hat{x}[%]");

  // Loop over bins (starting from 1 )
  for( unsigned j = 1 ; j < hist.GetNbinsX() +1; ++j ){

    double mean_i = 0 ;
    double weight = 0 ;
    double sectors = 0;
    // Loop over sectors to compute the mean variance per sector in each bin
    for( unsigned int i = 0 ; i < h_per_sector.size() ; ++i ){
      // We have to be careful with empty sectors:
      if( h_per_sector[i]->GetBinContent(j) != 0 ){
			  if( h_per_sector[i]->GetBinError(j) != 0 ) {
					mean_i += h_per_sector[i]->GetBinContent(j)/pow(h_per_sector[i]->GetBinError(j),2);
			  	weight += 1./pow(h_per_sector[i]->GetBinError(j),2) ;
			  }
				sectors+= 1 ;
      }
    }

    // Compute weighted average:
    if( weight != 0 ) mean_i /= weight ;

    // Compute RMS:
    double error_2 = 0;

    for( unsigned int i = 0 ; i < h_per_sector.size() ; ++i ){
      if( h_per_sector[i]->GetBinContent(j) != 0 ){
				// Compute error
				error_2 += pow( h_per_sector[i]->GetBinContent(j) - mean_i, 2) ;
				//and subtract stat. error
				if( h_per_sector[i]->GetBinError(j) > 0 ) error_2 -= pow(h_per_sector[i]->GetBinError(j),2);
      }
    }
    if( sectors > 1 ) error_2 /= sectors - 1 ;
    if( error_2 < 0 ) error_2 = 0 ; // If stat. is bigger than syst, just assign 0

    if( mean_i != 0 ) hist_syst_sector->SetBinContent(j,sqrt(error_2)/mean_i*100); // Store sector to sector uncertainty

    error_2 += pow(hist.GetBinError(j),2); // Add statistical error from final histogram

    if( hist.GetBinContent(j) > 0 ) hist.SetBinError(j,sqrt(error_2));
    else hist.SetBinError(j,0);
  }
  return hist_syst_sector;
}


TH2D * systematics::SectorVariationError( TH2D & hist, const std::vector<TH2D*> h_per_sector ) {

  TH2D * hist_syst_sector = (TH2D*)h_per_sector[0]->Clone();
  hist_syst_sector->Reset();
  hist_syst_sector->SetName("h_syst_sectors");
  hist_syst_sector->GetYaxis()->SetTitle("#sigma_{Sector}/#hat{x}[%]");

  // Loop over bins (starting from 1 )
  for( unsigned j = 1 ; j < hist.GetNbinsX() +1; ++j ){

    double mean_i = 0 ;
    double weight = 0 ;
    double sectors = 0;
    // Loop over sectors to compute the mean variance per sector in each bin
    for( unsigned int i = 0 ; i < h_per_sector.size() ; ++i ){
      // We have to be careful with empty sectors:
      if( h_per_sector[i]->GetBinContent(j) != 0 ){
			  if( h_per_sector[i]->GetBinError(j) != 0 ) {
					mean_i += h_per_sector[i]->GetBinContent(j)/pow(h_per_sector[i]->GetBinError(j),2);
			  	weight += 1./pow(h_per_sector[i]->GetBinError(j),2) ;
			  }
				sectors+= 1 ;
      }
    }

    // Compute weighted average:
    if( weight != 0 ) mean_i /= weight ;

    // Compute RMS:
    double error_2 = 0;

    for( unsigned int i = 0 ; i < h_per_sector.size() ; ++i ){
      if( h_per_sector[i]->GetBinContent(j) != 0 ){
				// Compute error
				error_2 += pow( h_per_sector[i]->GetBinContent(j) - mean_i, 2) ;
				//and subtract stat. error
				if( h_per_sector[i]->GetBinError(j) > 0 ) error_2 -= pow(h_per_sector[i]->GetBinError(j),2);
      }
    }
    if( sectors > 1 ) error_2 /= sectors - 1 ;
    if( error_2 < 0 ) error_2 = 0 ; // If stat. is bigger than syst, just assign 0

    if( mean_i != 0 ) hist_syst_sector->SetBinContent(j,sqrt(error_2)/mean_i*100); // Store sector to sector uncertainty

    error_2 += pow(hist.GetBinError(j),2); // Add statistical error from final histogram

    if( hist.GetBinContent(j) > 0 ) hist.SetBinError(j,sqrt(error_2));
    else hist.SetBinError(j,0);
  }
  return hist_syst_sector;
}
