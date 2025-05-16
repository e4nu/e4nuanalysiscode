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
	// hist : histogram which we want to add the syst. unc. to.
	// hist_w_error: the content is the per-cent relative uncertanty to add.
	// Binning must coincide.

  TH1D * hist_syst = (TH1D*)hist_w_error.Clone();
  hist_syst->Reset();
  hist_syst->SetName("h_syst");
  hist_syst->GetYaxis()->SetTitle("#sigma/#hat{x} [%]");

  double NBins = hist.GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
		// We are storing the error in a histogram (hist_w_error). the content is the error itself to add to the stat uncertanty from hist.
    double stat_error = hist.GetBinError(i);
    double syst_error = hist_w_error.GetBinContent(i)/100;
		
		// The error calculation assumes that it is a multiplicative factor
		// The error is added as a relative
		double newerror = TMath::Sqrt( TMath::Power(stat_error,2.) + TMath::Power(syst_error*hist.GetBinContent(i),2.));
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
    if( error_2 < 0 ) error_2 = 0 ;

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
