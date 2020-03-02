#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>

using namespace std;

void OverlayECal_FigExtData11() {
 
	TString target = "56Fe";

	double xmin = 1.8;
	double xmax = 5;
	int N_E_bins;

	TString version = "v3_0_6/";

	int FontStyle = 132;
	double TextSize = 0.07;
	int Ndivisions = 5;

	TLatex *lat1 = new TLatex();
	TLatex *lat2 = new TLatex();

	lat1->SetNDC();
	lat2->SetNDC();

	TString PathToFile = "/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/4_461/Data_Final/NoxBCut/";
	TFile* file_in = new TFile(PathToFile+target+"_4_461_Data_Final_Plots_FSI_em.root");

	gStyle->SetOptStat(0);
	gStyle->SetPadTopMargin(0.01);
	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetNdivisions(Ndivisions,"Y");
	gStyle->SetNdivisions(Ndivisions,"X");
	gStyle->SetHistLineWidth(5);
	TGaxis::SetMaxDigits(3);

	//------------- e- only-------------------------------------------------------------

	TH1D *h1_Erec_0pi_e = (TH1D*)file_in->Get("h1_E_rec_0pi");
	TH1D *h1_E_rec = (TH1D*)file_in->Get("h1_E_rec");
	TH1D *h1_Erec_0pi_sub_e = (TH1D*)file_in->Get("h_Erec_subtruct_piplpimi_noprot_4pi");
	TH1D *h1_Erec_1pi_e_weight = (TH1D*)file_in->Get("h1_E_rec_1pi_weight");
	TH1D *h1_Erec_2pi_e_weight = (TH1D*)file_in->Get("h1_E_rec_2pi_weight");
	TH1D *h1_Erec_3pi_e_weight = (TH1D*)file_in->Get("h1_E_rec_3pi_weight");
	TH1D *h1_Erec_4pi_e_weight = (TH1D*)file_in->Get("h1_E_rec_4pi_weight");
	TH1D *h1_Erec_1pi_e = (TH1D*)file_in->Get("h1_E_rec_1pi");

	h1_Erec_0pi_e->SetLineColor(38);
	h1_Erec_0pi_sub_e->SetLineColor(38);
	h1_Erec_1pi_e_weight->SetLineColor(38);
	h1_Erec_2pi_e_weight->SetLineColor(38);
	h1_Erec_3pi_e_weight->SetLineColor(38);
	h1_Erec_4pi_e_weight->SetLineColor(38);
	h1_Erec_1pi_e->SetLineColor(38);
	h1_E_rec->SetLineColor(38); 

	N_E_bins=h1_Erec_0pi_e->GetNbinsX();
	double Scale = 1E-6;

	// ---------------------------------------------------------------------------------------------------- 

	// Account for the different bin widths

	for(int i=1;i<=N_E_bins;i++) {

		h1_Erec_0pi_e->SetBinContent(i,h1_Erec_0pi_e->GetBinContent(i)/h1_Erec_0pi_e->GetBinWidth(i));
		h1_Erec_0pi_sub_e->SetBinContent(i,h1_Erec_0pi_sub_e->GetBinContent(i)/h1_Erec_0pi_sub_e->GetBinWidth(i));
		h1_Erec_1pi_e_weight->SetBinContent(i,h1_Erec_1pi_e_weight->GetBinContent(i)/h1_Erec_1pi_e_weight->GetBinWidth(i));
		h1_Erec_2pi_e_weight->SetBinContent(i,h1_Erec_2pi_e_weight->GetBinContent(i)/h1_Erec_2pi_e_weight->GetBinWidth(i));
		h1_Erec_3pi_e_weight->SetBinContent(i,h1_Erec_3pi_e_weight->GetBinContent(i)/h1_Erec_3pi_e_weight->GetBinWidth(i));
		h1_Erec_4pi_e_weight->SetBinContent(i,h1_Erec_4pi_e_weight->GetBinContent(i)/h1_Erec_4pi_e_weight->GetBinWidth(i));
		h1_Erec_1pi_e->SetBinContent(i,h1_Erec_1pi_e->GetBinContent(i)/h1_Erec_1pi_e->GetBinWidth(i));
		h1_E_rec->SetBinContent(i,h1_E_rec->GetBinContent(i)/h1_E_rec->GetBinWidth(i));


		h1_Erec_0pi_e->SetBinError(i,h1_Erec_0pi_e->GetBinError(i)/h1_Erec_0pi_e->GetBinWidth(i));
		h1_Erec_0pi_sub_e->SetBinError(i,h1_Erec_0pi_sub_e->GetBinError(i)/h1_Erec_0pi_sub_e->GetBinWidth(i));
		h1_Erec_1pi_e_weight->SetBinError(i,h1_Erec_1pi_e_weight->GetBinError(i)/h1_Erec_1pi_e_weight->GetBinWidth(i));
		h1_Erec_2pi_e_weight->SetBinError(i,h1_Erec_2pi_e_weight->GetBinError(i)/h1_Erec_2pi_e_weight->GetBinWidth(i));
		h1_Erec_3pi_e_weight->SetBinError(i,h1_Erec_3pi_e_weight->GetBinError(i)/h1_Erec_3pi_e_weight->GetBinWidth(i));
		h1_Erec_4pi_e_weight->SetBinError(i,h1_Erec_4pi_e_weight->GetBinError(i)/h1_Erec_4pi_e_weight->GetBinWidth(i));
		h1_Erec_1pi_e->SetBinError(i,h1_Erec_1pi_e->GetBinError(i)/h1_Erec_1pi_e->GetBinWidth(i));
		h1_E_rec->SetBinError(i,h1_E_rec->GetBinError(i)/h1_E_rec->GetBinWidth(i));

	}

	// ---------------------------------------------------------------------------------------------------- 

	TCanvas *c1=new TCanvas("c1","",1024,768);
	c1->SetTopMargin(0.01);
	c1->SetRightMargin(0.01);
	c1->SetLeftMargin(0.15);
	c1->SetBottomMargin(0.16);
	c1->cd();

	// Scale to get rid of the exponent

	double Scale = 1E6;

	for(int i=1;i<=N_E_bins;i++) {
	
		h1_E_rec->SetBinContent(i,h1_E_rec->GetBinContent(i)/Scale);
		h1_Erec_0pi_e->SetBinContent(i,h1_Erec_0pi_e->GetBinContent(i)/Scale);
		h1_Erec_0pi_sub_e->SetBinContent(i,h1_Erec_0pi_sub_e->GetBinContent(i)/Scale);

		h1_Erec_1pi_e_weight->SetBinContent(i,h1_Erec_1pi_e_weight->GetBinContent(i)/ Scale);
		h1_Erec_2pi_e_weight->SetBinContent(i,h1_Erec_2pi_e_weight->GetBinContent(i)/ Scale);
		h1_Erec_3pi_e_weight->SetBinContent(i,h1_Erec_3pi_e_weight->GetBinContent(i)/ Scale);
		h1_Erec_4pi_e_weight->SetBinContent(i,h1_Erec_4pi_e_weight->GetBinContent(i)/ Scale);
		h1_Erec_1pi_e->SetBinContent(i,h1_Erec_1pi_e->GetBinContent(i)/ Scale);

		h1_E_rec->SetBinError(i,h1_E_rec->GetBinError(i)/ Scale);
		h1_Erec_0pi_e->SetBinError(i,h1_Erec_0pi_e->GetBinError(i)/ Scale);
		h1_Erec_0pi_sub_e->SetBinError(i,h1_Erec_0pi_sub_e->GetBinError(i)/ Scale);
		
		h1_Erec_1pi_e_weight->SetBinError(i,h1_Erec_1pi_e_weight->GetBinError(i)/ Scale);
		h1_Erec_2pi_e_weight->SetBinError(i,h1_Erec_2pi_e_weight->GetBinError(i)/ Scale);
		h1_Erec_3pi_e_weight->SetBinError(i,h1_Erec_3pi_e_weight->GetBinError(i)/ Scale);
		h1_Erec_4pi_e_weight->SetBinError(i,h1_Erec_4pi_e_weight->GetBinError(i)/ Scale);
		h1_Erec_1pi_e->SetBinError(i,h1_Erec_1pi_e->GetBinError(i)/ Scale);

	}

	// ---------------------------------------------------------------------------------------------------- 

	gStyle->SetErrorX(kFALSE);

	h1_Erec_1pi_e_weight->SetAxisRange(xmin, xmax, "X"); 
	h1_Erec_1pi_e_weight->GetXaxis()->SetTitle("E^{QE}[GeV]");
	h1_Erec_1pi_e_weight->GetXaxis()->SetTitleOffset(1.);

	h1_Erec_1pi_e_weight->GetXaxis()->CenterTitle();
	h1_Erec_1pi_e_weight->GetXaxis()->SetTitleFont(FontStyle);
	h1_Erec_1pi_e_weight->GetXaxis()->SetLabelFont(FontStyle);
	h1_Erec_1pi_e_weight->GetXaxis()->SetTitleSize(TextSize);
	h1_Erec_1pi_e_weight->GetXaxis()->SetLabelSize(TextSize);
	h1_Erec_1pi_e_weight->GetXaxis()->SetNdivisions(Ndivisions);

	h1_Erec_1pi_e_weight->GetYaxis()->CenterTitle();
	h1_Erec_1pi_e_weight->GetYaxis()->SetTitleFont(FontStyle);
	h1_Erec_1pi_e_weight->GetYaxis()->SetLabelFont(FontStyle);
	h1_Erec_1pi_e_weight->GetYaxis()->SetTitleSize(TextSize);
	h1_Erec_1pi_e_weight->GetYaxis()->SetLabelSize(TextSize);
	h1_Erec_1pi_e_weight->GetYaxis()->SetNdivisions(Ndivisions);

	h1_Erec_1pi_e_weight->SetLineColor(46);

	h1_Erec_1pi_e_weight->GetYaxis()->SetTitle("Weighted Events / GeV");
	h1_Erec_1pi_e_weight->SetMarkerStyle(20);

	h1_Erec_1pi_e->SetMarkerStyle(20);
	h1_Erec_2pi_e_weight->SetMarkerStyle(20);
	h1_Erec_3pi_e_weight->SetMarkerStyle(20);
	h1_Erec_1pi_e_weight->SetMarkerSize(1.2);
	h1_Erec_1pi_e->SetMarkerSize(1.2);
	h1_Erec_2pi_e_weight->SetMarkerSize(1.2);
	h1_Erec_3pi_e_weight->SetMarkerSize(1.2);
	h1_Erec_1pi_e_weight->SetMarkerColor(46);
	h1_Erec_1pi_e->SetMarkerColor(38);
	h1_Erec_2pi_e_weight->SetMarkerColor(42);
	h1_Erec_3pi_e_weight->SetMarkerColor(28);

	h1_Erec_1pi_e_weight->GetYaxis()->SetRangeUser(0.,1.1*h1_Erec_1pi_e_weight->GetMaximum(););
	h1_Erec_1pi_e_weight->Draw("e");
	h1_Erec_1pi_e->Draw("e same");

	h1_Erec_2pi_e_weight->SetLineColor(42);
	h1_Erec_2pi_e_weight->Draw("P e Same");

	h1_Erec_3pi_e_weight->SetLineColor(28);
	h1_Erec_3pi_e_weight->Draw("P e Same");

	//h1_Erec_4pi_e_weight->SetLineColor(6);

	lat1->SetTextSize(TextSize);
	lat1->SetTextFont(FontStyle);
	
	lat1->SetTextColor(38);
	lat1->DrawLatex(0.58, 0.85, "Detected 1 #pi^{#pm} or #gamma");
	lat1->SetTextSize(TextSize);
	lat1->SetTextColor(46);
	lat1->DrawLatex(0.58, 0.78, "Undetected 1#pi^{#pm}/#gamma (-)");
	lat1->SetTextColor(42);
	lat1->DrawLatex(0.58, 0.71, "Undetected 2#pi^{#pm}/#gamma (+)");
	lat1->SetTextColor(28);
	lat1->DrawLatex(0.58, 0.64, "Undetected 3#pi^{#pm}/#gamma (+)");
	lat1->SetTextColor(6);

	lat1->SetTextColor(1);
	lat1->SetTextSize(TextSize);
	lat1->SetTextFont(FontStyle);
	lat1->SetTextSize(TextSize);
	lat1->DrawLatex(0.17, 0.85, "^{56}Fe");

	c1->SaveAs("../../myPlots/pdf/NoxBCut/"+version+target+"/FigExtData11_EQE_Inclusive_Subtractions_"+target+".pdf");

	// -------------------------------------------------------------------------------------------

	TCanvas *c2 = new TCanvas("c2","",1024,768);
	c2->cd();

	h1_E_rec->SetAxisRange(xmin, xmax, "X"); 
	h1_E_rec->GetXaxis()->SetTitle("E^{QE} [GeV]");
	h1_E_rec->GetYaxis()->SetTitle("Weighted Events / GeV");
	h1_E_rec->GetXaxis()->SetTitleOffset(1.2);
	h1_E_rec->SetMinimum(0);
	h1_E_rec->UseCurrentStyle();
	h1_E_rec->SetMarkerStyle(20);
	h1_E_rec->SetMarkerSize(1.2);
	h1_E_rec->SetMarkerColor(38);
	//h1_E_rec->SetMinimum(3*h1_Erec_3pi_e_weight->GetMinimum());

	h1_E_rec->GetXaxis()->CenterTitle();
	h1_E_rec->GetXaxis()->SetTitleFont(FontStyle);
	h1_E_rec->GetXaxis()->SetLabelFont(FontStyle);
	h1_E_rec->GetXaxis()->SetTitleSize(TextSize);
	h1_E_rec->GetXaxis()->SetLabelSize(TextSize);

	h1_E_rec->GetYaxis()->CenterTitle();
	h1_E_rec->GetYaxis()->SetTitleFont(FontStyle);
	h1_E_rec->GetYaxis()->SetLabelFont(FontStyle);
	h1_E_rec->GetYaxis()->SetTitleSize(TextSize);
	h1_E_rec->GetYaxis()->SetLabelSize(TextSize);

	//h1_E_rec->SetAxisRange(3*h1_Erec_3pi_e_weight->GetMinimum(),1.05*h1_E_rec->GetMaximum(),"Y");
	h1_E_rec->Draw("e");
	
	h1_Erec_0pi_e->SetAxisRange(xmin, xmax, "X"); 
	h1_Erec_0pi_e->SetLineColor(46);
	h1_Erec_0pi_e->GetXaxis()->SetTitleOffset(1.2);
	h1_Erec_0pi_e->SetMarkerStyle(20);
	h1_Erec_0pi_e->SetMarkerSize(1.2);
	h1_Erec_0pi_e->SetMarkerColor(46);
	h1_Erec_0pi_e->Draw("e same");

	h1_Erec_0pi_sub_e->SetLineColor(8);
	h1_Erec_0pi_sub_e->SetMarkerStyle(20);
	h1_Erec_0pi_sub_e->SetMarkerSize(1.2);
	h1_Erec_0pi_sub_e->SetMarkerColor(8);

	h1_Erec_0pi_sub_e->Draw("e Same");

	lat1->SetTextColor(38);
	lat1->DrawLatex(0.35, 0.73, "No cuts");

	lat1->SetTextColor(46);
	lat1->DrawLatex(0.3, 0.44, "No detected #pi^{#pm}/#gamma");
	 
	lat1->SetTextColor(8);
	lat1->DrawLatex(0.26, 0.24, "Subtract undetected #pi^{#pm}/#gamma");
		
	lat1->SetTextColor(1);
	lat1->SetTextSize(0.09);
	lat1->SetTextSize(TextSize);
	lat1->DrawLatex(0.17, 0.85, "^{56}Fe");

	c2->SaveAs("../../myPlots/pdf/NoxBCut/"+version+target+"/FigExtData11_EQE_Exclusive_Subtractions_"+target+".pdf");
 
}
