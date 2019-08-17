#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>

#include <iostream>
#include <vector>

using namespace std;

#include  "./Secondary_Code/CenterAxisTitle.cpp"
#include "./Secondary_Code/SetOffsetAndSize.cpp"
#include "./Secondary_Code/ToString.cpp"

void OverlayPlots() {


	SetOffsetAndSize();

	int FontStyle = 132;
	
	TString version = "v3_0_6/";

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> LabelsOfSamples; std::vector<TString> E;
	std::vector<TString> LabelE; std::vector<int> Colors; std::vector<TString> FSIModel;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots;

//	xBCut.push_back("xBCut");
//	xBCut.push_back("NoxBCut");

//	nucleus.push_back("3He"); LabelsOfSamples.push_back("^{3}He");
//	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He");
	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C");
//	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV");
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV");
//	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");
 
	Colors.push_back(kBlack); /*Colors.push_back(kBlue);*/ Colors.push_back(kRed); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);

	FSIModel.push_back("Data"); FSILabel.push_back("Data");
//	FSIModel.push_back("hN2018");FSILabel.push_back("hN2018");
//	FSIModel.push_back("hA2018_RadCorr"); FSILabel.push_back("Rad RFG");
//	FSIModel.push_back("hA2018_RadCorr_G4"); FSILabel.push_back("Rad G4");
//	FSIModel.push_back("hA2018_RadCorr_CFGM"); FSILabel.push_back("Rad CFG");
//	FSIModel.push_back("hA2018_RadCorr_LFGM"); FSILabel.push_back("Rad LFG");
//	FSIModel.push_back("hA2018_RadCorr_Benhar"); FSILabel.push_back("Rad + Benhar SF");
//	FSIModel.push_back("hA2018_RadCorr_pi0decay"); FSILabel.push_back("Rad pi0 decay");
	FSIModel.push_back("hA2018_Truth_NoRadCorr"); FSILabel.push_back("No Rad RFG");
//	FSIModel.push_back("hA2018_Truth_NoRadCorr"); FSILabel.push_back("Genie v3");
//	FSIModel.push_back("hA2015_Truth_NoRadCorr"); FSILabel.push_back("Genie v2");
//	FSIModel.push_back("hA2018_NoRadCorr_CFGM"); FSILabel.push_back("No Rad CFG");
//	FSIModel.push_back("hA2018_NoRadCorr_LFGM"); FSILabel.push_back("No Rad LFG");
//	FSIModel.push_back("hA2018_NoRadCorr_Benhar"); FSILabel.push_back("No Rad + Benhar SF");
//	FSIModel.push_back("hA2018_NoRadCorr"); FSILabel.push_back("hA2018");
//	FSIModel.push_back("hN2018_NoRadCorr"); FSILabel.push_back("hN2018");

//	NameOfPlots.push_back("ThetaAngleRecoNuBeamPlot"); 
//	NameOfPlots.push_back("CosThetaAngleRecoNuBeamPlot");
//	NameOfPlots.push_back("xB");
//	NameOfPlots.push_back("nu"); 
//	NameOfPlots.push_back("phi_ElectronOut"); 
//	NameOfPlots.push_back("Costheta_ElectronOut");
//	NameOfPlots.push_back("EePrime"); 
//	NameOfPlots.push_back("DeltaPhiTPlot"); 
//	NameOfPlots.push_back("DeltaAlphaTPlot"); 
//	NameOfPlots.push_back("phi_finalStateNucleonPlot");
//	NameOfPlots.push_back("costheta_finalStateNucleon"); 
//	NameOfPlots.push_back("Ep"); 
//	NameOfPlots.push_back("nu"); 
//	NameOfPlots.push_back("W");
//	NameOfPlots.push_back("Q2Inclusive");
//	NameOfPlots.push_back("EePrimeInclusive");
//	NameOfPlots.push_back("ElectronCosThetaInclusive"); 
//	NameOfPlots.push_back("ElectronPhiInclusive");
//	NameOfPlots.push_back("PionMultiPlot"); 
//	NameOfPlots.push_back("PionMultiQEPlot");
//	NameOfPlots.push_back("EcalReso");
//	NameOfPlots.push_back("EQEReso");

xBCut.push_back("NoxBCut");
	NameOfPlots.push_back("Q2");
	NameOfPlots.push_back("PionMultiPlot");  
	NameOfPlots.push_back("MissMomentum"); 
	NameOfPlots.push_back("epRecoEnergy_slice_0"); 
	NameOfPlots.push_back("eRecoEnergy_slice_0"); 
	NameOfPlots.push_back("epRecoEnergy_slice_1"); 
	NameOfPlots.push_back("eRecoEnergy_slice_1");
	NameOfPlots.push_back("epRecoEnergy_slice_2"); 
	NameOfPlots.push_back("eRecoEnergy_slice_2");
	NameOfPlots.push_back("epRecoEnergy_slice_3"); 
	NameOfPlots.push_back("eRecoEnergy_slice_3");

/*xBCut.push_back("xBCut");
	NameOfPlots.push_back("MissMomentum"); 
	NameOfPlots.push_back("epRecoEnergy_slice_0"); 
	NameOfPlots.push_back("eRecoEnergy_slice_0"); 
	NameOfPlots.push_back("epRecoEnergy_slice_1"); 
	NameOfPlots.push_back("eRecoEnergy_slice_1");
	NameOfPlots.push_back("epRecoEnergy_slice_2"); 
	NameOfPlots.push_back("eRecoEnergy_slice_2");
	NameOfPlots.push_back("epRecoEnergy_slice_3"); 
	NameOfPlots.push_back("eRecoEnergy_slice_3");*/

	std::vector<TCanvas*> PlotCanvas;

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> Plots_Clones;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	TString WhatModelsAreIncluded = "";
	for (int LoopOverFSIModels = 0 ; LoopOverFSIModels < NFSIModels ; LoopOverFSIModels ++) { WhatModelsAreIncluded += "_"+FSIModel[LoopOverFSIModels]; };

	TString RecoCalorimetry = "(e,e'p)";
	TString FSI = "FSI";

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					PlotCanvas.push_back(new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 205,34,1024,768));

					TPad* pad1 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],0,0.2,1,1, 21); 
					pad1->SetFillColor(kWhite); pad1->Draw();
					TPad* pad2 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],0,0.0,1,0.2,22); 
					pad2->SetFillColor(kWhite); pad2->Draw(); 
					pad1->cd();

					Plots.clear();

					TLegend* leg = new TLegend(0.15,0.92,0.9,0.99);
					leg->SetNColumns(3);

					double max = -99.;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );
						Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
						Plots[WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						//Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(0);
						//Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(0);
						Plots[WhichFSIModel]->GetYaxis()->SetTickSize(0.02);
						Plots[WhichFSIModel]->SetLineWidth(3);
						CenterAxisTitle(Plots[WhichFSIModel]);

//						if (FSIModel[WhichFSIModel] == "Data") { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "lep");}
//						else { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l"); }
leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l");

						double NBins = Plots[WhichFSIModel]->GetNbinsX(); 
				
						for (int i = 1; i <= NBins; i++) { 
					
							double content = Plots[WhichFSIModel]->GetBinContent(i);
							double error = Plots[WhichFSIModel]->GetBinError(i);
							double width = Plots[WhichFSIModel]->GetBinWidth(i);
							double newcontent = content / width;
							double newerror = error / width;				
							Plots[WhichFSIModel]->SetBinContent(i,newcontent);
							Plots[WhichFSIModel]->SetBinError(i,newerror);

						}

						double ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral();
						Plots[WhichFSIModel]->Scale(ScalingFactor);
if (NameOfPlots[WhichPlot] == "MissMomentum" || NameOfPlots[WhichPlot] == "Q2") { for (int i = 0; i < 2;i++) { Plots[WhichFSIModel]->Rebin(); } }

						double localmax = Plots[WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }
						Plots[0]->GetYaxis()->SetRangeUser(0,1.1*max);

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

if (NameOfPlots[WhichPlot] == "Q2" && E[WhichEnergy] == "1_161") { Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.,1.); }

//						if (FSIModel[WhichFSIModel] != "Data") { Plots[WhichFSIModel]->Draw("e hist same"); }
//						else { Plots[WhichFSIModel]->SetMarkerStyle(20); Plots[WhichFSIModel]->Draw("e same"); }
Plots[WhichFSIModel]->Draw("hist same");

					} // End of the loop over the FSI Models 

					leg->SetBorderSize(0);
					leg->SetTextFont(FontStyle);
					leg->SetTextSize(0.06);
//					leg->SetTextSize(0.07);			
					leg->Draw();

					TLatex latex;
					latex.SetTextFont(FontStyle);
					latex.SetTextSize(0.07);
/*					latex.DrawLatexNDC(0.22,0.93,LabelsOfSamples[WhichNucleus]+RecoCalorimetry+LabelE[WhichEnergy]);*/

					pad2->cd(); pad2->SetGrid();
					Plots_Clones.clear();
		
					for (int WhichFSIModel = 0 ; WhichFSIModel < NFSIModels ; WhichFSIModel ++ ) {
		
						Plots_Clones.push_back( (TH1D*)(Plots[WhichFSIModel]->Clone())) ;
						Plots_Clones[WhichFSIModel]->GetXaxis()->SetTitleSize(0.0);
						Plots_Clones[WhichFSIModel]->GetXaxis()->SetLabelSize(0.0);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetRangeUser(0,2);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetLabelSize(0.12);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleSize(0.2);

						// Residual

						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.2);
						Plots_Clones[WhichFSIModel]->Add(Plots[0],-1);
						Plots_Clones[WhichFSIModel]->Divide(Plots[0]);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitle("#frac{MC-Data}{Data}");			
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetRangeUser(-1.1,1.1);

						Plots_Clones[WhichFSIModel]->GetYaxis()->SetNdivisions(5);

//						if (WhichFSIModel == 0) { Plots_Clones[WhichFSIModel]->Draw("hist same"); }
//						else { Plots_Clones[WhichFSIModel]->Draw("e1 hist same"); }

Plots_Clones[WhichFSIModel]->Draw("hist same");

					}

					PlotCanvas[PlotCanvas.size()-1]->SaveAs("myPlots/eps/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"+E[WhichEnergy]+"/"+nucleus[WhichNucleus]+"_" 
						+E[WhichEnergy]+"_" +NameOfPlots[WhichPlot]+WhatModelsAreIncluded+".eps"); 

					PlotCanvas[PlotCanvas.size()-1]->SaveAs("myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"+E[WhichEnergy]+"/"+nucleus[WhichNucleus]+"_" 
						+E[WhichEnergy]+"_" +NameOfPlots[WhichPlot]+WhatModelsAreIncluded+".pdf");

				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
