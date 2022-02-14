void Correction_Panel()
{
//=========Macro generated from canvas: CommonPanel_NoxBCut/CommonPanel_NoxBCut
//=========  (Sun Feb 13 17:04:13 2022) by ROOT version 6.24/02
   TCanvas *CommonPanel_NoxBCut = new TCanvas("CommonPanel_NoxBCut", "CommonPanel_NoxBCut",88,64,1832,1016);
   gStyle->SetOptStat(0);
   CommonPanel_NoxBCut->Range(0,0,1,1);
   CommonPanel_NoxBCut->SetFillColor(0);
   CommonPanel_NoxBCut->SetBorderMode(0);
   CommonPanel_NoxBCut->SetBorderSize(2);
   CommonPanel_NoxBCut->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: One_161
   TPad *One_161 = new TPad("One_161", "",0.02,0.1,0.36,0.44);
   One_161->Draw();
   One_161->cd();
   One_161->Range(0.3613388,0.2365854,1.23748,1.7);
   One_161->SetFillColor(0);
   One_161->SetBorderMode(0);
   One_161->SetBorderSize(2);
   One_161->SetLeftMargin(0.15);
   One_161->SetRightMargin(0);
   One_161->SetTopMargin(0);
   One_161->SetBottomMargin(0.18);
   One_161->SetFrameBorderMode(0);
   One_161->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1[38] = {
   0.42,
   0.46,
   0.5,
   0.54,
   0.58,
   0.62,
   0.66,
   0.7,
   0.74,
   0.78,
   0.82,
   0.86,
   0.9,
   0.94,
   0.98,
   1.02,
   1.06,
   1.09,
   1.11,
   1.13,
   1.15,
   1.17,
   1.19,
   1.21,
   1.23,
   1.25,
   1.27,
   1.29,
   1.31,
   1.33,
   1.35,
   1.37,
   1.39,
   1.41,
   1.43,
   1.45,
   1.47,
   1.49};
   Double_t Graph0_fy1[38] = {
   0,
   1.010248,
   0.9450611,
   0.9566903,
   0.9362088,
   0.9349875,
   0.9318715,
   0.9551971,
   0.9413819,
   0.9444444,
   0.93405,
   0.9342322,
   0.9272029,
   0.879539,
   0.8275254,
   0.7491577,
   0.6874865,
   0.6944359,
   0.7102881,
   0.7870498,
   0.9367996,
   1.077697,
   1.161891,
   1.223094,
   1.202515,
   1.209453,
   2.029143,
   1.691029,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   TGraph *graph = new TGraph(38,Graph0_fx1,Graph0_fy1);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph01 = new TH1F("Graph_Graph01","",100,0.313,1.597);
   Graph_Graph01->SetMinimum(0.5);
   Graph_Graph01->SetMaximum(1.7);
   Graph_Graph01->SetDirectory(0);
   Graph_Graph01->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph01->SetLineColor(ci);
   Graph_Graph01->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph01->GetXaxis()->SetRange(15,72);
   Graph_Graph01->GetXaxis()->CenterTitle(true);
   Graph_Graph01->GetXaxis()->SetNdivisions(6);
   Graph_Graph01->GetXaxis()->SetLabelFont(132);
   Graph_Graph01->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph01->GetXaxis()->SetTitleSize(0);
   Graph_Graph01->GetXaxis()->SetTickLength(0.02);
   Graph_Graph01->GetXaxis()->SetTitleOffset(1);
   Graph_Graph01->GetXaxis()->SetTitleFont(132);
   Graph_Graph01->GetYaxis()->SetTitle("Radiation Correction");
   Graph_Graph01->GetYaxis()->CenterTitle(true);
   Graph_Graph01->GetYaxis()->SetNdivisions(8);
   Graph_Graph01->GetYaxis()->SetLabelFont(132);
   Graph_Graph01->GetYaxis()->SetLabelSize(0.08);
   Graph_Graph01->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph01->GetYaxis()->SetTickLength(0.02);
   Graph_Graph01->GetYaxis()->SetTitleOffset(1);
   Graph_Graph01->GetYaxis()->SetTitleFont(132);
   Graph_Graph01->GetZaxis()->SetLabelFont(42);
   Graph_Graph01->GetZaxis()->SetTitleOffset(1);
   Graph_Graph01->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph01);
   
   graph->Draw("ap");
   TLine *line = new TLine(1.161,0.5,1.161,1.7);
   line->SetLineStyle(2);
   line->Draw();
   TLatex *   tex = new TLatex(0.8,0.9,"(g)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.09);
   tex->SetLineWidth(2);
   tex->Draw();
   One_161->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: Two_261
   TPad *Two_261 = new TPad("Two_261", "Two_261",0.36,0.1,0.68,0.44);
   Two_261->Draw();
   Two_261->cd();
   Two_261->Range(0.6558,0.2365854,2.42646,1.7);
   Two_261->SetFillColor(0);
   Two_261->SetBorderMode(0);
   Two_261->SetBorderSize(2);
   Two_261->SetLeftMargin(0);
   Two_261->SetRightMargin(0);
   Two_261->SetTopMargin(0);
   Two_261->SetBottomMargin(0.18);
   Two_261->SetFrameBorderMode(0);
   Two_261->SetFrameBorderMode(0);
   
   Double_t Graph0_fx2[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   2.505,
   2.535,
   2.565,
   2.595,
   2.625,
   2.655,
   2.685,
   2.715,
   2.745,
   2.775,
   2.805,
   2.835,
   2.865,
   2.895,
   2.925,
   2.955,
   2.985};
   Double_t Graph0_fy2[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   1.055256,
   0.8996129,
   0.8845184,
   0.8799289,
   0.8930914,
   0.8924878,
   0.8978884,
   0.8906136,
   0.8851391,
   0.8895888,
   0.8807489,
   0.8741892,
   0.8652025,
   0.8824504,
   0.8661349,
   0.871832,
   0.8784866,
   0.8899812,
   0.7804046,
   0.6428642,
   0.5885779,
   0.5839221,
   0.8528351,
   1.254687,
   1.388233,
   1.387228,
   0.5587223,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(54,Graph0_fx2,Graph0_fy2);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph02 = new TH1F("Graph_Graph02","",100,0,3.279);
   Graph_Graph02->SetMinimum(0.5);
   Graph_Graph02->SetMaximum(1.7);
   Graph_Graph02->SetDirectory(0);
   Graph_Graph02->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph02->SetLineColor(ci);
   Graph_Graph02->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph02->GetXaxis()->SetRange(21,74);
   Graph_Graph02->GetXaxis()->CenterTitle(true);
   Graph_Graph02->GetXaxis()->SetNdivisions(6);
   Graph_Graph02->GetXaxis()->SetLabelFont(132);
   Graph_Graph02->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph02->GetXaxis()->SetTitleSize(0);
   Graph_Graph02->GetXaxis()->SetTickLength(0.02);
   Graph_Graph02->GetXaxis()->SetTitleOffset(1);
   Graph_Graph02->GetXaxis()->SetTitleFont(132);
   Graph_Graph02->GetYaxis()->SetTitle("Radiation Correction");
   Graph_Graph02->GetYaxis()->CenterTitle(true);
   Graph_Graph02->GetYaxis()->SetNdivisions(8);
   Graph_Graph02->GetYaxis()->SetLabelFont(132);
   Graph_Graph02->GetYaxis()->SetLabelSize(0.08);
   Graph_Graph02->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph02->GetYaxis()->SetTickLength(0.02);
   Graph_Graph02->GetYaxis()->SetTitleOffset(1);
   Graph_Graph02->GetYaxis()->SetTitleFont(132);
   Graph_Graph02->GetZaxis()->SetLabelFont(42);
   Graph_Graph02->GetZaxis()->SetTitleOffset(1);
   Graph_Graph02->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph02);
   
   graph->Draw("ap");
   
   Double_t Graph1_fx3[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   2.505,
   2.535,
   2.565,
   2.595,
   2.625,
   2.655,
   2.685,
   2.715,
   2.745,
   2.775,
   2.805,
   2.835,
   2.865,
   2.895,
   2.925,
   2.955,
   2.985};
   Double_t Graph1_fy3[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0.9193358,
   0.8628011,
   0.885365,
   0.8949566,
   0.8826929,
   0.8713615,
   0.8969707,
   0.8837793,
   0.8741324,
   0.8669068,
   0.8727566,
   0.8654152,
   0.8679702,
   0.8885612,
   0.8689399,
   0.8849983,
   0.9073521,
   0.9311382,
   0.8474227,
   0.741705,
   0.6857746,
   0.664357,
   0.9467782,
   1.284986,
   1.491228,
   1.637389,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(54,Graph1_fx3,Graph1_fy3);
   graph->SetName("Graph1");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#cc66cc");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph13 = new TH1F("Graph_Graph13","",100,0,3.279);
   Graph_Graph13->SetMinimum(0.5);
   Graph_Graph13->SetMaximum(1.7);
   Graph_Graph13->SetDirectory(0);
   Graph_Graph13->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph13->SetLineColor(ci);
   Graph_Graph13->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph13->GetXaxis()->SetRange(21,74);
   Graph_Graph13->GetXaxis()->CenterTitle(true);
   Graph_Graph13->GetXaxis()->SetNdivisions(6);
   Graph_Graph13->GetXaxis()->SetLabelFont(132);
   Graph_Graph13->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph13->GetXaxis()->SetTitleSize(0);
   Graph_Graph13->GetXaxis()->SetTickLength(0.02);
   Graph_Graph13->GetXaxis()->SetTitleOffset(1);
   Graph_Graph13->GetXaxis()->SetTitleFont(132);
   Graph_Graph13->GetYaxis()->SetTitle("Radiation Correction");
   Graph_Graph13->GetYaxis()->CenterTitle(true);
   Graph_Graph13->GetYaxis()->SetNdivisions(8);
   Graph_Graph13->GetYaxis()->SetLabelFont(132);
   Graph_Graph13->GetYaxis()->SetLabelSize(0.08);
   Graph_Graph13->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph13->GetYaxis()->SetTickLength(0.02);
   Graph_Graph13->GetYaxis()->SetTitleOffset(1);
   Graph_Graph13->GetYaxis()->SetTitleFont(132);
   Graph_Graph13->GetZaxis()->SetLabelFont(42);
   Graph_Graph13->GetZaxis()->SetTitleOffset(1);
   Graph_Graph13->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph13);
   
   graph->Draw("p");
   
   Double_t Graph2_fx4[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   2.505,
   2.535,
   2.565,
   2.595,
   2.625,
   2.655,
   2.685,
   2.715,
   2.745,
   2.775,
   2.805,
   2.835,
   2.865,
   2.895,
   2.925,
   2.955,
   2.985};
   Double_t Graph2_fy4[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0.9466655,
   0.8977445,
   0.8545473,
   0.8789924,
   0.8852509,
   0.8964149,
   0.8908612,
   0.8765246,
   0.8896123,
   0.8856629,
   0.8828254,
   0.857903,
   0.8393462,
   0.873602,
   0.8084213,
   0.798274,
   0.7758828,
   0.7523847,
   0.5693482,
   0.4318349,
   0.3963113,
   0.4055719,
   0.6080127,
   1.164875,
   1.361257,
   1.396638,
   2.449875,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(54,Graph2_fx4,Graph2_fy4);
   graph->SetName("Graph2");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#66cc66");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph24 = new TH1F("Graph_Graph24","",100,0,3.279);
   Graph_Graph24->SetMinimum(0.5);
   Graph_Graph24->SetMaximum(1.7);
   Graph_Graph24->SetDirectory(0);
   Graph_Graph24->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph24->SetLineColor(ci);
   Graph_Graph24->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph24->GetXaxis()->SetRange(21,74);
   Graph_Graph24->GetXaxis()->CenterTitle(true);
   Graph_Graph24->GetXaxis()->SetNdivisions(6);
   Graph_Graph24->GetXaxis()->SetLabelFont(132);
   Graph_Graph24->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph24->GetXaxis()->SetTitleSize(0);
   Graph_Graph24->GetXaxis()->SetTickLength(0.02);
   Graph_Graph24->GetXaxis()->SetTitleOffset(1);
   Graph_Graph24->GetXaxis()->SetTitleFont(132);
   Graph_Graph24->GetYaxis()->SetTitle("Radiation Correction");
   Graph_Graph24->GetYaxis()->CenterTitle(true);
   Graph_Graph24->GetYaxis()->SetNdivisions(8);
   Graph_Graph24->GetYaxis()->SetLabelFont(132);
   Graph_Graph24->GetYaxis()->SetLabelSize(0.08);
   Graph_Graph24->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph24->GetYaxis()->SetTickLength(0.02);
   Graph_Graph24->GetYaxis()->SetTitleOffset(1);
   Graph_Graph24->GetYaxis()->SetTitleFont(132);
   Graph_Graph24->GetZaxis()->SetLabelFont(42);
   Graph_Graph24->GetZaxis()->SetTitleOffset(1);
   Graph_Graph24->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph24);
   
   graph->Draw("p");
   line = new TLine(2.261,0.5,2.261,1.7);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.8,0.9,"(h)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.09);
   tex->SetLineWidth(2);
   tex->Draw();
   Two_261->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: Four_461
   TPad *Four_461 = new TPad("Four_461", "Four_461",0.68,0.1,1,0.44);
   Four_461->Draw();
   Four_461->cd();
   Four_461->Range(1.489725,0.2365854,4.765741,1.7);
   Four_461->SetFillColor(0);
   Four_461->SetBorderMode(0);
   Four_461->SetBorderSize(2);
   Four_461->SetLeftMargin(0);
   Four_461->SetRightMargin(0.04);
   Four_461->SetTopMargin(0);
   Four_461->SetBottomMargin(0.18);
   Four_461->SetFrameBorderMode(0);
   Four_461->SetFrameBorderMode(0);
   
   Double_t Graph0_fx5[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   4.925,
   4.975,
   5.025};
   Double_t Graph0_fy5[38] = {
   0,
   0,
   0,
   0,
   0,
   4.936194,
   1.010585,
   0.9371774,
   0.9272884,
   0.8802602,
   0.9044852,
   0.8963439,
   0.9013297,
   0.8935865,
   0.9239848,
   0.9016797,
   0.9139214,
   0.8884134,
   0.8950043,
   0.9569537,
   1.004665,
   1.033289,
   1.064011,
   0.9000366,
   0.5604075,
   0.9791487,
   1.347526,
   1.512727,
   0.8839779,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(38,Graph0_fx5,Graph0_fy5);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph05 = new TH1F("Graph_Graph05","",100,0,5.5175);
   Graph_Graph05->SetMinimum(0.5);
   Graph_Graph05->SetMaximum(1.7);
   Graph_Graph05->SetDirectory(0);
   Graph_Graph05->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph05->SetLineColor(ci);
   Graph_Graph05->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph05->GetXaxis()->SetRange(28,84);
   Graph_Graph05->GetXaxis()->CenterTitle(true);
   Graph_Graph05->GetXaxis()->SetNdivisions(6);
   Graph_Graph05->GetXaxis()->SetLabelFont(132);
   Graph_Graph05->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph05->GetXaxis()->SetTitleSize(0);
   Graph_Graph05->GetXaxis()->SetTickLength(0.02);
   Graph_Graph05->GetXaxis()->SetTitleOffset(1);
   Graph_Graph05->GetXaxis()->SetTitleFont(132);
   Graph_Graph05->GetYaxis()->SetTitle("Radiation Correction");
   Graph_Graph05->GetYaxis()->CenterTitle(true);
   Graph_Graph05->GetYaxis()->SetNdivisions(8);
   Graph_Graph05->GetYaxis()->SetLabelFont(132);
   Graph_Graph05->GetYaxis()->SetLabelSize(0.08);
   Graph_Graph05->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph05->GetYaxis()->SetTickLength(0.02);
   Graph_Graph05->GetYaxis()->SetTitleOffset(1);
   Graph_Graph05->GetYaxis()->SetTitleFont(132);
   Graph_Graph05->GetZaxis()->SetLabelFont(42);
   Graph_Graph05->GetZaxis()->SetTitleOffset(1);
   Graph_Graph05->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph05);
   
   graph->Draw("ap");
   
   Double_t Graph1_fx6[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   4.925,
   4.975,
   5.025};
   Double_t Graph1_fy6[38] = {
   0,
   0,
   0,
   0,
   0,
   0.8407089,
   0.833643,
   0.9445971,
   0.9058475,
   0.9085406,
   0.9063746,
   0.8925519,
   0.9008977,
   0.9111861,
   0.9167367,
   0.9203456,
   0.8981504,
   0.9032482,
   0.8992437,
   0.9615844,
   1.038489,
   1.052523,
   1.062335,
   1.011841,
   0.6695448,
   1.03837,
   1.386461,
   1.587464,
   4.000619,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(38,Graph1_fx6,Graph1_fy6);
   graph->SetName("Graph1");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#cc66cc");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph16 = new TH1F("Graph_Graph16","",100,0,5.5175);
   Graph_Graph16->SetMinimum(0.5);
   Graph_Graph16->SetMaximum(1.7);
   Graph_Graph16->SetDirectory(0);
   Graph_Graph16->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph16->SetLineColor(ci);
   Graph_Graph16->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph16->GetXaxis()->SetRange(28,84);
   Graph_Graph16->GetXaxis()->CenterTitle(true);
   Graph_Graph16->GetXaxis()->SetNdivisions(6);
   Graph_Graph16->GetXaxis()->SetLabelFont(132);
   Graph_Graph16->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph16->GetXaxis()->SetTitleSize(0);
   Graph_Graph16->GetXaxis()->SetTickLength(0.02);
   Graph_Graph16->GetXaxis()->SetTitleOffset(1);
   Graph_Graph16->GetXaxis()->SetTitleFont(132);
   Graph_Graph16->GetYaxis()->SetTitle("Radiation Correction");
   Graph_Graph16->GetYaxis()->CenterTitle(true);
   Graph_Graph16->GetYaxis()->SetNdivisions(8);
   Graph_Graph16->GetYaxis()->SetLabelFont(132);
   Graph_Graph16->GetYaxis()->SetLabelSize(0.08);
   Graph_Graph16->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph16->GetYaxis()->SetTickLength(0.02);
   Graph_Graph16->GetYaxis()->SetTitleOffset(1);
   Graph_Graph16->GetYaxis()->SetTitleFont(132);
   Graph_Graph16->GetZaxis()->SetLabelFont(42);
   Graph_Graph16->GetZaxis()->SetTitleOffset(1);
   Graph_Graph16->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph16);
   
   graph->Draw("p");
   
   Double_t Graph2_fx7[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   4.925,
   4.975,
   5.025};
   Double_t Graph2_fy7[38] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0.8460809,
   0.9435963,
   0.8690503,
   0.8833891,
   0.882144,
   0.8677245,
   0.9025255,
   0.9172765,
   0.902195,
   0.9104946,
   0.9285209,
   0.8897215,
   0.8853352,
   0.9446481,
   0.9981915,
   1.005326,
   0.9325224,
   0.7580866,
   0.3667641,
   0.8182022,
   1.299435,
   1.470901,
   1.43247,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(38,Graph2_fx7,Graph2_fy7);
   graph->SetName("Graph2");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#66cc66");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph27 = new TH1F("Graph_Graph27","",100,0,5.5175);
   Graph_Graph27->SetMinimum(0.5);
   Graph_Graph27->SetMaximum(1.7);
   Graph_Graph27->SetDirectory(0);
   Graph_Graph27->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph27->SetLineColor(ci);
   Graph_Graph27->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph27->GetXaxis()->SetRange(28,84);
   Graph_Graph27->GetXaxis()->CenterTitle(true);
   Graph_Graph27->GetXaxis()->SetNdivisions(6);
   Graph_Graph27->GetXaxis()->SetLabelFont(132);
   Graph_Graph27->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph27->GetXaxis()->SetTitleSize(0);
   Graph_Graph27->GetXaxis()->SetTickLength(0.02);
   Graph_Graph27->GetXaxis()->SetTitleOffset(1);
   Graph_Graph27->GetXaxis()->SetTitleFont(132);
   Graph_Graph27->GetYaxis()->SetTitle("Radiation Correction");
   Graph_Graph27->GetYaxis()->CenterTitle(true);
   Graph_Graph27->GetYaxis()->SetNdivisions(8);
   Graph_Graph27->GetYaxis()->SetLabelFont(132);
   Graph_Graph27->GetYaxis()->SetLabelSize(0.08);
   Graph_Graph27->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph27->GetYaxis()->SetTickLength(0.02);
   Graph_Graph27->GetYaxis()->SetTitleOffset(1);
   Graph_Graph27->GetYaxis()->SetTitleFont(132);
   Graph_Graph27->GetZaxis()->SetLabelFont(42);
   Graph_Graph27->GetZaxis()->SetTitleOffset(1);
   Graph_Graph27->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph27);
   
   graph->Draw("p");
   line = new TLine(4.461,0.5,4.461,1.7);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.8,0.9,"(i)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.09);
   tex->SetLineWidth(2);
   tex->Draw();
   Four_461->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: PadOne_161
   TPad *PadOne_161 = new TPad("PadOne_161", "",0.02,0.44,0.36,0.72);
   PadOne_161->Draw();
   PadOne_161->cd();
   PadOne_161->Range(0.3558241,-0.5,1.24333,17);
   PadOne_161->SetFillColor(0);
   PadOne_161->SetBorderMode(0);
   PadOne_161->SetBorderSize(2);
   PadOne_161->SetLeftMargin(0.15);
   PadOne_161->SetRightMargin(0);
   PadOne_161->SetTopMargin(0);
   PadOne_161->SetBottomMargin(0);
   PadOne_161->SetFrameBorderMode(0);
   PadOne_161->SetFrameBorderMode(0);
   
   Double_t Graph0_fx8[38] = {
   0.42,
   0.46,
   0.5,
   0.54,
   0.58,
   0.62,
   0.66,
   0.7,
   0.74,
   0.78,
   0.82,
   0.86,
   0.9,
   0.94,
   0.98,
   1.02,
   1.06,
   1.09,
   1.11,
   1.13,
   1.15,
   1.17,
   1.19,
   1.21,
   1.23,
   1.25,
   1.27,
   6.908798e-310,
   0,
   4.656784e-310,
   4.656789e-310,
   6.908798e-310,
   0,
   0,
   4.656784e-310,
   6.952837e-310,
   6.908802e-310,
   6.908783e-310};
   Double_t Graph0_fy8[38] = {
   0,
   57.73503,
   3.969716,
   2.65411,
   2.869703,
   3.144954,
   3.747712,
   2.752961,
   2.586166,
   1.518242,
   1.656629,
   0.3887273,
   0.621417,
   0.001119933,
   1.182577,
   0.891005,
   0.9111645,
   3.411766,
   3.411766,
   3.411766,
   3.411766,
   3.411766,
   3.411766,
   3.411766,
   3.411766,
   13.80455,
   24.416,
   6.908783e-310,
   2e-12,
   1.877449e-322,
   0,
   4.656789e-310,
   0,
   7.364134e-273,
   4.656789e-310,
   4.656789e-310,
   4.940656e-324,
   4.656789e-310};
   graph = new TGraph(38,Graph0_fx8,Graph0_fy8);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph08 = new TH1F("Graph_Graph08","",100,0,1.397);
   Graph_Graph08->SetMinimum(-0.5);
   Graph_Graph08->SetMaximum(17);
   Graph_Graph08->SetDirectory(0);
   Graph_Graph08->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph08->SetLineColor(ci);
   Graph_Graph08->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph08->GetXaxis()->SetRange(36,89);
   Graph_Graph08->GetXaxis()->CenterTitle(true);
   Graph_Graph08->GetXaxis()->SetNdivisions(6);
   Graph_Graph08->GetXaxis()->SetLabelFont(132);
   Graph_Graph08->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph08->GetXaxis()->SetTitleSize(0);
   Graph_Graph08->GetXaxis()->SetTickLength(0.02);
   Graph_Graph08->GetXaxis()->SetTitleOffset(1);
   Graph_Graph08->GetXaxis()->SetTitleFont(132);
   Graph_Graph08->GetYaxis()->SetTitle("Acceptance Uncertainty [%]");
   Graph_Graph08->GetYaxis()->CenterTitle(true);
   Graph_Graph08->GetYaxis()->SetNdivisions(8);
   Graph_Graph08->GetYaxis()->SetLabelFont(132);
   Graph_Graph08->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph08->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph08->GetYaxis()->SetTickLength(0.02);
   Graph_Graph08->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph08->GetYaxis()->SetTitleFont(132);
   Graph_Graph08->GetZaxis()->SetLabelFont(42);
   Graph_Graph08->GetZaxis()->SetTitleOffset(1);
   Graph_Graph08->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph08);
   
   graph->Draw("ap");
   line = new TLine(1.166,-0.4,1.166,17);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.8,0.9,"(d)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
   PadOne_161->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: NewTwo_261
   TPad *NewTwo_261 = new TPad("NewTwo_261", "NewTwo_261",0.36,0.44,0.68,0.72);
   NewTwo_261->Draw();
   NewTwo_261->cd();
   NewTwo_261->Range(0.6534,-0.5,2.423025,17);
   NewTwo_261->SetFillColor(0);
   NewTwo_261->SetBorderMode(0);
   NewTwo_261->SetBorderSize(2);
   NewTwo_261->SetLeftMargin(0);
   NewTwo_261->SetRightMargin(0);
   NewTwo_261->SetTopMargin(0);
   NewTwo_261->SetBottomMargin(0);
   NewTwo_261->SetFrameBorderMode(0);
   NewTwo_261->SetFrameBorderMode(0);
   
   Double_t Graph0_fx9[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   6.952837e-310,
   6.952837e-310,
   6.908802e-310,
   6.952837e-310,
   6.908802e-310,
   4.656789e-310,
   6.908798e-310,
   0,
   4.656784e-310,
   4.656789e-310,
   6.908798e-310,
   0,
   0,
   4.656784e-310,
   6.952837e-310,
   6.908802e-310,
   6.908783e-310};
   Double_t Graph0_fy9[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   57.73503,
   7.116107,
   5.467084,
   4.518096,
   3.780981,
   2.31185,
   1.777874,
   1.127145,
   0.3512244,
   0.1875748,
   0.1356921,
   0.2855485,
   0.1247536,
   0.6516234,
   0.420176,
   0.2294521,
   0.286587,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   8.215423,
   0,
   0,
   1.72923e-322,
   4.656789e-310,
   6.908783e-310,
   -6.434914,
   2.385,
   2.766768e-322,
   4.656789e-310,
   6.952837e-310,
   2.717361e-322,
   6.952837e-310,
   6.908783e-310,
   2.48515e-321,
   -1,
   1,
   1,
   4.656789e-310,
   4.656789e-310};
   graph = new TGraph(54,Graph0_fx9,Graph0_fy9);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph09 = new TH1F("Graph_Graph09","",100,0,2.7225);
   Graph_Graph09->SetMinimum(-0.5);
   Graph_Graph09->SetMaximum(17);
   Graph_Graph09->SetDirectory(0);
   Graph_Graph09->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph09->SetLineColor(ci);
   Graph_Graph09->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph09->GetXaxis()->SetRange(25,89);
   Graph_Graph09->GetXaxis()->CenterTitle(true);
   Graph_Graph09->GetXaxis()->SetNdivisions(6);
   Graph_Graph09->GetXaxis()->SetLabelFont(132);
   Graph_Graph09->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph09->GetXaxis()->SetTitleSize(0);
   Graph_Graph09->GetXaxis()->SetTickLength(0.02);
   Graph_Graph09->GetXaxis()->SetTitleOffset(1);
   Graph_Graph09->GetXaxis()->SetTitleFont(132);
   Graph_Graph09->GetYaxis()->SetTitle("Acceptance Uncertainty [%]");
   Graph_Graph09->GetYaxis()->CenterTitle(true);
   Graph_Graph09->GetYaxis()->SetNdivisions(8);
   Graph_Graph09->GetYaxis()->SetLabelFont(132);
   Graph_Graph09->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph09->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph09->GetYaxis()->SetTickLength(0.02);
   Graph_Graph09->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph09->GetYaxis()->SetTitleFont(132);
   Graph_Graph09->GetZaxis()->SetLabelFont(42);
   Graph_Graph09->GetZaxis()->SetTitleOffset(1);
   Graph_Graph09->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph09);
   
   graph->Draw("ap");
   
   Double_t Graph1_fx10[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   6.952837e-310,
   6.952837e-310,
   6.908802e-310,
   6.952837e-310,
   6.908802e-310,
   4.656789e-310,
   6.908798e-310,
   0,
   4.656784e-310,
   4.656789e-310,
   6.908798e-310,
   0,
   0,
   4.656784e-310,
   6.952837e-310,
   6.908802e-310,
   6.908783e-310};
   Double_t Graph1_fy10[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   57.73503,
   10.27673,
   6.73464,
   4.302078,
   3.447274,
   2.548063,
   1.635999,
   1.435263,
   0.8543941,
   0.3626916,
   0.2899328,
   0.4796319,
   0.4617094,
   0.7853907,
   0.5246562,
   0.1255086,
   1.293683,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   7.661252,
   0,
   0,
   1.72923e-322,
   4.656789e-310,
   6.908783e-310,
   -10.22367,
   2.385,
   2.766768e-322,
   4.656789e-310,
   6.952837e-310,
   2.717361e-322,
   6.952837e-310,
   6.908783e-310,
   2.48515e-321,
   -1,
   1,
   1,
   4.656789e-310,
   4.656789e-310};
   graph = new TGraph(54,Graph1_fx10,Graph1_fy10);
   graph->SetName("Graph1");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#cc66cc");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph110 = new TH1F("Graph_Graph110","",100,0,2.7225);
   Graph_Graph110->SetMinimum(-0.5);
   Graph_Graph110->SetMaximum(17);
   Graph_Graph110->SetDirectory(0);
   Graph_Graph110->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph110->SetLineColor(ci);
   Graph_Graph110->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph110->GetXaxis()->SetRange(25,89);
   Graph_Graph110->GetXaxis()->CenterTitle(true);
   Graph_Graph110->GetXaxis()->SetNdivisions(6);
   Graph_Graph110->GetXaxis()->SetLabelFont(132);
   Graph_Graph110->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph110->GetXaxis()->SetTitleSize(0);
   Graph_Graph110->GetXaxis()->SetTickLength(0.02);
   Graph_Graph110->GetXaxis()->SetTitleOffset(1);
   Graph_Graph110->GetXaxis()->SetTitleFont(132);
   Graph_Graph110->GetYaxis()->SetTitle("Acceptance Uncertainty [%]");
   Graph_Graph110->GetYaxis()->CenterTitle(true);
   Graph_Graph110->GetYaxis()->SetNdivisions(8);
   Graph_Graph110->GetYaxis()->SetLabelFont(132);
   Graph_Graph110->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph110->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph110->GetYaxis()->SetTickLength(0.02);
   Graph_Graph110->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph110->GetYaxis()->SetTitleFont(132);
   Graph_Graph110->GetZaxis()->SetLabelFont(42);
   Graph_Graph110->GetZaxis()->SetTitleOffset(1);
   Graph_Graph110->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph110);
   
   graph->Draw("p");
   
   Double_t Graph2_fx11[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   6.952837e-310,
   6.952837e-310,
   6.908802e-310,
   6.952837e-310,
   6.908802e-310,
   4.656789e-310,
   6.908798e-310,
   0,
   4.656784e-310,
   4.656789e-310,
   6.908798e-310,
   0,
   0,
   4.656784e-310,
   6.952837e-310,
   6.908802e-310,
   6.908783e-310};
   Double_t Graph2_fy11[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   57.73503,
   10.57272,
   4.630694,
   4.464939,
   3.99381,
   2.274445,
   1.631184,
   0.2878288,
   0.245061,
   0.2958873,
   0.2808782,
   0.3227102,
   0.04039297,
   0.3575605,
   0.2914832,
   0.3802301,
   0.1544711,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   12.34424,
   0,
   0,
   1.72923e-322,
   4.656789e-310,
   6.908783e-310,
   -13.74519,
   2.385,
   2.766768e-322,
   4.656789e-310,
   6.952837e-310,
   2.717361e-322,
   6.952837e-310,
   6.908783e-310,
   2.48515e-321,
   -1,
   1,
   1,
   4.656789e-310,
   4.656789e-310};
   graph = new TGraph(54,Graph2_fx11,Graph2_fy11);
   graph->SetName("Graph2");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#66cc66");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph211 = new TH1F("Graph_Graph211","",100,0,2.7225);
   Graph_Graph211->SetMinimum(-0.5);
   Graph_Graph211->SetMaximum(17);
   Graph_Graph211->SetDirectory(0);
   Graph_Graph211->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph211->SetLineColor(ci);
   Graph_Graph211->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph211->GetXaxis()->SetRange(25,89);
   Graph_Graph211->GetXaxis()->CenterTitle(true);
   Graph_Graph211->GetXaxis()->SetNdivisions(6);
   Graph_Graph211->GetXaxis()->SetLabelFont(132);
   Graph_Graph211->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph211->GetXaxis()->SetTitleSize(0);
   Graph_Graph211->GetXaxis()->SetTickLength(0.02);
   Graph_Graph211->GetXaxis()->SetTitleOffset(1);
   Graph_Graph211->GetXaxis()->SetTitleFont(132);
   Graph_Graph211->GetYaxis()->SetTitle("Acceptance Uncertainty [%]");
   Graph_Graph211->GetYaxis()->CenterTitle(true);
   Graph_Graph211->GetYaxis()->SetNdivisions(8);
   Graph_Graph211->GetYaxis()->SetLabelFont(132);
   Graph_Graph211->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph211->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph211->GetYaxis()->SetTickLength(0.02);
   Graph_Graph211->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph211->GetYaxis()->SetTitleFont(132);
   Graph_Graph211->GetZaxis()->SetLabelFont(42);
   Graph_Graph211->GetZaxis()->SetTitleOffset(1);
   Graph_Graph211->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph211);
   
   graph->Draw("p");
   line = new TLine(2.2611,-0.4,2.2611,17);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.8,0.9,"(e)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
   NewTwo_261->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: NewFour_461
   TPad *NewFour_461 = new TPad("NewFour_461", "NewFour_461",0.68,0.44,1,0.72);
   NewFour_461->Draw();
   NewFour_461->cd();
   NewFour_461->Range(1.447875,-0.5,4.743578,17);
   NewFour_461->SetFillColor(0);
   NewFour_461->SetBorderMode(0);
   NewFour_461->SetBorderSize(2);
   NewFour_461->SetLeftMargin(0);
   NewFour_461->SetRightMargin(0.04);
   NewFour_461->SetTopMargin(0);
   NewFour_461->SetBottomMargin(0);
   NewFour_461->SetFrameBorderMode(0);
   NewFour_461->SetFrameBorderMode(0);
   
   Double_t Graph0_fx12[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   6.952837e-310,
   6.908802e-310,
   6.908783e-310};
   Double_t Graph0_fy12[38] = {
   0,
   0,
   0,
   0,
   0,
   68.79974,
   13.02349,
   3.210566,
   1.852148,
   3.017496,
   2.314896,
   2.945206,
   3.181873,
   2.806228,
   1.647044,
   1.902682,
   1.496802,
   1.368032,
   1.877642,
   2.324841,
   1.057803,
   4.036195,
   4.036195,
   4.036195,
   4.036195,
   4.036195,
   4.036195,
   4.036195,
   4.036195,
   4.036195,
   4.036195,
   4.036195,
   0,
   0,
   0,
   4.65679e-310,
   4.940656e-324,
   4.65679e-310};
   graph = new TGraph(38,Graph0_fx12,Graph0_fy12);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph012 = new TH1F("Graph_Graph012","",100,0,5.3625);
   Graph_Graph012->SetMinimum(-0.5);
   Graph_Graph012->SetMaximum(17);
   Graph_Graph012->SetDirectory(0);
   Graph_Graph012->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph012->SetLineColor(ci);
   Graph_Graph012->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph012->GetXaxis()->SetRange(28,86);
   Graph_Graph012->GetXaxis()->CenterTitle(true);
   Graph_Graph012->GetXaxis()->SetNdivisions(6);
   Graph_Graph012->GetXaxis()->SetLabelFont(132);
   Graph_Graph012->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph012->GetXaxis()->SetTitleSize(0);
   Graph_Graph012->GetXaxis()->SetTickLength(0.02);
   Graph_Graph012->GetXaxis()->SetTitleOffset(1);
   Graph_Graph012->GetXaxis()->SetTitleFont(132);
   Graph_Graph012->GetYaxis()->SetTitle("Acceptance Uncertainty [%]");
   Graph_Graph012->GetYaxis()->CenterTitle(true);
   Graph_Graph012->GetYaxis()->SetNdivisions(8);
   Graph_Graph012->GetYaxis()->SetLabelFont(132);
   Graph_Graph012->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph012->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph012->GetYaxis()->SetTickLength(0.02);
   Graph_Graph012->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph012->GetYaxis()->SetTitleFont(132);
   Graph_Graph012->GetZaxis()->SetLabelFont(42);
   Graph_Graph012->GetZaxis()->SetTitleOffset(1);
   Graph_Graph012->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph012);
   
   graph->Draw("ap");
   
   Double_t Graph1_fx13[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   6.952837e-310,
   6.908802e-310,
   6.908783e-310};
   Double_t Graph1_fy13[38] = {
   0,
   0,
   0,
   0,
   0,
   57.73503,
   17.27365,
   11.86781,
   6.922317,
   5.448867,
   2.396836,
   2.074712,
   1.765266,
   1.319505,
   1.297572,
   1.011903,
   0.09240666,
   0.3157385,
   0.1278618,
   0.1362183,
   1.264836,
   3.737632,
   3.737632,
   3.737632,
   3.737632,
   3.737632,
   3.737632,
   3.737632,
   3.737632,
   3.737632,
   3.737632,
   3.737632,
   0,
   0,
   0,
   4.65679e-310,
   4.940656e-324,
   4.65679e-310};
   graph = new TGraph(38,Graph1_fx13,Graph1_fy13);
   graph->SetName("Graph1");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#cc66cc");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph113 = new TH1F("Graph_Graph113","",100,0,5.3625);
   Graph_Graph113->SetMinimum(-0.5);
   Graph_Graph113->SetMaximum(17);
   Graph_Graph113->SetDirectory(0);
   Graph_Graph113->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph113->SetLineColor(ci);
   Graph_Graph113->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph113->GetXaxis()->SetRange(28,86);
   Graph_Graph113->GetXaxis()->CenterTitle(true);
   Graph_Graph113->GetXaxis()->SetNdivisions(6);
   Graph_Graph113->GetXaxis()->SetLabelFont(132);
   Graph_Graph113->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph113->GetXaxis()->SetTitleSize(0);
   Graph_Graph113->GetXaxis()->SetTickLength(0.02);
   Graph_Graph113->GetXaxis()->SetTitleOffset(1);
   Graph_Graph113->GetXaxis()->SetTitleFont(132);
   Graph_Graph113->GetYaxis()->SetTitle("Acceptance Uncertainty [%]");
   Graph_Graph113->GetYaxis()->CenterTitle(true);
   Graph_Graph113->GetYaxis()->SetNdivisions(8);
   Graph_Graph113->GetYaxis()->SetLabelFont(132);
   Graph_Graph113->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph113->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph113->GetYaxis()->SetTickLength(0.02);
   Graph_Graph113->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph113->GetYaxis()->SetTitleFont(132);
   Graph_Graph113->GetZaxis()->SetLabelFont(42);
   Graph_Graph113->GetZaxis()->SetTitleOffset(1);
   Graph_Graph113->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph113);
   
   graph->Draw("p");
   
   Double_t Graph2_fx14[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   6.952837e-310,
   6.908802e-310,
   6.908783e-310};
   Double_t Graph2_fy14[38] = {
   0,
   0,
   0,
   0,
   0,
   57.73503,
   51.86582,
   11.36635,
   4.749323,
   1.210928,
   1.574308,
   1.23987,
   1.761288,
   0.9680436,
   1.635504,
   0.7016507,
   1.440976,
   2.036099,
   1.528925,
   2.047327,
   0.7554545,
   3.785882,
   3.785882,
   3.785882,
   3.785882,
   3.785882,
   3.785882,
   3.785882,
   3.785882,
   3.785882,
   3.785882,
   3.785882,
   0,
   0,
   0,
   4.65679e-310,
   4.940656e-324,
   4.65679e-310};
   graph = new TGraph(38,Graph2_fx14,Graph2_fy14);
   graph->SetName("Graph2");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#66cc66");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph214 = new TH1F("Graph_Graph214","",100,0,5.3625);
   Graph_Graph214->SetMinimum(-0.5);
   Graph_Graph214->SetMaximum(17);
   Graph_Graph214->SetDirectory(0);
   Graph_Graph214->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph214->SetLineColor(ci);
   Graph_Graph214->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph214->GetXaxis()->SetRange(28,86);
   Graph_Graph214->GetXaxis()->CenterTitle(true);
   Graph_Graph214->GetXaxis()->SetNdivisions(6);
   Graph_Graph214->GetXaxis()->SetLabelFont(132);
   Graph_Graph214->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph214->GetXaxis()->SetTitleSize(0);
   Graph_Graph214->GetXaxis()->SetTickLength(0.02);
   Graph_Graph214->GetXaxis()->SetTitleOffset(1);
   Graph_Graph214->GetXaxis()->SetTitleFont(132);
   Graph_Graph214->GetYaxis()->SetTitle("Acceptance Uncertainty [%]");
   Graph_Graph214->GetYaxis()->CenterTitle(true);
   Graph_Graph214->GetYaxis()->SetNdivisions(8);
   Graph_Graph214->GetYaxis()->SetLabelFont(132);
   Graph_Graph214->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph214->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph214->GetYaxis()->SetTickLength(0.02);
   Graph_Graph214->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph214->GetYaxis()->SetTitleFont(132);
   Graph_Graph214->GetZaxis()->SetLabelFont(42);
   Graph_Graph214->GetZaxis()->SetTitleOffset(1);
   Graph_Graph214->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph214);
   
   graph->Draw("p");
   line = new TLine(4.441,-0.4,4.441,17);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.8,0.9,"(f)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
   NewFour_461->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: NewPadOne_161
   TPad *NewPadOne_161 = new TPad("NewPadOne_161", "NewPadOne_161",0.02,0.72,0.36,1);
   NewPadOne_161->Draw();
   NewPadOne_161->cd();
   NewPadOne_161->Range(0.3613388,0.1,1.23748,14.03939);
   NewPadOne_161->SetFillColor(0);
   NewPadOne_161->SetBorderMode(0);
   NewPadOne_161->SetBorderSize(2);
   NewPadOne_161->SetLeftMargin(0.15);
   NewPadOne_161->SetRightMargin(0);
   NewPadOne_161->SetTopMargin(0.01);
   NewPadOne_161->SetBottomMargin(0);
   NewPadOne_161->SetFrameBorderMode(0);
   NewPadOne_161->SetFrameBorderMode(0);
   
   Double_t Graph0_fx15[38] = {
   0.42,
   0.46,
   0.5,
   0.54,
   0.58,
   0.62,
   0.66,
   0.7,
   0.74,
   0.78,
   0.82,
   0.86,
   0.9,
   0.94,
   0.98,
   1.02,
   1.06,
   1.09,
   1.11,
   1.13,
   1.15,
   1.17,
   1.19,
   1.21,
   1.23,
   1.25,
   1.27,
   1.29,
   1.31,
   1.33,
   1.35,
   1.37,
   1.39,
   1.41,
   1.43,
   1.45,
   1.47,
   1.49};
   Double_t Graph0_fy15[38] = {
   0,
   7.075893,
   11.69871,
   10.28445,
   9.500216,
   8.750488,
   8.107828,
   7.395401,
   7.180181,
   6.94997,
   6.842801,
   6.872849,
   6.737767,
   6.55723,
   6.689176,
   6.068722,
   6.05409,
   6.065336,
   6.063967,
   5.136935,
   5.240409,
   7.191909,
   5.222169,
   4.304375,
   7.765543,
   9.92206,
   6.032649,
   6.727366,
   2.658637,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(38,Graph0_fx15,Graph0_fy15);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph015 = new TH1F("Graph_Graph015","",100,0.313,1.597);
   Graph_Graph015->SetMinimum(0.1);
   Graph_Graph015->SetMaximum(13.9);
   Graph_Graph015->SetDirectory(0);
   Graph_Graph015->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph015->SetLineColor(ci);
   Graph_Graph015->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph015->GetXaxis()->SetRange(15,72);
   Graph_Graph015->GetXaxis()->CenterTitle(true);
   Graph_Graph015->GetXaxis()->SetNdivisions(6);
   Graph_Graph015->GetXaxis()->SetLabelFont(132);
   Graph_Graph015->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph015->GetXaxis()->SetTitleSize(0);
   Graph_Graph015->GetXaxis()->SetTickLength(0.02);
   Graph_Graph015->GetXaxis()->SetTitleOffset(1);
   Graph_Graph015->GetXaxis()->SetTitleFont(132);
   Graph_Graph015->GetYaxis()->SetTitle("Acceptance Correction");
   Graph_Graph015->GetYaxis()->CenterTitle(true);
   Graph_Graph015->GetYaxis()->SetNdivisions(8);
   Graph_Graph015->GetYaxis()->SetLabelFont(132);
   Graph_Graph015->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph015->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph015->GetYaxis()->SetTickLength(0.02);
   Graph_Graph015->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph015->GetYaxis()->SetTitleFont(132);
   Graph_Graph015->GetZaxis()->SetLabelFont(42);
   Graph_Graph015->GetZaxis()->SetTitleOffset(1);
   Graph_Graph015->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph015);
   
   graph->Draw("ap");
   line = new TLine(1.161,0.1,1.161,13.9);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.8,0.9,"(a)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.45,0.6,"^{12}C");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.23,0.9,"1.159 GeV");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
   NewPadOne_161->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: ExtraNewTwo_261
   TPad *ExtraNewTwo_261 = new TPad("ExtraNewTwo_261", "ExtraNewTwo_261",0.36,0.72,0.68,1);
   ExtraNewTwo_261->Draw();
   ExtraNewTwo_261->cd();
   ExtraNewTwo_261->Range(0.6558,0.1,2.42646,14.03939);
   ExtraNewTwo_261->SetFillColor(0);
   ExtraNewTwo_261->SetBorderMode(0);
   ExtraNewTwo_261->SetBorderSize(2);
   ExtraNewTwo_261->SetLeftMargin(0);
   ExtraNewTwo_261->SetRightMargin(0);
   ExtraNewTwo_261->SetTopMargin(0.01);
   ExtraNewTwo_261->SetBottomMargin(0);
   ExtraNewTwo_261->SetFrameBorderMode(0);
   ExtraNewTwo_261->SetFrameBorderMode(0);
   
   Double_t Graph0_fx16[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   2.505,
   2.535,
   2.565,
   2.595,
   2.625,
   2.655,
   2.685,
   2.715,
   2.745,
   2.775,
   2.805,
   2.835,
   2.865,
   2.895,
   2.925,
   2.955,
   2.985};
   Double_t Graph0_fy16[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   1.972934,
   4.015866,
   3.89163,
   3.936015,
   4.107296,
   4.25948,
   4.496864,
   4.732779,
   4.97421,
   5.291952,
   5.623029,
   5.938017,
   6.14838,
   6.415941,
   6.559432,
   6.90741,
   6.782865,
   7.303998,
   7.502884,
   7.486979,
   7.724768,
   7.22288,
   4.277636,
   6.793331,
   4.408682,
   6.554409,
   5.337473,
   3.217457,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(54,Graph0_fx16,Graph0_fy16);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph016 = new TH1F("Graph_Graph016","",100,0,3.279);
   Graph_Graph016->SetMinimum(0.1);
   Graph_Graph016->SetMaximum(13.9);
   Graph_Graph016->SetDirectory(0);
   Graph_Graph016->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph016->SetLineColor(ci);
   Graph_Graph016->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph016->GetXaxis()->SetRange(21,74);
   Graph_Graph016->GetXaxis()->CenterTitle(true);
   Graph_Graph016->GetXaxis()->SetNdivisions(6);
   Graph_Graph016->GetXaxis()->SetLabelFont(132);
   Graph_Graph016->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph016->GetXaxis()->SetTitleSize(0);
   Graph_Graph016->GetXaxis()->SetTickLength(0.02);
   Graph_Graph016->GetXaxis()->SetTitleOffset(1);
   Graph_Graph016->GetXaxis()->SetTitleFont(132);
   Graph_Graph016->GetYaxis()->SetTitle("Acceptance Correction");
   Graph_Graph016->GetYaxis()->CenterTitle(true);
   Graph_Graph016->GetYaxis()->SetNdivisions(8);
   Graph_Graph016->GetYaxis()->SetLabelFont(132);
   Graph_Graph016->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph016->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph016->GetYaxis()->SetTickLength(0.02);
   Graph_Graph016->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph016->GetYaxis()->SetTitleFont(132);
   Graph_Graph016->GetZaxis()->SetLabelFont(42);
   Graph_Graph016->GetZaxis()->SetTitleOffset(1);
   Graph_Graph016->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph016);
   
   graph->Draw("ap");
   
   Double_t Graph1_fx17[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   2.505,
   2.535,
   2.565,
   2.595,
   2.625,
   2.655,
   2.685,
   2.715,
   2.745,
   2.775,
   2.805,
   2.835,
   2.865,
   2.895,
   2.925,
   2.955,
   2.985};
   Double_t Graph1_fy17[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   2.136609,
   3.652936,
   3.615403,
   3.624344,
   3.749257,
   3.868758,
   4.026677,
   4.281995,
   4.519169,
   4.874677,
   5.191637,
   5.574812,
   5.969105,
   6.341089,
   6.804264,
   7.273703,
   7.328168,
   7.750076,
   8.021397,
   7.997495,
   8.801497,
   8.643009,
   5.847297,
   7.915361,
   4.820453,
   7.553698,
   5.567097,
   10.65857,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(54,Graph1_fx17,Graph1_fy17);
   graph->SetName("Graph1");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#cc66cc");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph117 = new TH1F("Graph_Graph117","",100,0,3.279);
   Graph_Graph117->SetMinimum(0.1);
   Graph_Graph117->SetMaximum(13.9);
   Graph_Graph117->SetDirectory(0);
   Graph_Graph117->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph117->SetLineColor(ci);
   Graph_Graph117->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph117->GetXaxis()->SetRange(21,74);
   Graph_Graph117->GetXaxis()->CenterTitle(true);
   Graph_Graph117->GetXaxis()->SetNdivisions(6);
   Graph_Graph117->GetXaxis()->SetLabelFont(132);
   Graph_Graph117->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph117->GetXaxis()->SetTitleSize(0);
   Graph_Graph117->GetXaxis()->SetTickLength(0.02);
   Graph_Graph117->GetXaxis()->SetTitleOffset(1);
   Graph_Graph117->GetXaxis()->SetTitleFont(132);
   Graph_Graph117->GetYaxis()->SetTitle("Acceptance Correction");
   Graph_Graph117->GetYaxis()->CenterTitle(true);
   Graph_Graph117->GetYaxis()->SetNdivisions(8);
   Graph_Graph117->GetYaxis()->SetLabelFont(132);
   Graph_Graph117->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph117->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph117->GetYaxis()->SetTickLength(0.02);
   Graph_Graph117->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph117->GetYaxis()->SetTitleFont(132);
   Graph_Graph117->GetZaxis()->SetLabelFont(42);
   Graph_Graph117->GetZaxis()->SetTitleOffset(1);
   Graph_Graph117->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph117);
   
   graph->Draw("p");
   
   Double_t Graph2_fx18[54] = {
   0.045,
   0.135,
   0.225,
   0.315,
   0.405,
   0.495,
   0.585,
   0.675,
   0.765,
   0.855,
   0.945,
   1.035,
   1.125,
   1.215,
   1.305,
   1.395,
   1.485,
   1.575,
   1.665,
   1.755,
   1.845,
   1.935,
   2.025,
   2.085,
   2.115,
   2.145,
   2.175,
   2.205,
   2.235,
   2.265,
   2.295,
   2.325,
   2.355,
   2.385,
   2.415,
   2.445,
   2.475,
   2.505,
   2.535,
   2.565,
   2.595,
   2.625,
   2.655,
   2.685,
   2.715,
   2.745,
   2.775,
   2.805,
   2.835,
   2.865,
   2.895,
   2.925,
   2.955,
   2.985};
   Double_t Graph2_fy18[54] = {
   0,
   0,
   0,
   0,
   0,
   0,
   1.952487,
   4.18263,
   3.922122,
   3.904914,
   4.002121,
   4.106637,
   4.386604,
   4.600437,
   4.856285,
   5.198433,
   5.584305,
   5.921814,
   6.142343,
   6.359284,
   6.654678,
   7.00322,
   6.782439,
   7.732957,
   7.924661,
   8.316862,
   8.873245,
   8.541245,
   5.673307,
   7.230272,
   6.504399,
   9.355305,
   5.657155,
   6.350804,
   1.44658,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(54,Graph2_fx18,Graph2_fy18);
   graph->SetName("Graph2");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#66cc66");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph218 = new TH1F("Graph_Graph218","",100,0,3.279);
   Graph_Graph218->SetMinimum(0.1);
   Graph_Graph218->SetMaximum(13.9);
   Graph_Graph218->SetDirectory(0);
   Graph_Graph218->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph218->SetLineColor(ci);
   Graph_Graph218->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph218->GetXaxis()->SetRange(21,74);
   Graph_Graph218->GetXaxis()->CenterTitle(true);
   Graph_Graph218->GetXaxis()->SetNdivisions(6);
   Graph_Graph218->GetXaxis()->SetLabelFont(132);
   Graph_Graph218->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph218->GetXaxis()->SetTitleSize(0);
   Graph_Graph218->GetXaxis()->SetTickLength(0.02);
   Graph_Graph218->GetXaxis()->SetTitleOffset(1);
   Graph_Graph218->GetXaxis()->SetTitleFont(132);
   Graph_Graph218->GetYaxis()->SetTitle("Acceptance Correction");
   Graph_Graph218->GetYaxis()->CenterTitle(true);
   Graph_Graph218->GetYaxis()->SetNdivisions(8);
   Graph_Graph218->GetYaxis()->SetLabelFont(132);
   Graph_Graph218->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph218->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph218->GetYaxis()->SetTickLength(0.02);
   Graph_Graph218->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph218->GetYaxis()->SetTitleFont(132);
   Graph_Graph218->GetZaxis()->SetLabelFont(42);
   Graph_Graph218->GetZaxis()->SetTitleOffset(1);
   Graph_Graph218->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph218);
   
   graph->Draw("p");
   line = new TLine(2.261,0.1,2.261,13.9);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.8,0.9,"(b)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.45,0.19,"^{56}Fe");
tex->SetNDC();

   ci = TColor::GetColor("#cc66cc");
   tex->SetTextColor(ci);
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.35,0.48,"^{4}He");
tex->SetNDC();

   ci = TColor::GetColor("#66cc66");
   tex->SetTextColor(ci);
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.1,0.9,"2.257 GeV");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
   ExtraNewTwo_261->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: ExtraNewFour_461
   TPad *ExtraNewFour_461 = new TPad("ExtraNewFour_461", "ExtraNewFour_461",0.68,0.72,1,1);
   ExtraNewFour_461->Draw();
   ExtraNewFour_461->cd();
   ExtraNewFour_461->Range(1.489725,0.1,4.765741,14.03939);
   ExtraNewFour_461->SetFillColor(0);
   ExtraNewFour_461->SetBorderMode(0);
   ExtraNewFour_461->SetBorderSize(2);
   ExtraNewFour_461->SetLeftMargin(0);
   ExtraNewFour_461->SetRightMargin(0.04);
   ExtraNewFour_461->SetTopMargin(0.01);
   ExtraNewFour_461->SetBottomMargin(0);
   ExtraNewFour_461->SetFrameBorderMode(0);
   ExtraNewFour_461->SetFrameBorderMode(0);
   
   Double_t Graph0_fx19[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   4.925,
   4.975,
   5.025};
   Double_t Graph0_fy19[38] = {
   0,
   0,
   0,
   0,
   0,
   19.78706,
   4.099523,
   1.371919,
   1.410425,
   1.463576,
   1.692234,
   1.951311,
   2.232224,
   2.601952,
   3.03945,
   3.451237,
   3.879521,
   4.279665,
   4.58084,
   4.748089,
   4.757689,
   4.443601,
   5.261392,
   5.099435,
   5.925565,
   3.50415,
   7.394255,
   2.478617,
   2.366617,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(38,Graph0_fx19,Graph0_fy19);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph019 = new TH1F("Graph_Graph019","",100,0,5.5175);
   Graph_Graph019->SetMinimum(0.1);
   Graph_Graph019->SetMaximum(13.9);
   Graph_Graph019->SetDirectory(0);
   Graph_Graph019->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph019->SetLineColor(ci);
   Graph_Graph019->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph019->GetXaxis()->SetRange(28,84);
   Graph_Graph019->GetXaxis()->CenterTitle(true);
   Graph_Graph019->GetXaxis()->SetNdivisions(6);
   Graph_Graph019->GetXaxis()->SetLabelFont(132);
   Graph_Graph019->GetXaxis()->SetLabelSize(0.09);
   Graph_Graph019->GetXaxis()->SetTitleSize(0);
   Graph_Graph019->GetXaxis()->SetTickLength(0.02);
   Graph_Graph019->GetXaxis()->SetTitleOffset(1);
   Graph_Graph019->GetXaxis()->SetTitleFont(132);
   Graph_Graph019->GetYaxis()->SetTitle("Acceptance Correction");
   Graph_Graph019->GetYaxis()->CenterTitle(true);
   Graph_Graph019->GetYaxis()->SetNdivisions(8);
   Graph_Graph019->GetYaxis()->SetLabelFont(132);
   Graph_Graph019->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph019->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph019->GetYaxis()->SetTickLength(0.02);
   Graph_Graph019->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph019->GetYaxis()->SetTitleFont(132);
   Graph_Graph019->GetZaxis()->SetLabelFont(42);
   Graph_Graph019->GetZaxis()->SetTitleOffset(1);
   Graph_Graph019->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph019);
   
   graph->Draw("ap");
   
   Double_t Graph1_fx20[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   4.925,
   4.975,
   5.025};
   Double_t Graph1_fy20[38] = {
   0,
   0,
   0,
   0,
   0,
   -2.148768,
   4.364832,
   1.714682,
   1.545323,
   1.663759,
   1.828234,
   2.130858,
   2.450889,
   2.829176,
   3.21088,
   3.640609,
   4.082938,
   4.501426,
   4.776271,
   4.937823,
   5.056278,
   4.81974,
   5.309961,
   5.452607,
   5.915678,
   3.742887,
   7.529018,
   2.468161,
   2.473505,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(38,Graph1_fx20,Graph1_fy20);
   graph->SetName("Graph1");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#cc66cc");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph120 = new TH1F("Graph_Graph120","",100,0,5.5175);
   Graph_Graph120->SetMinimum(0.1);
   Graph_Graph120->SetMaximum(13.9);
   Graph_Graph120->SetDirectory(0);
   Graph_Graph120->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph120->SetLineColor(ci);
   Graph_Graph120->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph120->GetXaxis()->SetRange(28,84);
   Graph_Graph120->GetXaxis()->CenterTitle(true);
   Graph_Graph120->GetXaxis()->SetNdivisions(6);
   Graph_Graph120->GetXaxis()->SetLabelFont(132);
   Graph_Graph120->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph120->GetXaxis()->SetTitleSize(0);
   Graph_Graph120->GetXaxis()->SetTickLength(0.02);
   Graph_Graph120->GetXaxis()->SetTitleOffset(1);
   Graph_Graph120->GetXaxis()->SetTitleFont(132);
   Graph_Graph120->GetYaxis()->SetTitle("Acceptance Correction");
   Graph_Graph120->GetYaxis()->CenterTitle(true);
   Graph_Graph120->GetYaxis()->SetNdivisions(8);
   Graph_Graph120->GetYaxis()->SetLabelFont(132);
   Graph_Graph120->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph120->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph120->GetYaxis()->SetTickLength(0.02);
   Graph_Graph120->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph120->GetYaxis()->SetTitleFont(132);
   Graph_Graph120->GetZaxis()->SetLabelFont(42);
   Graph_Graph120->GetZaxis()->SetTitleOffset(1);
   Graph_Graph120->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph120);
   
   graph->Draw("p");
   
   Double_t Graph2_fx21[38] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.1,
   2.3,
   2.5,
   2.7,
   2.9,
   3.1,
   3.3,
   3.5,
   3.7,
   3.9,
   4.1,
   4.225,
   4.275,
   4.325,
   4.375,
   4.425,
   4.475,
   4.525,
   4.575,
   4.625,
   4.675,
   4.725,
   4.775,
   4.825,
   4.875,
   4.925,
   4.975,
   5.025};
   Double_t Graph2_fy21[38] = {
   0,
   0,
   0,
   0,
   0,
   73.89358,
   4.894189,
   1.962705,
   1.689432,
   1.66891,
   1.752777,
   2.041255,
   2.364319,
   2.747459,
   3.156859,
   3.531258,
   3.924415,
   4.34064,
   4.64755,
   4.721483,
   4.622182,
   4.222096,
   5.506883,
   5.021066,
   6.329356,
   3.436007,
   7.164978,
   3.110794,
   3.05767,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   graph = new TGraph(38,Graph2_fx21,Graph2_fy21);
   graph->SetName("Graph2");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#66cc66");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph221 = new TH1F("Graph_Graph221","",100,0,5.5175);
   Graph_Graph221->SetMinimum(0.1);
   Graph_Graph221->SetMaximum(13.9);
   Graph_Graph221->SetDirectory(0);
   Graph_Graph221->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph221->SetLineColor(ci);
   Graph_Graph221->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} E^{cal} [GeV]");
   Graph_Graph221->GetXaxis()->SetRange(28,84);
   Graph_Graph221->GetXaxis()->CenterTitle(true);
   Graph_Graph221->GetXaxis()->SetNdivisions(6);
   Graph_Graph221->GetXaxis()->SetLabelFont(132);
   Graph_Graph221->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph221->GetXaxis()->SetTitleSize(0);
   Graph_Graph221->GetXaxis()->SetTickLength(0.02);
   Graph_Graph221->GetXaxis()->SetTitleOffset(1);
   Graph_Graph221->GetXaxis()->SetTitleFont(132);
   Graph_Graph221->GetYaxis()->SetTitle("Acceptance Correction");
   Graph_Graph221->GetYaxis()->CenterTitle(true);
   Graph_Graph221->GetYaxis()->SetNdivisions(8);
   Graph_Graph221->GetYaxis()->SetLabelFont(132);
   Graph_Graph221->GetYaxis()->SetLabelSize(0.09);
   Graph_Graph221->GetYaxis()->SetTitleSize(0.09);
   Graph_Graph221->GetYaxis()->SetTickLength(0.02);
   Graph_Graph221->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph221->GetYaxis()->SetTitleFont(132);
   Graph_Graph221->GetZaxis()->SetLabelFont(42);
   Graph_Graph221->GetZaxis()->SetTitleOffset(1);
   Graph_Graph221->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph221);
   
   graph->Draw("p");
   line = new TLine(4.461,0.1,4.461,13.9);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.8,0.9,"(c)");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.1,0.9,"4.453 GeV");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();
   ExtraNewFour_461->Modified();
   CommonPanel_NoxBCut->cd();
  
// ------------>Primitives in pad: Ecal
   TPad *Ecal = new TPad("Ecal", "Ecal",0,0.03,1,0.13);
   Ecal->Draw();
   Ecal->cd();
   Ecal->Range(0,0,1,1);
   Ecal->SetFillColor(0);
   Ecal->SetBorderMode(0);
   Ecal->SetBorderSize(2);
   Ecal->SetFrameBorderMode(0);
      tex = new TLatex(0.42,0.5,"(e,e'p)_{1p0#pi} E_{cal} [GeV]");
tex->SetNDC();
   tex->SetTextFont(132);
   tex->SetTextSize(0.45);
   tex->SetLineWidth(2);
   tex->Draw();
   Ecal->Modified();
   CommonPanel_NoxBCut->cd();
   CommonPanel_NoxBCut->Modified();
   CommonPanel_NoxBCut->cd();
   CommonPanel_NoxBCut->SetSelected(CommonPanel_NoxBCut);
}
