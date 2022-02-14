void DataXSec_Inclusive_Validation()
{
//=========Macro generated from canvas: DataXSec_Inclusive_Validation/DataXSec_Inclusive_Validation
//=========  (Sun Feb 13 16:57:22 2022) by ROOT version 6.24/02
   TCanvas *DataXSec_Inclusive_Validation = new TCanvas("DataXSec_Inclusive_Validation", "DataXSec_Inclusive_Validation",205,71,1024,768);
   gStyle->SetOptStat(0);
   DataXSec_Inclusive_Validation->Range(-0.09876002,-1862.781,0.95724,9094.754);
   DataXSec_Inclusive_Validation->SetFillColor(0);
   DataXSec_Inclusive_Validation->SetBorderMode(0);
   DataXSec_Inclusive_Validation->SetBorderSize(2);
   DataXSec_Inclusive_Validation->SetLeftMargin(0.18);
   DataXSec_Inclusive_Validation->SetBottomMargin(0.17);
   DataXSec_Inclusive_Validation->SetFrameBorderMode(0);
   DataXSec_Inclusive_Validation->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1001[65] = {
   0.163,
   0.173,
   0.183,
   0.193,
   0.203,
   0.213,
   0.223,
   0.233,
   0.243,
   0.253,
   0.263,
   0.273,
   0.283,
   0.293,
   0.303,
   0.313,
   0.323,
   0.333,
   0.343,
   0.353,
   0.363,
   0.373,
   0.383,
   0.393,
   0.403,
   0.413,
   0.423,
   0.433,
   0.443,
   0.453,
   0.463,
   0.473,
   0.483,
   0.493,
   0.503,
   0.513,
   0.523,
   0.533,
   0.543,
   0.553,
   0.563,
   0.573,
   0.583,
   0.593,
   0.603,
   0.613,
   0.623,
   0.633,
   0.643,
   0.653,
   0.663,
   0.673,
   0.683,
   0.693,
   0.703,
   0.713,
   0.723,
   0.733,
   0.743,
   0.753,
   0.763,
   0.773,
   0.783,
   0.793,
   0.803};
   Double_t Graph0_fy1001[65] = {
   258,
   294,
   398.4,
   465.6,
   612,
   735.6,
   848.4,
   951.6,
   1168,
   1342,
   1450,
   1457,
   1589,
   1763,
   1751,
   1750,
   1843,
   1756,
   1829,
   1712,
   1704,
   1621,
   1571,
   1486,
   1426,
   1355,
   1354,
   1343,
   1440,
   1355,
   1339,
   1418,
   1398,
   1448,
   1475,
   1607,
   1571,
   1634,
   1692,
   1693,
   1598,
   1654,
   1654,
   1678,
   1640,
   1625,
   1679,
   1613,
   1598,
   1543,
   1480,
   1529,
   1496,
   1512,
   1386,
   1337,
   1435,
   1408,
   1296,
   1414,
   1432,
   1495,
   1402,
   1362,
   1517};
   Double_t Graph0_fex1001[65] = {
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
   0,
   0};
   Double_t Graph0_fey1001[65] = {
   18.84,
   19.56,
   17.76,
   18.96,
   27.72,
   38.52,
   41.28,
   43.8,
   47.28,
   39.24,
   38.76,
   47.16,
   58.56,
   61.08,
   60.72,
   61.08,
   44.52,
   40.56,
   55.32,
   53.52,
   52.92,
   51.12,
   66.72,
   50.52,
   49.2,
   47.88,
   46.8,
   47.28,
   37.44,
   32.4,
   40.8,
   45.6,
   45.6,
   46.44,
   35.28,
   41.88,
   55.8,
   56.52,
   57.12,
   44.88,
   40.92,
   55.2,
   55.68,
   55.92,
   41.4,
   44.64,
   57.84,
   55.8,
   50.28,
   38.52,
   52.44,
   52.44,
   52.32,
   36.24,
   44.04,
   46.56,
   48.24,
   33.6,
   45.24,
   47.16,
   40.08,
   36.48,
   47.04,
   46.32,
   48.84};
   TGraphErrors *gre = new TGraphErrors(65,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
   gre->SetName("Graph0");
   gre->SetTitle("");
   gre->SetFillStyle(1000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000cc");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#0000cc");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph01001 = new TH1F("Graph_Graph01001","",100,0.099,0.867);
   Graph_Graph01001->SetMinimum(0);
   Graph_Graph01001->SetMaximum(7999);
   Graph_Graph01001->SetDirectory(0);
   Graph_Graph01001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph01001->SetLineColor(ci);
   Graph_Graph01001->GetXaxis()->SetTitle("Energy Transfer [GeV]");
   Graph_Graph01001->GetXaxis()->SetRange(0,98);
   Graph_Graph01001->GetXaxis()->CenterTitle(true);
   Graph_Graph01001->GetXaxis()->SetNdivisions(8);
   Graph_Graph01001->GetXaxis()->SetLabelFont(132);
   Graph_Graph01001->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph01001->GetXaxis()->SetLabelSize(0.07);
   Graph_Graph01001->GetXaxis()->SetTitleSize(0.07);
   Graph_Graph01001->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph01001->GetXaxis()->SetTitleFont(132);
   Graph_Graph01001->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#OmegadE_{e}} #left[#frac{#mub}{sr GeV ^{12}C}#right]");
   Graph_Graph01001->GetYaxis()->CenterTitle(true);
   Graph_Graph01001->GetYaxis()->SetNdivisions(7);
   Graph_Graph01001->GetYaxis()->SetLabelFont(132);
   Graph_Graph01001->GetYaxis()->SetLabelSize(0.07);
   Graph_Graph01001->GetYaxis()->SetTitleSize(0.07);
   Graph_Graph01001->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph01001->GetYaxis()->SetTitleFont(132);
   Graph_Graph01001->GetZaxis()->SetLabelFont(42);
   Graph_Graph01001->GetZaxis()->SetTitleOffset(1);
   Graph_Graph01001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph01001);
   
   gre->Draw("ap");
   
   Double_t Graph1_fx1002[60] = {
   0.073,
   0.083,
   0.093,
   0.103,
   0.113,
   0.123,
   0.133,
   0.143,
   0.153,
   0.163,
   0.173,
   0.183,
   0.193,
   0.203,
   0.213,
   0.223,
   0.233,
   0.243,
   0.253,
   0.263,
   0.273,
   0.283,
   0.293,
   0.303,
   0.313,
   0.323,
   0.333,
   0.343,
   0.353,
   0.363,
   0.373,
   0.383,
   0.393,
   0.403,
   0.413,
   0.423,
   0.433,
   0.443,
   0.453,
   0.463,
   0.473,
   0.483,
   0.493,
   0.503,
   0.513,
   0.523,
   0.533,
   0.543,
   0.553,
   0.563,
   0.573,
   0.583,
   0.593,
   0.603,
   0.613,
   0.623,
   0.633,
   0.643,
   0.653,
   0.663};
   Double_t Graph1_fy1002[60] = {
   746.4,
   1145,
   1596,
   2221,
   2861,
   3482,
   4166,
   5167,
   5298,
   5730,
   6383,
   6638,
   6790,
   6601,
   6155,
   6013,
   5813,
   5322,
   4813,
   4403,
   3913,
   3964,
   3614,
   3622,
   3341,
   3436,
   3547,
   3400,
   3731,
   3704,
   3906,
   4057,
   4132,
   4081,
   4344,
   4112,
   4666,
   4508,
   4403,
   4510,
   4362,
   4434,
   4246,
   4082,
   3935,
   3660,
   3527,
   3524,
   3344,
   3390,
   3170,
   3115,
   3018,
   2969,
   2845,
   2942,
   2676,
   2795,
   2868,
   2998};
   Double_t Graph1_fex1002[60] = {
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
   Double_t Graph1_fey1002[60] = {
   83.04,
   79.32,
   93.12,
   108.2,
   124.2,
   135.8,
   129.6,
   132.7,
   194.3,
   201.4,
   209.9,
   181.2,
   152,
   166.2,
   201.5,
   199,
   172,
   122.3,
   143.9,
   146.9,
   136.7,
   113.9,
   93.96,
   122,
   117.5,
   118.8,
   86.76,
   112.3,
   129.4,
   127.8,
   100.3,
   128.2,
   143.6,
   142.3,
   104.2,
   142.6,
   152.8,
   127,
   119.6,
   154.6,
   150.6,
   104.4,
   139.7,
   135,
   95.88,
   118.8,
   116.3,
   89.64,
   125.3,
   123.1,
   94.92,
   119.5,
   101.2,
   108.8,
   117.7,
   91.68,
   109.9,
   78,
   93.96,
   126.4};
   gre = new TGraphErrors(60,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#ff6600");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#ff6600");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph11002 = new TH1F("Graph_Graph11002","Graph",100,0.014,0.722);
   Graph_Graph11002->SetMinimum(35.496);
   Graph_Graph11002->SetMaximum(7569.864);
   Graph_Graph11002->SetDirectory(0);
   Graph_Graph11002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph11002->SetLineColor(ci);
   Graph_Graph11002->GetXaxis()->SetLabelFont(42);
   Graph_Graph11002->GetXaxis()->SetTitleOffset(1);
   Graph_Graph11002->GetXaxis()->SetTitleFont(42);
   Graph_Graph11002->GetYaxis()->SetLabelFont(42);
   Graph_Graph11002->GetYaxis()->SetTitleFont(42);
   Graph_Graph11002->GetZaxis()->SetLabelFont(42);
   Graph_Graph11002->GetZaxis()->SetTitleOffset(1);
   Graph_Graph11002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph11002);
   
   gre->Draw("p");
   
   Double_t Graph2_fx1[30] = {
   0.06166667,
   0.085,
   0.1083333,
   0.1316667,
   0.155,
   0.1783333,
   0.2016667,
   0.225,
   0.2483333,
   0.2716667,
   0.295,
   0.3183333,
   0.3416667,
   0.365,
   0.3883333,
   0.4116667,
   0.435,
   0.4583333,
   0.4816667,
   0.505,
   0.5283333,
   0.5516667,
   0.575,
   0.5983333,
   0.6216667,
   0.645,
   0.6683333,
   0.6916667,
   0.715,
   0.7383333};
   Double_t Graph2_fy1[30] = {
   19.20415,
   69.71942,
   192.5982,
   486.0877,
   1003.904,
   1688.295,
   2403.859,
   2934.617,
   3236.039,
   3279.874,
   3086.163,
   2818.557,
   2467.733,
   2181.898,
   2096.731,
   2255.653,
   2572.243,
   3013.66,
   3440.326,
   3754.272,
   3940.33,
   3841.248,
   3678.847,
   3430.307,
   3197.491,
   2961.892,
   2794.204,
   2685.38,
   2637.788,
   2656.157};
   TGraph *graph = new TGraph(30,Graph2_fx1,Graph2_fy1);
   graph->SetName("Graph2");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#009900");
   graph->SetLineColor(ci);
   graph->SetLineWidth(2);
   
   TH1F *Graph_Graph21 = new TH1F("Graph_Graph21","",100,0,0.806);
   Graph_Graph21->SetMinimum(17.28374);
   Graph_Graph21->SetMaximum(4332.443);
   Graph_Graph21->SetDirectory(0);
   Graph_Graph21->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph21->SetLineColor(ci);
   Graph_Graph21->GetXaxis()->SetTitle("Energy Transfer [GeV]");
   Graph_Graph21->GetXaxis()->SetRange(7,94);
   Graph_Graph21->GetXaxis()->CenterTitle(true);
   Graph_Graph21->GetXaxis()->SetNdivisions(8);
   Graph_Graph21->GetXaxis()->SetLabelFont(132);
   Graph_Graph21->GetXaxis()->SetLabelSize(0.07);
   Graph_Graph21->GetXaxis()->SetTitleSize(0.07);
   Graph_Graph21->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph21->GetXaxis()->SetTitleFont(132);
   Graph_Graph21->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#Omega dE} [#mub/sr/GeV]");
   Graph_Graph21->GetYaxis()->CenterTitle(true);
   Graph_Graph21->GetYaxis()->SetNdivisions(8);
   Graph_Graph21->GetYaxis()->SetLabelFont(132);
   Graph_Graph21->GetYaxis()->SetLabelSize(0.07);
   Graph_Graph21->GetYaxis()->SetTitleSize(0.07);
   Graph_Graph21->GetYaxis()->SetTitleOffset(1);
   Graph_Graph21->GetYaxis()->SetTitleFont(132);
   Graph_Graph21->GetZaxis()->SetLabelFont(42);
   Graph_Graph21->GetZaxis()->SetTitleOffset(1);
   Graph_Graph21->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph21);
   
   graph->Draw("c ");
   
   Double_t Graph3_fx2[10] = {
   0.137,
   0.211,
   0.285,
   0.359,
   0.433,
   0.507,
   0.581,
   0.655,
   0.729,
   0.803};
   Double_t Graph3_fy2[10] = {
   101.3181,
   824.7487,
   1805.736,
   1998.903,
   1586.055,
   1902.319,
   2272.557,
   2122,
   1937.355,
   1917.47};
   graph = new TGraph(10,Graph3_fx2,Graph3_fy2);
   graph->SetName("Graph3");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#0000cc");
   graph->SetLineColor(ci);
   graph->SetLineWidth(2);
   
   TH1F *Graph_Graph32 = new TH1F("Graph_Graph32","",100,0.0704,0.8696);
   Graph_Graph32->SetMinimum(91.18633);
   Graph_Graph32->SetMaximum(2489.68);
   Graph_Graph32->SetDirectory(0);
   Graph_Graph32->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph32->SetLineColor(ci);
   Graph_Graph32->GetXaxis()->SetLabelFont(42);
   Graph_Graph32->GetXaxis()->SetTitleOffset(1);
   Graph_Graph32->GetXaxis()->SetTitleFont(42);
   Graph_Graph32->GetYaxis()->SetLabelFont(42);
   Graph_Graph32->GetYaxis()->SetTitleFont(42);
   Graph_Graph32->GetZaxis()->SetLabelFont(42);
   Graph_Graph32->GetZaxis()->SetTitleOffset(1);
   Graph_Graph32->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph32);
   
   graph->Draw("c ");
   
   Double_t Graph4_fx3[10] = {
   0.035,
   0.105,
   0.175,
   0.245,
   0.315,
   0.385,
   0.455,
   0.525,
   0.595,
   0.665};
   Double_t Graph4_fy3[10] = {
   198.1191,
   2984.638,
   7251.16,
   5815.064,
   3496.535,
   3936.68,
   5279.607,
   4765.568,
   3803.887,
   3220.239};
   graph = new TGraph(10,Graph4_fx3,Graph4_fy3);
   graph->SetName("Graph4");
   graph->SetTitle("");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#ff6600");
   graph->SetLineColor(ci);
   graph->SetLineWidth(2);
   
   TH1F *Graph_Graph43 = new TH1F("Graph_Graph43","",100,0,0.728);
   Graph_Graph43->SetMinimum(178.3072);
   Graph_Graph43->SetMaximum(7956.464);
   Graph_Graph43->SetDirectory(0);
   Graph_Graph43->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph43->SetLineColor(ci);
   Graph_Graph43->GetXaxis()->SetLabelFont(42);
   Graph_Graph43->GetXaxis()->SetTitleOffset(1);
   Graph_Graph43->GetXaxis()->SetTitleFont(42);
   Graph_Graph43->GetYaxis()->SetLabelFont(42);
   Graph_Graph43->GetYaxis()->SetTitleFont(42);
   Graph_Graph43->GetZaxis()->SetLabelFont(42);
   Graph_Graph43->GetZaxis()->SetTitleOffset(1);
   Graph_Graph43->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph43);
   
   graph->Draw("c ");
   
   TH1F *BinCentCorrDatae4vPlot__1 = new TH1F("BinCentCorrDatae4vPlot__1","(e,e') ^{12}C, #theta = 37.5^{o}",187,0,5.984);
   BinCentCorrDatae4vPlot__1->SetBinContent(3,70.78925);
   BinCentCorrDatae4vPlot__1->SetBinContent(4,258.0344);
   BinCentCorrDatae4vPlot__1->SetBinContent(5,987.4088);
   BinCentCorrDatae4vPlot__1->SetBinContent(6,2055.191);
   BinCentCorrDatae4vPlot__1->SetBinContent(7,3070.521);
   BinCentCorrDatae4vPlot__1->SetBinContent(8,3650.882);
   BinCentCorrDatae4vPlot__1->SetBinContent(9,3714.26);
   BinCentCorrDatae4vPlot__1->SetBinContent(10,3121.552);
   BinCentCorrDatae4vPlot__1->SetBinContent(11,2556.01);
   BinCentCorrDatae4vPlot__1->SetBinContent(12,2268.246);
   BinCentCorrDatae4vPlot__1->SetBinContent(13,2545.049);
   BinCentCorrDatae4vPlot__1->SetBinContent(14,2635.459);
   BinCentCorrDatae4vPlot__1->SetBinContent(15,3401.744);
   BinCentCorrDatae4vPlot__1->SetBinContent(16,3452.932);
   BinCentCorrDatae4vPlot__1->SetBinContent(17,3380.339);
   BinCentCorrDatae4vPlot__1->SetBinContent(18,3471.404);
   BinCentCorrDatae4vPlot__1->SetBinContent(19,3292.409);
   BinCentCorrDatae4vPlot__1->SetBinContent(20,3112.284);
   BinCentCorrDatae4vPlot__1->SetBinContent(21,2709.086);
   BinCentCorrDatae4vPlot__1->SetBinContent(22,2763.446);
   BinCentCorrDatae4vPlot__1->SetBinContent(23,2893.911);
   BinCentCorrDatae4vPlot__1->SetBinContent(24,1871.229);
   BinCentCorrDatae4vPlot__1->SetBinError(3,41.51107);
   BinCentCorrDatae4vPlot__1->SetBinError(4,82.39747);
   BinCentCorrDatae4vPlot__1->SetBinError(5,84.30856);
   BinCentCorrDatae4vPlot__1->SetBinError(6,138.921);
   BinCentCorrDatae4vPlot__1->SetBinError(7,202.9588);
   BinCentCorrDatae4vPlot__1->SetBinError(8,239.2455);
   BinCentCorrDatae4vPlot__1->SetBinError(9,243.0985);
   BinCentCorrDatae4vPlot__1->SetBinError(10,205.2812);
   BinCentCorrDatae4vPlot__1->SetBinError(11,168.6894);
   BinCentCorrDatae4vPlot__1->SetBinError(12,160.1445);
   BinCentCorrDatae4vPlot__1->SetBinError(13,175.215);
   BinCentCorrDatae4vPlot__1->SetBinError(14,173.8716);
   BinCentCorrDatae4vPlot__1->SetBinError(15,225.5724);
   BinCentCorrDatae4vPlot__1->SetBinError(16,227.1454);
   BinCentCorrDatae4vPlot__1->SetBinError(17,222.9371);
   BinCentCorrDatae4vPlot__1->SetBinError(18,233.2894);
   BinCentCorrDatae4vPlot__1->SetBinError(19,218.889);
   BinCentCorrDatae4vPlot__1->SetBinError(20,205.9399);
   BinCentCorrDatae4vPlot__1->SetBinError(21,179.5176);
   BinCentCorrDatae4vPlot__1->SetBinError(22,182.5302);
   BinCentCorrDatae4vPlot__1->SetBinError(23,191.4944);
   BinCentCorrDatae4vPlot__1->SetBinError(24,127.854);
   BinCentCorrDatae4vPlot__1->SetEntries(57904);

   ci = TColor::GetColor("#009900");
   BinCentCorrDatae4vPlot__1->SetLineColor(ci);

   ci = TColor::GetColor("#009900");
   BinCentCorrDatae4vPlot__1->SetMarkerColor(ci);
   BinCentCorrDatae4vPlot__1->SetMarkerStyle(20);
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetTitle("Energy Transfer [GeV]");
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetRange(3,23);
   BinCentCorrDatae4vPlot__1->GetXaxis()->CenterTitle(true);
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetNdivisions(9);
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetLabelFont(132);
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetLabelOffset(0.015);
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetLabelSize(0.07);
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetTitleSize(0.07);
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetTitleOffset(1.1);
   BinCentCorrDatae4vPlot__1->GetXaxis()->SetTitleFont(132);
   BinCentCorrDatae4vPlot__1->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#OmegadE_{e}} [#frac{nb}{sr GeV}]");
   BinCentCorrDatae4vPlot__1->GetYaxis()->CenterTitle(true);
   BinCentCorrDatae4vPlot__1->GetYaxis()->SetNdivisions(7);
   BinCentCorrDatae4vPlot__1->GetYaxis()->SetLabelFont(132);
   BinCentCorrDatae4vPlot__1->GetYaxis()->SetLabelSize(0.07);
   BinCentCorrDatae4vPlot__1->GetYaxis()->SetTitleSize(0.07);
   BinCentCorrDatae4vPlot__1->GetYaxis()->SetTitleOffset(1.1);
   BinCentCorrDatae4vPlot__1->GetYaxis()->SetTitleFont(132);
   BinCentCorrDatae4vPlot__1->GetZaxis()->SetLabelFont(42);
   BinCentCorrDatae4vPlot__1->GetZaxis()->SetTitleOffset(1);
   BinCentCorrDatae4vPlot__1->GetZaxis()->SetTitleFont(42);
   BinCentCorrDatae4vPlot__1->Draw("e1x0 same");
   TLatex *   tex = new TLatex(0.5,0.22,"SLAC 1.299 GeV");
tex->SetNDC();

   ci = TColor::GetColor("#0000cc");
   tex->SetTextColor(ci);
   tex->SetTextFont(132);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.32,0.8,"SLAC 0.961 GeV");
tex->SetNDC();

   ci = TColor::GetColor("#ff6600");
   tex->SetTextColor(ci);
   tex->SetTextFont(132);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7,0.57,"#splitline{JLab}{1.159 GeV}");
tex->SetNDC();
   tex->SetTextAlign(12);

   ci = TColor::GetColor("#009900");
   tex->SetTextColor(ci);
   tex->SetTextFont(132);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.86,"#theta = 37.5^{o}");
tex->SetNDC();
   tex->SetTextAlign(12);
   tex->SetTextFont(132);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
   DataXSec_Inclusive_Validation->Modified();
   DataXSec_Inclusive_Validation->cd();
   DataXSec_Inclusive_Validation->SetSelected(DataXSec_Inclusive_Validation);
}
