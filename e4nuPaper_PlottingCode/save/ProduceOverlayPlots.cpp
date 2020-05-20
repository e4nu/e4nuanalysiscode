{
gROOT->ProcessLine(".L /home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp");
gROOT->ProcessLine(".L /home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp");
gROOT->ProcessLine(".L /home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp");

// Main e4nu Paper Plots

//gROOT->ProcessLine(".L OverlayEQE_Fig2.cpp"); gROOT->ProcessLine("OverlayEQE_Fig2()");

//gROOT->ProcessLine(".L OverlayPmissFig3a_e4nuPaper.cpp"); gROOT->ProcessLine("OverlayPmissFig3a_e4nuPaper()");

//gROOT->ProcessLine(".L OverlayPmissFig3b_e4nuPaper.cpp"); gROOT->ProcessLine("OverlayPmissFig3b_e4nuPaper()");

//gROOT->ProcessLine(".L OverlayECalFig4_e4nuPaper.cpp"); gROOT->ProcessLine("OverlayECalFig4_e4nuPaper()");

// ----------------------------------------------------------------------------------------

// Extended Data Figures

gROOT->ProcessLine(".L ExtDataTable1.cpp"); gROOT->ProcessLine("ExtDataTable1()"); // do it again

//gROOT->ProcessLine(".L ECalExtDataFig6_e4nuPaper.cpp"); gROOT->ProcessLine("ECalExtDataFig6_e4nuPaper()"); // QE-like case

//gROOT->ProcessLine(".L OverlayQ2VsNu_FigExtData6.cpp"); gROOT->ProcessLine("OverlayQ2VsNu_FigExtData6()");

//gROOT->ProcessLine(".L OverlayMultiplicities_FigExtData7.cpp"); gROOT->ProcessLine("OverlayMultiplicities_FigExtData7()");

//gROOT->ProcessLine(".L OverlayReso_FigExtData8.cpp"); gROOT->ProcessLine("OverlayReso_FigExtData8()");

//gROOT->ProcessLine(".L OverlayECalVsEQE_FigExtData9.cpp"); gROOT->ProcessLine("OverlayECalVsEQE_FigExtData9()");

//gROOT->ProcessLine(".L OverlayPmiss_FigExtData10.cpp"); gROOT->ProcessLine("OverlayPmiss_FigExtData10()");

//gROOT->ProcessLine(".L OverlayECal_FigExtData11.cpp"); gROOT->ProcessLine("OverlayECal_FigExtData11()");


// ----------------------------------------------------------------------------------------

// Playground Code

//gROOT->ProcessLine(".L OverlayPlots.cpp"); gROOT->ProcessLine("OverlayPlots()");

//gROOT->ProcessLine(".L ../GenieSyst.cpp"); gROOT->ProcessLine("GenieSyst()");

//gROOT->ProcessLine(".L BreakDown.cpp"); gROOT->ProcessLine("BreakDown()");

//gROOT->ProcessLine(".L OverlayReso.cpp"); gROOT->ProcessLine("OverlayReso()");

//gROOT->ProcessLine(".L ../Create2DPlots.cpp++"); gROOT->ProcessLine("Create2DPlots()");

//gROOT->ProcessLine(".L OverlaySlices.cpp++"); gROOT->ProcessLine("OverlaySlices()"); // Not now

//gROOT->ProcessLine(".q");
};
