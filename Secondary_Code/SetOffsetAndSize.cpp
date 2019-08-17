#include <TStyle.h>
#include <TROOT.h>

void SetOffsetAndSize(){

	gStyle->SetTitleSize(0.06,"t");

	gStyle->SetTitleOffset(0.8,"x");
	gStyle->SetTitleSize(0.05,"x");

	gStyle->SetTitleOffset(0.7,"y");
	gStyle->SetTitleSize(0.06,"y");

	gStyle->SetStatX(0.9);                
	gStyle->SetStatY(0.9);  
	gStyle->SetStatH(0.2);

	gStyle->SetHistLineWidth(2);

	gStyle->SetOptStat(0);
	//gStyle->SetOptStat("m");

	gStyle->SetTitleAlign(23);

	gROOT->ForceStyle();

};
