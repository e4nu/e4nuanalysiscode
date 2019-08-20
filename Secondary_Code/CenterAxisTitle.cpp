#include <TH1D.h>
#include <TH2D.h>

void CenterAxisTitle(TH1D* histo){

	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->CenterTitle();

};

void CenterAxisTitle(TH2D* histo2){

	histo2->GetXaxis()->CenterTitle();
	histo2->GetYaxis()->CenterTitle();

};
