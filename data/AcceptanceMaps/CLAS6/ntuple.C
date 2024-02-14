#define ntuple_cxx
#include "ntuple.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH3D.h>

#include "../Constants.h"

using namespace std;
using namespace Constants;

void ntuple::Loop() {

	TFile* FileTH3D = new TFile("AcceptanceMap_"+particle+"_TH3D.root","recreate");

	TH3D* h3 = new TH3D("h3","",DeltaPNBins,DeltaPMin,DeltaPMax,YPTarENBins,YPTarEMin,YPTarEMin,ZvertENBins,ZvertEMin,ZvertEMax);


	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);	nbytes += nb;

		h3->SetBinContent(ID,IY,Zvert,Acc);

	}

	FileTH3D->cd();
	FileTH3D->Write();

}
