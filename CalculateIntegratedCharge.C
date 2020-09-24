#define CalculateIntegratedCharge_cxx
#include "CalculateIntegratedCharge.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

void CalculateIntegratedCharge::Loop() {

	 if (fChain == 0) return;

	 Long64_t nentries = fChain->GetEntriesFast();
	 Long64_t nbytes = 0, nb = 0;

	double IntegratedCharge = 0;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);	 nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		// ---------------------------------------------------------------------------------------------------------------------

		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// ------------------------------------------------------------------------------------------------------------------------

		double fTorusCurrent = 0;

		if((runnb>18283 && runnb<18289) || (runnb>18300 && runnb<18304) || (runnb>18317 && runnb<18329)) fTorusCurrent=750; //setting appropriate torrus magnet current
		else if ((runnb>18293 && runnb<18301) || (runnb>18305 && runnb<18317) || (runnb>18328 && runnb<18336))  fTorusCurrent=1500;
		else fTorusCurrent=2250;

		//if (fbeam_en == "1161" && fTorusCurrent > 760) { continue; }                                                              
                //if (fbeam_en == "1161" && fTorusCurrent < 760) { continue; }

		IntegratedCharge += q_l;

	} // End of the loop over the events

	// ---------------------------------------------------------------------------------------------------------------------------------

	cout << "Integrated charge = " << IntegratedCharge << endl;

}
