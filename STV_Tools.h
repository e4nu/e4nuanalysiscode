#ifndef STV_TOOLS_H
#define STV_TOOLS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TMath.h"
#include <TVector3.h>
#include <TLorentzVector.h>

class STV_Tools {

	private:

		double fkMiss;
		double fEMiss;
		double fPMissMinus;
		double fPMiss;		
		double fPt;
		double fDeltaAlphaT;
		double fDeltaPhiT;
		double fECal;
		double fEQE;
		double fQ2;										

	public:

		// Default constructor
		STV_Tools(TVector3 MuonVector,TVector3 ProtonVector, double MuonEnergy, double ProtonEnergy);

		// Default destructor
		//~STV_Tools(){}

		double ReturnkMiss();
		double ReturnEMiss();
		double ReturnPMissMinus();
		double ReturnPMiss();
		double ReturnPt();
		double ReturnDeltaAlphaT();
		double ReturnDeltaPhiT();
		double ReturnECal();
		double ReturnEQE();
		double ReturnQ2();		

};

#endif
