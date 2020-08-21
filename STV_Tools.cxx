
//Class created by Afroditi Papadopoulou (apapadop@mit.edu)

// _________________________________________________________________________________________________________________________________________________

#ifndef STV_TOOLS_CXX
#define STV_TOOLS_CXX

#include <TMath.h>

#include "STV_Tools.h"

// __________________________________________________________________________________________________________________________________________________

STV_Tools::STV_Tools(TVector3 MuonVector,TVector3 ProtonVector, double MuonEnergy, double ProtonEnergy) {

	double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938272, NeutronMass_GeV = 0.939565; // GeV
	double DeltaM2 = TMath::Power(NeutronMass_GeV,2.) - TMath::Power(ProtonMass_GeV,2.);	
	double BE = 0.04; // GeV	
			
	// STV Calculation		
			
	TVector3 MuonVectorTrans;
	MuonVectorTrans.SetXYZ(MuonVector.X(),MuonVector.Y(),0.);
	double MuonVectorTransMag = MuonVectorTrans.Mag();
	
	TLorentzVector MuonLorentzVector(MuonVector,MuonEnergy);
//	double MuonKE = MuonEnergy - MuonMass_GeV;	
			
	TVector3 ProtonVectorTrans;
	ProtonVectorTrans.SetXYZ(ProtonVector.X(),ProtonVector.Y(),0.);
	double ProtonVectorTransMag = ProtonVectorTrans.Mag();
	
	TLorentzVector ProtonLorentzVector(ProtonVector,ProtonEnergy);
	double ProtonKE = ProtonEnergy - ProtonMass_GeV;		

	TVector3 PtVector = MuonVectorTrans + ProtonVectorTrans;

	fPt = PtVector.Mag();

	fDeltaAlphaT = TMath::ACos( (- MuonVectorTrans * PtVector) / ( MuonVectorTransMag * fPt ) ) * 180./TMath::Pi();
	if (fDeltaAlphaT > 180.) { fDeltaAlphaT -= 180.; }
	if (fDeltaAlphaT < 0.) { fDeltaAlphaT += 180.; }

	fDeltaPhiT = TMath::ACos( (- MuonVectorTrans * ProtonVectorTrans) / ( MuonVectorTransMag * ProtonVectorTransMag ) ) * 180./TMath::Pi();
	if (fDeltaPhiT > 180.) { fDeltaPhiT -= 180.; }
	if (fDeltaPhiT < 0.) { fDeltaPhiT += 180.; }

	// -------------------------------------------------------------------------------------------------------------------------

	// Calorimetric Energy Reconstruction

	fECal = MuonEnergy + ProtonKE + BE; // GeV

	// QE Energy Reconstruction

	double EQENum = 2 * (NeutronMass_GeV - BE) * MuonEnergy - (BE*BE - 2 * NeutronMass_GeV *BE + MuonMass_GeV * MuonMass_GeV + DeltaM2);
	double EQEDen = 2 * ( NeutronMass_GeV - BE - MuonEnergy + MuonVector.Mag() * MuonVector.CosTheta() );
	fEQE = EQENum / EQEDen;

	// Reconstructed Q2

	TLorentzVector nuLorentzVector(0.,0.,fECal,fECal);
	TLorentzVector qLorentzVector = nuLorentzVector - MuonLorentzVector;
	fQ2 = - qLorentzVector.Mag2();
	
	// -------------------------------------------------------------------------------------------------------------------------
	
	// Light Cone Variables
	
	TLorentzVector MissLorentzVector = nuLorentzVector - MuonLorentzVector - ProtonLorentzVector;
	
	fEMiss = TMath::Abs(MissLorentzVector.E());
	fPMiss = (MissLorentzVector.Vect()).Mag();	
	fPMissMinus = fEMiss - MissLorentzVector.Z();
	
	double fkMissNum = ( TMath::Power(fPt,2.) + TMath::Power(ProtonMass_GeV,2.) );
	double fkMissDen = ( fPMissMinus * (2*ProtonMass_GeV - fPMissMinus) );
	
	double fkMiss2 = TMath::Power(ProtonMass_GeV,2.) * fkMissNum / fkMissDen - TMath::Power(ProtonMass_GeV,2.); // Jackson's note				

	fkMiss = sqrt(fkMiss2);
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnkMiss() {

	return fkMiss;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnEMiss() {

	return fEMiss;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPMissMinus() {

	return fPMissMinus;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPMiss() {

	return fPMiss;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPt() {

	return fPt;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnDeltaAlphaT() {

	return fDeltaAlphaT;

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnDeltaPhiT() {

	return fDeltaPhiT;

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnECal() {

	return fECal;

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnEQE() {

	return fEQE;

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnQ2() {

	return fQ2;

}
// __________________________________________________________________________________________________________________________________________________

#endif
