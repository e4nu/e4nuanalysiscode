#include <TH3D.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

using namespace std;

double acceptance_c(double p, double cost, double phi, int id,TFile* file_acceptance) {

	//Redefinition of the phi angle
	
	int redef = -30;
	//int redef = 0;

	TH3D * e_acc;
	TH3D * e_gen;

	TH3D * p_acc;
	TH3D * p_gen;

        // Electron

	if (id == 11) {

		if (redef == -30) {

			e_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
			e_gen = (TH3D*)file_acceptance->Get("Generated Particles");
		}

		if (redef == 0) {

			e_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
			e_gen = (TH3D*)file_acceptance->Get("Generated Particles");
		}		
	
		//Find number of generated events

		double e_pbin_gen = e_gen->GetXaxis()->FindBin(p);
		double e_tbin_gen = e_gen->GetYaxis()->FindBin(cost);
		double e_phibin_gen = e_gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double e_num_gen = e_gen->GetBinContent(e_pbin_gen, e_tbin_gen, e_phibin_gen);

		//Find number of accepted events

		double e_pbin_acc = e_acc->GetXaxis()->FindBin(p);
		double e_tbin_acc = e_acc->GetYaxis()->FindBin(cost);
		double e_phibin_acc = e_acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double e_num_acc = e_acc->GetBinContent(e_pbin_acc, e_tbin_acc, e_phibin_acc);

		double e_acc_ratio = (double)e_num_acc / (double)e_num_gen;
		double e_acc_err = (double)sqrt(e_acc_ratio*(1-e_acc_ratio)) / (double)e_num_gen;

		return e_acc_ratio;
	}

        // Proton

	if (id == 2212) {

		if (redef == -30){

			p_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
			p_gen = (TH3D*)file_acceptance->Get("Generated Particles");
		}

		if (redef == 0){

			p_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
			p_gen = (TH3D*)file_acceptance->Get("Generated Particles");
		}

		//Find number of generated events

		double p_pbin_gen = p_gen->GetXaxis()->FindBin(p);
		double p_tbin_gen = p_gen->GetYaxis()->FindBin(cost);
		double p_phibin_gen = p_gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double p_num_gen = p_gen->GetBinContent(p_pbin_gen, p_tbin_gen, p_phibin_gen);

		//Find number of accepted events

		double p_pbin_acc = p_acc->GetXaxis()->FindBin(p);
		double p_tbin_acc = p_acc->GetYaxis()->FindBin(cost);
		double p_phibin_acc = p_acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double p_num_acc = p_acc->GetBinContent(p_pbin_acc, p_tbin_acc, p_phibin_acc);
		double p_acc_ratio = (double)p_num_acc / (double)p_num_gen;
		double p_acc_prr = (double)sqrt(p_acc_ratio*(1-p_acc_ratio)) / (double)p_num_gen;

		return p_acc_ratio;
	}

	return 0.;

}
