#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "Riostream.h"

using namespace std;

void acc(){

////	TTree *T = new TTree("ntuple","data from dat file","IDe:IYe:IXe:IDp:IYp:IXp:Ngen:Nrec:Acc");
//	TTree *T = new TTree("ntuple","data from dat file","IDe:IYe:IXe:IDp:IYp:IXp:Ngen:Nrec:Acc");
//	Long64_t nlines = T->ReadFile("H_mid_acc_map_6D.dat");
//	//printf("found %lld points\n",nlines);

//	TFile* file = new TFile("AcceptanceMaps.root","recreate");
//	file->Write();


	ifstream in;
//	in.open("H_mid_acc_map_6D.dat"); TFile *f = new TFile("AcceptanceMaps.root","RECREATE");
//	in.open("H_e4nu_e_map.dat"); TFile *f = new TFile("AcceptanceMaps_e.root","RECREATE");
	in.open("H_e4nu_p_map.dat"); TFile *f = new TFile("AcceptanceMaps_p.root","RECREATE");

	int ID, IY, IX, Zvert; 
	float Ngen, Nrec, Acc;

	Int_t nlines = 0;
	TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","ID:IY:IX:Zvert:Ngen:Nrec:Acc");

	while (1) {
		// 6D
//		in >> IDe >> IYe >> IXe >> IDp >> IYp >> IXp >> Ngen >> Nrec >> Acc;
//		if (!in.good()) break;
//		ntuple->Fill(IDe , IYe , IXe , IDp , IYp , IXp , Ngen , Nrec , Acc);

		// 4D for e/p
		in >> ID >> IY >> IX >> Zvert >> Ngen >> Nrec >> Acc;
		if (!in.good()) break;
		ntuple->Fill(ID , IY , IX , Zvert, Ngen , Nrec , Acc);

		nlines++;
	}
	printf(" found %d points\n",nlines);

	in.close();

	f->Write();

}
