#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"

int ReadDataMacro()
{
    std::string data_file_name = "/storage/t3_data/tvjulia/e4nuanalysis_files/1pi_analysis/analised_data/e4nuanalysis_1pimanalysis_e_on_1000060120_4461MeV_4MaxBkgMult_clas6data.root";
    // TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    TFile *file_true_data = new TFile((data_file_name).c_str(), "ROOT");
    TTree *my_data_tree = (TTree *)file_true_data->Get("CLAS6Tree");

    if (!file_true_data)
    {
        std::cout << "ERROR: the " << file_true_data << " does not exist." << std::endl;
        return 0;
    }
    if (!my_data_tree)
        std::cout << " Tree does not exist " << std::endl;

    // OBSERVABLE DEFINITION:
    double RecoEvPion;

    my_gst_tree->SetBranchAddress("RecoEvPion", &RecoEvPion);

    long NEntries = my_data_tree->GetEntries();

    for (int j = 0; j < NEntries; ++j)
    {
        std::cout << RecoEvPion << std::endl;
    }

    return 0;
}