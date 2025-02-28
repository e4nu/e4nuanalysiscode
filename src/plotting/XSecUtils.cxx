#include "plotting/XSecUtils.h"
#include "plotting/Systematics.h"
#include "TLegend.h"
#include <iomanip>
#include <filesystem>
#include <sstream>
#include <iostream>
#include <string>
#include "TGaxis.h"
#include "THStack.h"

using namespace e4nu;
using namespace e4nu::plotting;

void plotting::Plot1DXSec(std::vector<std::string> MC_files_name, std::string data_file_name,
                          std::string acceptance_file_name, std::string radcorr_file, std::string observable,
                          std::string title, std::string data_name, std::vector<std::string> model,
                          std::string input_MC_location, std::string input_data_location, std::string output_location,
                          std::string output_file_name, bool plot_data, std::map<string, double> systematic_map, std::map<std::string,std::vector<double>> cuts,
                          std::string analysis_id, bool store_root, bool log_scale)
{

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  std::vector<TFile *> files_true_MC;
  for (unsigned int id = 0; id < MC_files_name.size(); ++id)
  {
    files_true_MC.push_back(new TFile((input_MC_location + MC_files_name[id] + "_true.root").c_str(), "ROOT"));
    if (!files_true_MC[id])
    {
      std::cout << "ERROR: the " << input_MC_location << MC_files_name[id] << "_true.root does not exist." << std::endl;
      return;
    }
  }
  TFile *file_data = nullptr;
  if (plot_data)
  {
    file_data = new TFile((input_data_location + data_file_name + ".root").c_str(), "READ");
    if (!file_data)
    {
      std::cout << "ERROR: file data doesn't exist" << std::endl;
      return;
    }
  }
  else
  {
    std::cout << "Not plotting data" << std::endl;
  }

  TFile *file_acceptance = new TFile((output_location + acceptance_file_name + ".root").c_str(), "READ");
  TFile *file_radcorr = nullptr;
  if (radcorr_file != "")
    file_radcorr = new TFile((output_location + radcorr_file + ".root").c_str(), "READ");
  if (!file_data && plot_data)
  {
    std::cout << "ERROR: the " << input_data_location << data_file_name << ".root does not exist." << std::endl;
    return;
  }
  if (!file_acceptance)
  {
    std::cout << "ERROR: the " << output_location << acceptance_file_name << ".root does not exist." << std::endl;
    return;
  }

  TH1D *h_acceptance = (TH1D *)file_acceptance->Get("Acceptance");
  TH1D *h_acceptance_0 = (TH1D *)file_acceptance->Get("Acceptance_0");
  TH1D *h_acceptance_1 = (TH1D *)file_acceptance->Get("Acceptance_1");
  TH1D *h_acceptance_2 = (TH1D *)file_acceptance->Get("Acceptance_2");
  TH1D *h_acceptance_3 = (TH1D *)file_acceptance->Get("Acceptance_3");
  TH1D *h_acceptance_4 = (TH1D *)file_acceptance->Get("Acceptance_4");
  TH1D *h_acceptance_5 = (TH1D *)file_acceptance->Get("Acceptance_5");

  TH1D *h_radcorr = nullptr;
  if (file_radcorr)
  {
    h_radcorr = (TH1D *)file_radcorr->Get("Acceptance");
  }

  // Get Tree for main model
  TTree *tree_true = (TTree *)files_true_MC[0]->Get("MCCLAS6Tree");
  // Get configured energy, used for plotting
  double BeamE;
  tree_true->SetBranchAddress("BeamE", &BeamE);
  tree_true->GetEntry(0);

  // For submodels only total prediction is plotted
  std::vector<TTree *> tree_submodels;
  for (unsigned int id = 1; id < MC_files_name.size(); ++id)
  {
    tree_submodels.push_back((TTree *)files_true_MC[id]->Get("MCCLAS6Tree"));
    if (!tree_submodels[id - 1])
    {
      std::cout << "ERROR: the threes do not exist." << std::endl;
      return;
    }
  }

  TTree *tree_data = nullptr;
  if (plot_data)
  {
    tree_data = (TTree *)file_data->Get("CLAS6Tree");
    if (!tree_data)
    {
      std::cout << "ERROR: tree data doesn't exist" << std::endl;
      return;
    }
  }

  if (!h_acceptance )
  {
    std::cout << "ERROR: Acceptance is not defined" << std::endl;
    return;
  }
  if (!tree_true || (!tree_data && plot_data))
  {
    std::cout << "ERROR: the threes do not exist." << std::endl;
    return;
  }

  // Create histogram for total and total xsec per sector
  TH1D *hist_true = (TH1D *)h_acceptance->Clone();
  hist_true->SetName("MC_True");
  hist_true->Reset();

  TH1D *hist_true_0 = (TH1D *)h_acceptance_0->Clone();
  hist_true_0->SetName("MC_True_Sector_0");
  hist_true_0->Reset();
  TH1D *hist_true_1 = (TH1D *)h_acceptance_1->Clone();
  hist_true_1->SetName("MC_True_Sector_1");
  hist_true_1->Reset();
  TH1D *hist_true_2 = (TH1D *)h_acceptance_2->Clone();
  hist_true_2->SetName("MC_True_Sector_2");
  hist_true_2->Reset();
  TH1D *hist_true_3 = (TH1D *)h_acceptance_3->Clone();
  hist_true_3->SetName("MC_True_Sector_3");
  hist_true_3->Reset();
  TH1D *hist_true_4 = (TH1D *)h_acceptance_4->Clone();
  hist_true_4->SetName("MC_True_Sector_4");
  hist_true_4->Reset();
  TH1D *hist_true_5 = (TH1D *)h_acceptance_5->Clone();
  hist_true_5->SetName("MC_True_Sector_5");
  hist_true_5->Reset();

  // Breakdown histograms for total (all sectors only):
  TH1D *hist_true_QEL = (TH1D *)h_acceptance->Clone();
  hist_true_QEL->SetName("MC_True_QEL");
  hist_true_QEL->Reset();
  TH1D *hist_true_RES_Delta = (TH1D *)h_acceptance->Clone();
  hist_true_RES_Delta->SetName("MC_True_RES_Delta");
  hist_true_RES_Delta->Reset();
  TH1D *hist_true_RES = (TH1D *)h_acceptance->Clone();
  hist_true_RES->SetName("MC_True_RES");
  hist_true_RES->Reset();
  TH1D *hist_true_SIS = (TH1D *)h_acceptance->Clone();
  hist_true_SIS->SetName("MC_True_SIS");
  hist_true_SIS->Reset();
  TH1D *hist_true_MEC = (TH1D *)h_acceptance->Clone();
  hist_true_MEC->SetName("MC_True_MEC");
  hist_true_MEC->Reset();
  TH1D *hist_true_DIS = (TH1D *)h_acceptance->Clone();
  hist_true_DIS->SetName("MC_True_DIS");
  hist_true_DIS->Reset();

  // Same per model - only total prediction
  std::vector<TH1D *> hists_true_submodel;
  for (unsigned int id = 1; id < MC_files_name.size(); ++id)
  {
    hists_true_submodel.push_back((TH1D *)h_acceptance->Clone());
    hists_true_submodel[id - 1]->SetName(("MC_True_Model_" + std::to_string(id)).c_str());
    hists_true_submodel[id - 1]->Reset();
    hists_true_submodel[id - 1]->SetLineWidth(3);
  }

  // total and per sector
  TH1D *hist_data = nullptr, *hist_data_0 = nullptr, *hist_data_1 = nullptr, *hist_data_2 = nullptr, *hist_data_3 = nullptr, *hist_data_4 = nullptr, *hist_data_5 = nullptr;

  if (plot_data)
  {
    hist_data = (TH1D *)h_acceptance->Clone();
    hist_data->SetName("Data");
    hist_data->Reset();

    hist_data_0 = (TH1D *)h_acceptance_0->Clone();
    hist_data_0->SetName("Data_Sector_0");
    hist_data_0->Reset();
    hist_data_1 = (TH1D *)h_acceptance_1->Clone();
    hist_data_1->SetName("Data_Sector_1");
    hist_data_1->Reset();
    hist_data_2 = (TH1D *)h_acceptance_2->Clone();
    hist_data_2->Reset();
    hist_data_2->SetName("Data_Sector_2");
    hist_data_3 = (TH1D *)h_acceptance_3->Clone();
    hist_data_3->SetName("Data_Sector_3");
    hist_data_3->Reset();
    hist_data_4 = (TH1D *)h_acceptance_4->Clone();
    hist_data_4->SetName("Data_Sector_4");
    hist_data_4->Reset();
    hist_data_5 = (TH1D *)h_acceptance_5->Clone();
    hist_data_5->SetName("Data_Sector_5");
    hist_data_5->Reset();
  }

  std::vector<TTree *> trees = {tree_true};
  if (plot_data)
    trees.push_back(tree_data);

  std::vector<TH1D *> hists = {hist_true, hist_data, hist_true_0, hist_data_0,
                               hist_true_1, hist_data_1, hist_true_2, hist_data_2,
                               hist_true_3, hist_data_3, hist_true_4, hist_data_4,
                               hist_true_5, hist_data_5};

  unsigned int size_primary_trees = trees.size();
  unsigned int size_primary_hists = hists.size();
  // Adding total predictions for alternative models
  for (unsigned int id = 1; id < MC_files_name.size(); ++id)
  {
    trees.push_back(tree_submodels[id - 1]);
    hists.push_back(hists_true_submodel[id - 1]);
  }

  // If data is plot, the position of its trees and histograms is 1.
  // Otherwise we set it to a big number so it is ignored
  unsigned int id_data = 9999;
  if (plot_data)
    id_data = 1;

  // OBSERVABLE DEFINITION specific for MC
  long NEntries;
  bool IsBkg;
  int ElectronSector;
  bool QEL, RES, DIS, MEC;
  double MCNormalization, DataNormalization;
  std::vector<double> mc_norm;
  int resid;

  for (unsigned int i = 0; i < trees.size(); ++i)
  {
    NEntries = trees[i]->GetEntries();
    if (!trees[i])
      continue;
    trees[i]->SetBranchAddress("TotWeight", &TotWeight);
    trees[i]->SetBranchAddress("IsBkg", &IsBkg);
    trees[i]->SetBranchAddress("ECal", &ECal);
    trees[i]->SetBranchAddress("pfl_theta", &pfl_theta);
    trees[i]->SetBranchAddress("pfl_phi", &pfl_phi);
    trees[i]->SetBranchAddress("pfl", &pfl);
    trees[i]->SetBranchAddress("proton_mom", &proton_mom);
    trees[i]->SetBranchAddress("proton_theta", &proton_theta);
    trees[i]->SetBranchAddress("proton_phi", &proton_phi);
    trees[i]->SetBranchAddress("pim_mom", &pim_mom);
    trees[i]->SetBranchAddress("pim_theta", &pim_theta);
    trees[i]->SetBranchAddress("pim_phi", &pim_phi);
    trees[i]->SetBranchAddress("pip_mom", &pip_mom);
    trees[i]->SetBranchAddress("pip_theta", &pip_theta);
    trees[i]->SetBranchAddress("pip_phi", &pip_phi);
    trees[i]->SetBranchAddress("RecoW", &RecoW);
    trees[i]->SetBranchAddress("Recoq3", &Recoq3);
    trees[i]->SetBranchAddress("RecoQELEnu", &RecoQELEnu);
    trees[i]->SetBranchAddress("RecoXBJK", &RecoXBJK);
    trees[i]->SetBranchAddress("RecoQ2", &RecoQ2);
    trees[i]->SetBranchAddress("RecoEnergyTransfer", &RecoEnergyTransfer);
    trees[i]->SetBranchAddress("AlphaT", &AlphaT);
    trees[i]->SetBranchAddress("HadAlphaT", &HadAlphaT);
    trees[i]->SetBranchAddress("DeltaPT", &DeltaPT);
    trees[i]->SetBranchAddress("HadDeltaPT", &HadDeltaPT);
    trees[i]->SetBranchAddress("HadDeltaPTx", &HadDeltaPTx);
    trees[i]->SetBranchAddress("HadDeltaPTy", &HadDeltaPTy);
    trees[i]->SetBranchAddress("DeltaPhiT", &DeltaPhiT);
    trees[i]->SetBranchAddress("HadDeltaPhiT", &HadDeltaPhiT);
    trees[i]->SetBranchAddress("ElectronSector", &ElectronSector);
    trees[i]->SetBranchAddress("HadSystemMass", &HadSystemMass);
    trees[i]->SetBranchAddress("MissingEnergy", &MissingEnergy);
    trees[i]->SetBranchAddress("MissingAngle", &MissingAngle);
    trees[i]->SetBranchAddress("MissingMomentum", &MissingMomentum);
    trees[i]->SetBranchAddress("InferedNucleonMom", &InferedNucleonMom);
    trees[i]->SetBranchAddress("HadronsAngle", &HadronsAngle);
    trees[i]->SetBranchAddress("AdlerAngleThetaP", &AdlerAngleThetaP);
    trees[i]->SetBranchAddress("AdlerAnglePhiP", &AdlerAnglePhiP);
    trees[i]->SetBranchAddress("AdlerAngleThetaPi", &AdlerAngleThetaPi);
    trees[i]->SetBranchAddress("AdlerAnglePhiPi", &AdlerAnglePhiPi);
    trees[i]->SetBranchAddress("Angleqvshad", &Angleqvshad);
    trees[i]->SetBranchAddress("RecoEvPion", &RecoEvPion);
    trees[i]->SetBranchAddress("RecoWPion", &RecoWPion);
    trees[i]->SetBranchAddress("ElectronPT", &ElectronPT);
    trees[i]->SetBranchAddress("PionPT", &PionPT);

    // Only fill true info for the first model:
    if (i == 0)
    {
      trees[i]->SetBranchAddress("QEL", &QEL);
      trees[i]->SetBranchAddress("RES", &RES);
      trees[i]->SetBranchAddress("MEC", &MEC);
      trees[i]->SetBranchAddress("DIS", &DIS);
      trees[i]->SetBranchAddress("resid", &resid);
    }
    // Only second tree corresponds to data
    if (i != id_data)
    {
      trees[i]->SetBranchAddress("MCNormalization", &MCNormalization);
    }
    else
    {
      trees[i]->SetBranchAddress("DataNormalization", &DataNormalization);
    }

    for (int j = 0; j < NEntries; ++j)
    {
      trees[i]->GetEntry(j);
      double content = 0;
      double w = TotWeight;
      if (i != id_data && j == 0)
        mc_norm.push_back(MCNormalization);

      content = GetObservable(observable);

      // Check if passes cuts
      bool do_fill =true ;
      for (auto it = cuts.begin(); it != cuts.end(); it++)
      {
        double min = it->second[0] ;
        double max = it->second[1] ;
        if( GetObservable(it->first) < min || GetObservable(it->first) > max ) {
          do_fill = false;
          continue;
        }
      }
      if( !do_fill ) continue ;

      unsigned int id_hist = i;
      // Fill the per Sector  histogram. Only for primary model
      if (hists[size_primary_trees * (ElectronSector + 1) + i])
      {
        if (i < size_primary_trees)
        {
          hists[size_primary_trees * (ElectronSector + 1) + i]->Fill(content, w);
        }
      }
      if (i > size_primary_trees - 1)
        id_hist = size_primary_hists + (i - size_primary_trees);

      if (hists[id_hist])
      {
        hists[id_hist]->Fill(content, w);
        hists[id_hist]->SetLineWidth(3);
      }

      if (i == 0)
      {
        if (QEL)
          hist_true_QEL->Fill(content, w);
        if (RES)
        {
          if (resid == 0)
            hist_true_RES_Delta->Fill(content, w);
          else
            hist_true_RES->Fill(content, w);
        }
        if (DIS)
        {
          if (RecoW < 1.7)
            hist_true_SIS->Fill(content, w);
          else
            hist_true_DIS->Fill(content, w);
        }
        if (MEC)
          hist_true_MEC->Fill(content, w);
      }
    }
  }

  // Store uncorrected data
  TH1D *hist_data_uncorr = nullptr, *hist_data_uncorr_0 = nullptr, *hist_data_uncorr_1 = nullptr, *hist_data_uncorr_2 = nullptr, *hist_data_uncorr_3 = nullptr, *hist_data_uncorr_4 = nullptr, *hist_data_uncorr_5 = nullptr;
  // Store corrected for acceptance but not for radiation
  TH1D *hist_data_uncorrrad = nullptr;
  // Corr event rate
  TH1D *hist_data_correventrate = nullptr, *hist_data_correventrate_0 = nullptr, *hist_data_correventrate_1 = nullptr, *hist_data_correventrate_2 = nullptr, *hist_data_correventrate_3 = nullptr, *hist_data_correventrate_4 = nullptr, *hist_data_correventrate_5 = nullptr;
  // Event rate with Systematics
  TH1D *hist_data_correventrate_wsyst = nullptr, *hist_data_correventrate_wsyst_0 = nullptr, *hist_data_correventrate_wsyst_1 = nullptr, *hist_data_correventrate_wsyst_2 = nullptr, *hist_data_correventrate_wsyst_3 = nullptr, *hist_data_correventrate_wsyst_4 = nullptr, *hist_data_correventrate_wsyst_5 = nullptr;
  // Event rate with Systematics
  TH1D *hist_xsec_wsyst = nullptr, *hist_xsec_wsyst_0 = nullptr, *hist_xsec_wsyst_1 = nullptr, *hist_xsec_wsyst_2 = nullptr, *hist_xsec_wsyst_3 = nullptr, *hist_xsec_wsyst_4 = nullptr, *hist_xsec_wsyst_5 = nullptr;

  // Store data event rate before acceptance correction:
  if (plot_data && hist_data)
  {
    hist_data_uncorr = (TH1D *)hist_data->Clone();
    hist_data_uncorr->SetName("Uncorrected Data");
    hist_data_uncorr_0 = (TH1D *)hist_data_0->Clone();
    hist_data_uncorr_0->SetName("Uncorrected Data Sector  0");
    hist_data_uncorr_1 = (TH1D *)hist_data_1->Clone();
    hist_data_uncorr_1->SetName("Uncorrected Data Sector  1");
    hist_data_uncorr_2 = (TH1D *)hist_data_2->Clone();
    hist_data_uncorr_2->SetName("Uncorrected Data Sector  2");
    hist_data_uncorr_3 = (TH1D *)hist_data_3->Clone();
    hist_data_uncorr_3->SetName("Uncorrected Data Sector  3");
    hist_data_uncorr_4 = (TH1D *)hist_data_4->Clone();
    hist_data_uncorr_4->SetName("Uncorrected Data Sector  4");
    hist_data_uncorr_5 = (TH1D *)hist_data_5->Clone();
    hist_data_uncorr_5->SetName("Uncorrected Data Sector  5");

    // Normaize for bin size
    NormalizeHist(hist_data_uncorr, 1);
    NormalizeHist(hist_data_uncorr_0, 1);
    NormalizeHist(hist_data_uncorr_1, 1);
    NormalizeHist(hist_data_uncorr_2, 1);
    NormalizeHist(hist_data_uncorr_3, 1);
    NormalizeHist(hist_data_uncorr_4, 1);
    NormalizeHist(hist_data_uncorr_5, 1);

    // Correct data for detector acceptance :
    // Notice this step already propagates the error from the acceptance to the corrected event rate
    // (Err_Corr_eventrate)^2 = (Err_Raw_EventRate)^2 * (Acc)^2 + (Raw_EventRate)^2 * (Err_Acc)^2
    CorrectData(hist_data, h_acceptance);
    CorrectData(hist_data_0, h_acceptance_0);
    CorrectData(hist_data_1, h_acceptance_1);
    CorrectData(hist_data_2, h_acceptance_2);
    CorrectData(hist_data_3, h_acceptance_3);
    CorrectData(hist_data_4, h_acceptance_4);
    CorrectData(hist_data_5, h_acceptance_5);

    hist_data_uncorrrad = (TH1D *)hist_data->Clone();
    hist_data_uncorrrad->SetName("Corrected for acceptance before rad corr data");
    NormalizeHist(hist_data_uncorrrad, 1);

    //Apply radiative correction
    if (h_radcorr)
    {
      CorrectData(hist_data, h_radcorr);
      CorrectData(hist_data_0, h_radcorr);
      CorrectData(hist_data_1, h_radcorr);
      CorrectData(hist_data_2, h_radcorr);
      CorrectData(hist_data_3, h_radcorr);
      CorrectData(hist_data_4, h_radcorr);
      CorrectData(hist_data_5, h_radcorr);
    }

    hist_data_correventrate = (TH1D *)hist_data->Clone();
    hist_data_correventrate->SetName("Corrected_Event_Rate_Data");
    hist_data_correventrate_0 = (TH1D *)hist_data_0->Clone();
    hist_data_correventrate_0->SetName("Corrected_Event_Rate_Data_Sector_0");
    hist_data_correventrate_1 = (TH1D *)hist_data_1->Clone();
    hist_data_correventrate_1->SetName("Corrected_Event_Rate_Data_Sector_1");
    hist_data_correventrate_2 = (TH1D *)hist_data_2->Clone();
    hist_data_correventrate_2->SetName("Corrected_Event_Rate_Data_Sector_2");
    hist_data_correventrate_3 = (TH1D *)hist_data_3->Clone();
    hist_data_correventrate_3->SetName("Corrected_Event_Rate_Data_Sector_3");
    hist_data_correventrate_4 = (TH1D *)hist_data_4->Clone();
    hist_data_correventrate_4->SetName("Corrected_Event_Rate_Data_Sector_4");
    hist_data_correventrate_5 = (TH1D *)hist_data_5->Clone();
    hist_data_correventrate_5->SetName("Corrected_Event_Rate_Data_Sector_5");

    // Store event rate with systematics
    hist_data_correventrate_wsyst = (TH1D *)hist_data->Clone();
    hist_data_correventrate_wsyst->SetName("Corrected Event Rate with Systematics Data");
    hist_data_correventrate_wsyst_0 = (TH1D *)hist_data_0->Clone();
    hist_data_correventrate_wsyst_0->SetName("Corrected Event Rate with Systematics Sector 0");
    hist_data_correventrate_wsyst_1 = (TH1D *)hist_data_1->Clone();
    hist_data_correventrate_wsyst_1->SetName("Corrected Event Rate with Systematics Sector  1");
    hist_data_correventrate_wsyst_2 = (TH1D *)hist_data_2->Clone();
    hist_data_correventrate_wsyst_2->SetName("Corrected Event Rate with Systematics Sector  2");
    hist_data_correventrate_wsyst_3 = (TH1D *)hist_data_3->Clone();
    hist_data_correventrate_wsyst_3->SetName("Corrected Event Rate with Systematics Sector  3");
    hist_data_correventrate_wsyst_4 = (TH1D *)hist_data_4->Clone();
    hist_data_correventrate_wsyst_4->SetName("Corrected Event Rate with Systematics Sector  4");
    hist_data_correventrate_wsyst_5 = (TH1D *)hist_data_5->Clone();
    hist_data_correventrate_wsyst_5->SetName("Corrected Event Rate with Systematics Sector  5");

    // Normaize by bin width
    NormalizeHist(hist_data_correventrate_wsyst, 1);
    NormalizeHist(hist_data_correventrate_wsyst_0, 1);
    NormalizeHist(hist_data_correventrate_wsyst_1, 1);
    NormalizeHist(hist_data_correventrate_wsyst_2, 1);
    NormalizeHist(hist_data_correventrate_wsyst_3, 1);
    NormalizeHist(hist_data_correventrate_wsyst_4, 1);
    NormalizeHist(hist_data_correventrate_wsyst_5, 1);

    // Normaize by bin width
    NormalizeHist(hist_data_correventrate, 1);
    NormalizeHist(hist_data_correventrate_0, 1);
    NormalizeHist(hist_data_correventrate_1, 1);
    NormalizeHist(hist_data_correventrate_2, 1);
    NormalizeHist(hist_data_correventrate_3, 1);
    NormalizeHist(hist_data_correventrate_4, 1);
    NormalizeHist(hist_data_correventrate_5, 1);

    // Normalize to cross-section
    NormalizeHist(hist_data, DataNormalization);
    NormalizeHist(hist_data_0, DataNormalization);
    NormalizeHist(hist_data_1, DataNormalization);
    NormalizeHist(hist_data_2, DataNormalization);
    NormalizeHist(hist_data_3, DataNormalization);
    NormalizeHist(hist_data_4, DataNormalization);
    NormalizeHist(hist_data_5, DataNormalization);

    // Add Systematics
    // 1 - Acceptance model dependence (already included)
    // 2 - Sector Sector Variation
    // 3 - Relative uncertanties from configuration

    // Adding Acceptance correction systematics from model dependence
    TH1D *hist_syst_acc = systematics::AddSystematic(*hist_data, *h_acceptance);

    TCanvas *cacc = new TCanvas("cacc", "cacc", 800, 600);
    hist_syst_acc->Draw("hist");
    cacc->SaveAs((output_location + "/XSecPerSector/" + output_file_name + "_syst_accmodel_" + observable + ".root").c_str());
    delete cacc;

    // Add sector variation ERROR. Store relative error in histogram
    // We use the bkg substracted, eff corrected distributions for the calculation
    TH1D *hist_syst_sector = systematics::SectorVariationError(*hist_data, {hist_data_0, hist_data_1, hist_data_2, hist_data_3, hist_data_4, hist_data_5});
    systematics::SectorVariationError(*hist_data_correventrate_wsyst, {hist_data_0, hist_data_1, hist_data_2, hist_data_3, hist_data_4, hist_data_5});

    TCanvas *csect = new TCanvas("csect", "csect", 800, 600);
    hist_syst_sector->Draw("hist");
    csect->SaveAs((output_location + "/XSecPerSector/" + output_file_name + "_syst_persector_" + observable + ".root").c_str());
    delete csect;

    // adding systematics from systematic map. Relative systematic added to all bins
    for (auto it = systematic_map.begin(); it != systematic_map.end(); ++it)
    {
      std::cout << " Adding " << it->second << " % systematic on " << it->first << std::endl;
      systematics::AddSystematic(*hist_data, it->second, it->first);
      systematics::AddSystematic(*hist_data_correventrate_wsyst, it->second, it->first);
    }

    // Hard coding some well known systematics
    systematics::AddSystematic(*hist_data, 5, "Radiative");
    systematics::AddSystematic(*hist_data, 3, "Normalization");
    systematics::AddSystematic(*hist_data, 1, "AnglDependence");
    systematics::AddSystematic(*hist_data, 1, "MaxMultiplicity");

    systematics::AddSystematic(*hist_data_correventrate_wsyst, 5, "Radiative");
    systematics::AddSystematic(*hist_data_correventrate_wsyst, 3, "Normalization");
    systematics::AddSystematic(*hist_data_correventrate_wsyst, 1, "AnglDependence");
    systematics::AddSystematic(*hist_data_correventrate_wsyst, 1, "MaxMultiplicity");

  } // end if data

  // Normalize MC to cross-section
  for (unsigned int id = 0; id < hists_true_submodel.size(); ++id)
  {
    NormalizeHist(hists_true_submodel[id], mc_norm[id + 1]);
    StandardFormat(hists_true_submodel[id], title, kBlack, 2 + id, observable);
  }

  NormalizeHist(hist_true, mc_norm[0]);
  NormalizeHist(hist_true_0, mc_norm[0]);
  NormalizeHist(hist_true_1, mc_norm[0]);
  NormalizeHist(hist_true_2, mc_norm[0]);
  NormalizeHist(hist_true_3, mc_norm[0]);
  NormalizeHist(hist_true_4, mc_norm[0]);
  NormalizeHist(hist_true_5, mc_norm[0]);
  NormalizeHist(hist_true_QEL, mc_norm[0]);
  NormalizeHist(hist_true_RES, mc_norm[0]);
  NormalizeHist(hist_true_RES_Delta, mc_norm[0]);
  NormalizeHist(hist_true_SIS, mc_norm[0]);
  NormalizeHist(hist_true_MEC, mc_norm[0]);
  NormalizeHist(hist_true_DIS, mc_norm[0]);

  // Store histograms for plotting
  std::vector<TH1D> mc_hists = {*hist_true};
  std::vector<TH1D *> mc_hists_xsec = {hist_true};
  for (unsigned int id = 0; id < hists_true_submodel.size(); ++id)
  {
    mc_hists.push_back(*hists_true_submodel[id]);
    mc_hists_xsec.push_back(hists_true_submodel[id]);
  }

  std::vector<TH1D> breakdown = {*hist_true_QEL, *hist_true_RES_Delta, *hist_true_RES, *hist_true_SIS, *hist_true_MEC, *hist_true_DIS};
  std::vector<TH1D *> breakdown_xsec = {hist_true_QEL, hist_true_RES_Delta, hist_true_RES, hist_true_SIS, hist_true_MEC, hist_true_DIS};
  std::vector<TH1D *> mc_per_sector = {hist_true_0, hist_true_1, hist_true_2, hist_true_3, hist_true_4, hist_true_5};
  std::vector<TH1D *> data_per_sector = {hist_data_0, hist_data_1, hist_data_2, hist_data_3, hist_data_4, hist_data_5};
  std::vector<TH1D *> data_per_sector_uncorr = {hist_data_uncorr_0, hist_data_uncorr_1, hist_data_uncorr_2, hist_data_uncorr_3, hist_data_uncorr_4, hist_data_uncorr_5};
  // Plot Total, XSector, Legend
  if (plot_data)
  {
    plotting::PlotEventRate(hist_data_uncorr, observable, title, data_name, input_data_location, output_location,
                            output_file_name + "_raw_event_rate", analysis_id, store_root);

    plotting::PlotEventRatePerSector(data_per_sector, observable, title, data_name, input_data_location, output_location,
                                     output_file_name + "_event_rate_corracc", analysis_id, store_root);

    plotting::PlotComparisonDataNormalized(mc_hists, breakdown, hist_data_correventrate, observable, title, data_name, model, input_MC_location,
                                           input_data_location, output_location, output_file_name + "_normalized_to_data_with_breakdown", systematic_map, true, analysis_id, store_root);

    plotting::PlotComparisonDataNormalized(mc_hists, breakdown, hist_data_correventrate, observable, title, data_name, model, input_MC_location,
                                           input_data_location, output_location, output_file_name + "_normalized_to_data_no_breakdown", systematic_map, false, analysis_id, store_root);

    plotting::PlotComparisonDataNormalized(mc_hists, breakdown, hist_data_correventrate_wsyst, observable, title, data_name, model, input_MC_location,
                                           input_data_location, output_location, output_file_name + "_normalized_to_data_wsyst_no_breakdown", systematic_map, false, analysis_id, store_root);

    plotting::PlotComparisonDataNormalized(mc_hists, breakdown, hist_data_correventrate_wsyst, observable, title, data_name, model, input_MC_location,
                                           input_data_location, output_location, output_file_name + "_normalized_to_data_wsyst_with_breakdown", systematic_map, true, analysis_id, store_root);

    plotting::PlotXsecDataTotal(hist_data, observable, title, data_name, input_data_location, output_location, output_file_name + "MYTEST",
                                systematic_map, analysis_id, store_root);

  }

  plotting::PlotTotalXSec(mc_hists_xsec, breakdown_xsec, hist_data, observable, title, data_name, model, input_MC_location,
                          input_data_location, output_location, output_file_name + "_with_breakdown", systematic_map, true, analysis_id, store_root, log_scale);

  plotting::PlotTotalXSec(mc_hists_xsec, breakdown_xsec, hist_data, observable, title, data_name, model, input_MC_location,
                          input_data_location, output_location, output_file_name + "_no_breakdown", systematic_map, false, analysis_id, store_root, log_scale);

  plotting::PlotPerSector(mc_per_sector, data_per_sector, observable, title, data_name, model, input_MC_location,
                          input_data_location, output_location, output_file_name, systematic_map, analysis_id, store_root);

  plotting::PlotLegend(mc_hists_xsec, breakdown_xsec, hist_data, observable, data_name, model, output_location, output_file_name, store_root);

}

void plotting::PlotXsecDataTotal(TH1D *data, std::string observable, std::string title, std::string data_name,
                                 std::string input_data_location, std::string output_location,
                                 std::string output_file_name, std::map<string, double> systematic_map,
                                 std::string analysis_id, bool store_root, bool log_scale)
{
  if( log_scale ) gPad->SetLogy();
  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 700, 500);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Format plots
  if (data)
  {
    StandardFormat(data, title, kBlack, 8, observable);
    data->SetLineStyle(1);
  }

  if (data)
    data->Draw(" err ");

  std::string output_name = output_file_name + "_onlydata_dxsec_d" + observable;

  std::filesystem::path totalxsec_path{(output_location + "/TotalXSec/").c_str()};
  if (!std::filesystem::exists(totalxsec_path))
    std::filesystem::create_directory(totalxsec_path);

  if (store_root)
  {
    TFile root_file((output_location + "/TotalXSec/" + output_name + ".root").c_str(), "recreate");
    data->Write();
    c1->Write();
  }
  c1->SaveAs((output_location + "/TotalXSec/" + output_name + ".pdf").c_str());
  delete c1;
}

void plotting::PlotComparisonDataNormalized(std::vector<TH1D> mc_hists, std::vector<TH1D> breakdown,
                                            TH1D *data, std::string observable, std::string title, std::string data_name, std::vector<std::string> model,
                                            std::string input_MC_location, std::string input_data_location, std::string output_location,
                                            std::string output_file_name, std::map<string, double> systematic_map, bool show_breakdown,
                                            std::string analysis_id, bool store_root, bool log_scale)
{
  if( log_scale ) gPad->SetLogy();
  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 700, 500);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);

  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Print out integral for debugging
  double data_integral = 0;
  if (data)
    data_integral = data->Integral("width");
  double mc_integral = mc_hists[0].Integral("width");

  // Normalize MC to data
  mc_hists[0].Scale(data_integral / mc_integral);
  for (unsigned int i = 1; i < mc_hists.size(); ++i)
    mc_hists[i].Scale(data_integral / mc_integral);
  if (show_breakdown)
  {
    for (unsigned int i = 0; i < breakdown.size(); ++i)
      breakdown[i].Scale(data_integral / mc_integral);
  }

  // Find absolute y max
  std::vector<TH1D *> temp_check = {&mc_hists[0]};
  if (data)
    temp_check.push_back(data);

  for (unsigned int id = 1; id < mc_hists.size(); ++id)
  {
    temp_check.push_back(&mc_hists[id]);
  }
  double y_max_total = GetMaximum(temp_check);

  // Format plots
  if (data)
  {
    StandardFormat(data, title, kBlack, 8, observable, y_max_total);
    data->SetLineStyle(1);
  }

  StandardFormat(&mc_hists[0], title, kBlack, 1, observable, log_scale, y_max_total, "Counts/Bin Width");
  if (breakdown.size() == 6)
  {
    StandardFormat(&breakdown[0], title, kGreen + 1 , 1, observable, log_scale, y_max_total, "Counts/Bin Width");
    StandardFormat(&breakdown[1], title, kRed - 4, 1, observable, log_scale, y_max_total, "Counts/Bin Width");
    StandardFormat(&breakdown[2], title, kMagenta - 3, 1, observable, log_scale, y_max_total, "Counts/Bin Width");
    StandardFormat(&breakdown[3], title, kCyan + 1, 1, observable, log_scale, y_max_total, "Counts/Bin Width");
    StandardFormat(&breakdown[4], title, kOrange, 1, observable, log_scale, y_max_total, "Counts/Bin Width");
    StandardFormat(&breakdown[5], title, kBlue, 1, observable, log_scale, y_max_total, "Counts/Bin Width");
  }

  // Draw total xsec (all sectors):
  mc_hists[0].Draw("hist");
  for (unsigned int i = 1; i < mc_hists.size(); ++i)
    mc_hists[i].Draw("hist same");
  if (show_breakdown)
  {
    for (unsigned int i = 0; i < breakdown.size(); ++i)
      breakdown[i].Draw("hist same");
  }

  if (data)
    data->Draw(" err same ");

  std::string output_name = output_file_name + "_dxsec_d" + observable;

  std::filesystem::path totalxsec_path{(output_location + "/TotalXSec/").c_str()};
  if (!std::filesystem::exists(totalxsec_path))
    std::filesystem::create_directory(totalxsec_path);

  if (store_root)
  {
    TFile root_file((output_location + "/TotalXSec/" + output_name + ".root").c_str(), "recreate");
    data->Write();
    for (unsigned int i = 0; i < mc_hists.size(); ++i)
      mc_hists[i].Write();
    for (unsigned int i = 0; i < breakdown.size(); ++i)
      breakdown[i].Write();
    c1->Write();
  }
  c1->SaveAs((output_location + "/TotalXSec/" + output_name + ".pdf").c_str());
  delete c1;
}

void plotting::PlotTotalXSec(std::vector<TH1D *> mc_hists, std::vector<TH1D *> breakdown, TH1D *data, std::string observable,
                             std::string title, std::string data_name, std::vector<std::string> model,
                             std::string input_MC_location, std::string input_data_location, std::string output_location,
                             std::string output_file_name, std::map<string, double> systematic_map, bool show_breakdown,
                             std::string analysis_id, bool store_root, bool log_scale)
{
  TCanvas *c1 = new TCanvas("c1", "Canvas with Two Pads", 800, 600);

  // Create the upper pad, taking the top half of the canvas
  TPad *pad1 = new TPad("pad1", "Top Pad", 0, 0.35, 1, 1);
  pad1->SetBottomMargin(0.025); // Minimize gap between pads
  pad1->SetLeftMargin(0.15);
  pad1->Draw();

  // Create the lower pad, taking the bottom half of the canvas
  TPad *sub_pad = new TPad("sub_pad", "Bottom Pad", 0, 0, 1, 0.35);
  sub_pad->SetTopMargin(0.05); // Minimize gap between pads
  sub_pad->SetBottomMargin(0.4); // Leave space for the x-axis labels
  sub_pad->SetBottomMargin(0.35);
  sub_pad->SetLeftMargin(0.15);
  sub_pad->Draw();

  if( log_scale ) pad1->SetLogy();

  // Fill xsec canvas
  pad1->cd();

  // Set correct plotting style
  // Find absolute y max
  std::vector<TH1D *> temp_check = {mc_hists[0]};
  if (data)
    temp_check.push_back(data);

  for (unsigned int id = 1; id < mc_hists.size(); ++id)
  {
    temp_check.push_back(mc_hists[id]);
  }

  // Format plots
  if (data)
  {
    StandardFormat(data, title, kBlack, 8, observable, log_scale);
    data->SetLineStyle(1);
  }

  if (breakdown.size() == 6)
  {
    StandardFormat(breakdown[0], title, ColorBlindPalette(0), 1, observable, log_scale);
    StandardFormat(breakdown[1], title, ColorBlindPalette(1), 1, observable, log_scale);
    StandardFormat(breakdown[2], title, ColorBlindPalette(9), 1, observable, log_scale);
    StandardFormat(breakdown[3], title, ColorBlindPalette(2), 1, observable, log_scale);
    StandardFormat(breakdown[4], title, ColorBlindPalette(6), 1, observable, log_scale);
    StandardFormat(breakdown[5], title, ColorBlindPalette(3), 1, observable, log_scale);

    breakdown[0]->SetFillColorAlpha(ColorBlindPalette(0),0.6);
    breakdown[1]->SetFillColorAlpha(ColorBlindPalette(1),0.6);
    breakdown[2]->SetFillColorAlpha(ColorBlindPalette(9),0.6);
    breakdown[3]->SetFillColorAlpha(ColorBlindPalette(2),0.6);
    breakdown[4]->SetFillColorAlpha(ColorBlindPalette(6),0.6);
    breakdown[5]->SetFillColorAlpha(ColorBlindPalette(3),0.6);

  }

  double max_hist = plotting::GetMaximum(mc_hists);
  double min_hist = 0 ;
  if( log_scale ) min_hist = 1E-4 ;
  for (unsigned int i = 0; i < mc_hists.size(); ++i){
    StandardFormat(mc_hists[i], title, kBlack, i+1, observable, log_scale);
  }
  // StandardFormat(mc_hists[0], title, kBlack, 1, observable, log_scale);
  // Remove top plot label
  mc_hists[0]->GetXaxis()->SetLabelSize(0.);
  mc_hists[0]->GetXaxis()->SetTitleSize(0.);
  mc_hists[0]->GetYaxis()->SetTitleOffset(0.8);
  mc_hists[0]->SetMarkerSize(0);
  mc_hists[0]->GetYaxis()->SetRangeUser(min_hist, max_hist);

  // Fill tstack Plot
  auto hs = new THStack("hs","");
  hs->Add(breakdown[0]);
  hs->Add(breakdown[1]);
  hs->Add(breakdown[2]);
  hs->Add(breakdown[3]);
  hs->Add(breakdown[4]);
  hs->Add(breakdown[5]);

  // print
  mc_hists[0]->Draw("hist err ");
  hs->Draw("hist err same");
  for (unsigned int i = 0; i < mc_hists.size(); ++i) {
      mc_hists[i]->Draw("hist err same");
      mc_hists[i]->SetMarkerSize(0);
  }

  if (data)
      data->Draw("err same");

  // Plot Ratio
  sub_pad->cd();
  TH1D * data_ratio = (TH1D*)data->Clone();
  std::vector<TH1D*> mc_ratio ;
  for (unsigned int i = 0; i < mc_hists.size(); ++i){
    mc_ratio.push_back((TH1D*)mc_hists[i]->Clone());
    StandardFormat(mc_ratio[i], title, kBlack, i+1, observable, log_scale);
  }
  StandardFormat(data_ratio, title, kBlack, 8, observable, log_scale);
  data_ratio->SetLineStyle(1);
  data_ratio->Divide(data);
  data_ratio->GetXaxis()->SetLabelSize(0.1);
  data_ratio->GetXaxis()->SetTitleSize(0.15);
  data_ratio->GetYaxis()->SetLabelSize(0.1);
  data_ratio->GetYaxis()->SetTitleSize(0.13);
  data_ratio->GetYaxis()->SetTitleOffset(0.4);
  data_ratio->SetMinimum(0);
  data_ratio->GetYaxis()->SetMaxDigits(5);
  data_ratio->GetYaxis()->SetTitle("Pred / Data ");


  for (unsigned int i = 0; i < mc_hists.size(); ++i){
    mc_ratio[i]->Divide(data);
    mc_ratio[i]->GetYaxis()->SetTitle("Pred / Data ");
  }
  // Find correct range
  double max_ratio = GetMaximum(mc_ratio)*(1-0.35);
  double min_ratio = GetMinimum(mc_ratio)*(1-0.15);
  data_ratio->GetYaxis()->SetRangeUser(min_ratio,max_ratio);
  if( log_scale ) sub_pad->SetLogy();
  data_ratio->Draw("err");
  for (unsigned int i = 0; i < mc_hists.size(); ++i){
    mc_ratio[i]->GetYaxis()->SetRangeUser(min_ratio,max_ratio);
    mc_ratio[i]->Draw("hist err same");
  }

  std::string output_name = output_file_name + "_dxsec_d" + observable ;
  if (store_root)
  {
    TFile root_file((output_location + "/TotalXSec/" + output_name + ".root").c_str(), "recreate");
    if (data)
      data->Write();
    for (unsigned int i = 0; i < mc_hists.size(); ++i)
      mc_hists[i]->Write();
    for (unsigned int i = 0; i < breakdown.size(); ++i)
      breakdown[i]->Write();
    c1->Write();
  }
  // Return to the canvas for further customization if needed
  c1->cd();


  std::filesystem::path totalxsec_path{(output_location + "/TotalXSec/").c_str()};
  if (!std::filesystem::exists(totalxsec_path))
    std::filesystem::create_directory(totalxsec_path);

  // Print out integral for debugging
  double data_integral;
  if (data)
    data_integral = data->Integral("width");
  double mc_integral = mc_hists[0]->Integral("width");
  if (data)
  {
    std::cout << " Total integrated cross section (data) " << data_integral << std::endl;
  }
  std::cout << " Total integrated cross section (mc) " << mc_integral << std::endl;

  c1->SaveAs((output_location + "/TotalXSec/" + output_name + ".pdf").c_str());
  delete c1;

}


void plotting::PlotTotal2DXSec(std::vector<TH2D *> mc_hists, std::vector<TH2D *> breakdown, TH2D *data, std::string x_observable, std::string y_observable,
                             std::vector<double> & y_cuts, std::string title, std::string data_name, std::vector<std::string> model,
                             std::string input_MC_location, std::string input_data_location, std::string output_location,
                             std::string output_file_name, std::map<string, double> systematic_map, bool show_breakdown,
                             std::string analysis_id, bool store_root, bool log_scale)
{

  TCanvas *c_sector = new TCanvas("c_sector", "c_sector", 200, 10, 700, 500);
  c_sector->cd();
  TPad *pad_sector = new TPad("pad1", "", 0, 0, 1, 1);
  pad_sector->Draw();
  pad_sector->cd();
  pad_sector->SetBottomMargin(0.15);
  pad_sector->SetLeftMargin(0.15);
  pad_sector->SetRightMargin(0.15);
  pad_sector->Divide(2, 1);

  // Set correct plotting style
  // Find absolute y max
  std::vector<TH2D *> temp_check = {mc_hists[0]};
  if (data)
    temp_check.push_back(data);

  for (unsigned int id = 1; id < mc_hists.size(); ++id)
  {
    temp_check.push_back(mc_hists[id]);
  }

  // Format plots
  if (data)
  {
    StandardFormat(data, title, kBlack, 8, x_observable, y_observable);
    data->SetLineStyle(1);
  }

  StandardFormat(mc_hists[0], title, kBlack, 1, x_observable, y_observable);

  if (breakdown.size() == 6)
  {
    StandardFormat(breakdown[0], title, ColorBlindPalette(0), 1, x_observable, y_observable);
    StandardFormat(breakdown[1], title, ColorBlindPalette(1), 1, x_observable, y_observable);
    StandardFormat(breakdown[2], title, ColorBlindPalette(9), 1, x_observable, y_observable);
    StandardFormat(breakdown[3], title, ColorBlindPalette(2), 1, x_observable, y_observable);
    StandardFormat(breakdown[4], title, ColorBlindPalette(6), 1, x_observable, y_observable);
    StandardFormat(breakdown[5], title, ColorBlindPalette(3), 1, x_observable, y_observable);

    breakdown[0]->SetFillColorAlpha(ColorBlindPalette(0),0.6);
    breakdown[1]->SetFillColorAlpha(ColorBlindPalette(1),0.6);
    breakdown[2]->SetFillColorAlpha(ColorBlindPalette(9),0.6);
    breakdown[3]->SetFillColorAlpha(ColorBlindPalette(2),0.6);
    breakdown[4]->SetFillColorAlpha(ColorBlindPalette(6),0.6);
    breakdown[5]->SetFillColorAlpha(ColorBlindPalette(3),0.6);

  }

  // Fill tstack Plot
  auto hs = new THStack("hs","");
  hs->Add(breakdown[0]);
  hs->Add(breakdown[1]);
  hs->Add(breakdown[2]);
  hs->Add(breakdown[3]);
  hs->Add(breakdown[4]);
  hs->Add(breakdown[5]);

  TPad *pad_sector_1 = (TPad *)pad_sector->cd(1);
  pad_sector_1->cd();
  pad_sector_1->SetBottomMargin(0.15);
  pad_sector_1->SetLeftMargin(0.15);
  pad_sector_1->SetRightMargin(0.15);
  mc_hists[0]->SetTitle("Prediction");
  mc_hists[0]->Draw("COLZ");

  TPad *pad_sector_2 = (TPad *)pad_sector->cd(2);
  pad_sector_2->cd();
  pad_sector_2->SetBottomMargin(0.15);
  pad_sector_2->SetLeftMargin(0.15);
  pad_sector_2->SetRightMargin(0.15);

  data->SetTitle("Data");
  data->GetZaxis()->SetRangeUser(mc_hists[0]->GetMinimum(),mc_hists[0]->GetMaximum());
  data->Draw("COLZ same");

  std::string output_name = output_file_name + "_dxsec_d" + x_observable + "_vs_" + y_observable +"_data_vs_mc";
  if (store_root)
  {
    TFile root_file((output_location + "/TotalXSec/" + output_name + ".root").c_str(), "recreate");
    if (data)
      data->Write();
    for (unsigned int i = 0; i < mc_hists.size(); ++i)
      mc_hists[i]->Write();
    for (unsigned int i = 0; i < breakdown.size(); ++i)
      breakdown[i]->Write();
    c_sector->Write();
  }
  // Return to the canvas for further customization if needed
  c_sector->cd();
  std::filesystem::path totalxsec_path{(output_location + "/TotalXSec/").c_str()};
  if (!std::filesystem::exists(totalxsec_path))
    std::filesystem::create_directory(totalxsec_path);

  c_sector->SaveAs((output_location + "/TotalXSec/" + output_name + ".pdf").c_str());
  delete c_sector;


  std::string output_name_2 = output_file_name + "_dxsec_d" + x_observable + "_vs_" + y_observable +"_projectionX";
  plotting::PlotSlicesGeneralized(mc_hists, breakdown, data, x_observable, y_observable, y_cuts, title, output_location, output_name_2, store_root, log_scale );
}

void plotting::PlotEventRate(TH1D *data, std::string observable, std::string title, std::string data_name, std::string input_data_location,
                             std::string output_location, std::string output_file_name, std::string analysis_id, bool store_root, bool log_scale)
{
  if( log_scale ) gPad->SetLogy();
  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 700, 500);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Format plots
  if (data)
  {
    StandardFormat(data, title, kBlack, 8, observable, log_scale, 0, "Counts/Bin Width");
    data->SetLineColor(kBlue);
    data->SetLineStyle(1);
  }

  if (data)
  {
    data->Draw(" err ");
  }

  double data_integral;
  if (data)
    data_integral = data->Integral("width");
  if (data)
    std::cout << " Total integrated event rate (data) " << data_integral << std::endl;

  std::string output_name = output_file_name + "_Nevents_" + observable;

  std::filesystem::path totalxsec_path{(output_location + "/EventRate/").c_str()};
  if (!std::filesystem::exists(totalxsec_path))
    std::filesystem::create_directory(totalxsec_path);

  if (store_root)
  {
    TFile root_file((output_location + "/EventRate/" + output_name + ".root").c_str(), "recreate");
    if (data)
      data->Write();
    c1->Write();
  }
  c1->SaveAs((output_location + "/EventRate/" + output_name + ".pdf").c_str());
  delete c1;
}

void plotting::PlotEventRatePerSector(std::vector<TH1D *> data_per_sector, std::string observable, std::string title, std::string data_name, std::string input_data_location,
                                      std::string output_location, std::string output_file_name, std::string analysis_id, bool store_root, bool log_scale)
{
  if( log_scale ) gPad->SetLogy();
  TGaxis::SetMaxDigits(3);
  // Format plots
  if (data_per_sector.size() != 0)
  {
    if (data_per_sector[0])
      StandardFormat(data_per_sector[0], title + " Sector  0", kOrange + 1, 8, observable, log_scale, 0, "Counts/Bin Width");
    if (data_per_sector[1])
      StandardFormat(data_per_sector[1], title + " Sector  1", kPink + 4, 8, observable, log_scale, 0, "Counts/Bin Width");
    if (data_per_sector[2])
      StandardFormat(data_per_sector[2], title + " Sector  2", kViolet + 5, 8, observable, log_scale, 0, "Counts/Bin Width");
    if (data_per_sector[3])
      StandardFormat(data_per_sector[3], title + " Sector  3", kAzure - 5, 8, observable, log_scale, 0, "Counts/Bin Width");
    if (data_per_sector[4])
      StandardFormat(data_per_sector[4], title + " Sector  4", kTeal - 7, 8, observable, log_scale, 0, "Counts/Bin Width");
    if (data_per_sector[5])
      StandardFormat(data_per_sector[5], title + " Sector  5", kGreen - 3, 8, observable, log_scale, 0, "Counts/Bin Width");
  }
  // Draw total xsec per sectors
  TCanvas *c_sector = new TCanvas("c_sector", "c_sector", 200, 10, 700, 500);
  c_sector->cd();
  TPad *pad_sector = new TPad("pad1", "", 0, 0, 1, 1);
  pad_sector->Draw();
  pad_sector->cd();
  pad_sector->SetBottomMargin(0.15);
  pad_sector->SetLeftMargin(0.15);
  pad_sector->Divide(3, 2);

  for (unsigned int i = 0; i < 6; ++i)
  {
    TPad *pad_sector_i = (TPad *)pad_sector->cd(i + 1);
    pad_sector_i->cd();
    pad_sector_i->SetBottomMargin(0.15);
    pad_sector_i->SetLeftMargin(0.15);

    if (data_per_sector.size() != 0 && data_per_sector[i])
    {
      data_per_sector[i]->SetMarkerSize(0.7);
      data_per_sector[i]->GetYaxis()->SetTitleOffset(1.2);
      if (i == 0)
        data_per_sector[i]->Draw(" err ");
      else
        data_per_sector[i]->Draw(" err same");
    }
  }
  TGaxis::SetMaxDigits(3);
  std::string output_name = output_file_name + "_event_rate_" + observable + "_persector";
  std::filesystem::path xsecpersector_path{(output_location + "/XSecPerSector/").c_str()};
  if (!std::filesystem::exists(xsecpersector_path))
    std::filesystem::create_directory(xsecpersector_path);
  if (store_root)
    c_sector->SaveAs((output_location + "/XSecPerSector/" + output_name + ".root").c_str());
  c_sector->SaveAs((output_location + "/XSecPerSector/" + output_name + ".pdf").c_str());

  for (unsigned int i = 0; i < 6; ++i)
  {
    TPad *pad_sector_i = (TPad *)pad_sector->cd(i + 1);
    pad_sector_i->cd();
    pad_sector_i->SetBottomMargin(0.15);
    pad_sector_i->SetLeftMargin(0.15);

    if (data_per_sector.size() != 0 && data_per_sector[i])
    {
      data_per_sector[i]->SetMarkerSize(0.7);
      data_per_sector[i]->Draw(" err ");
    }
  }

  output_name = output_file_name + "_dataonly_dxsec_d" + observable + "_persector";
  xsecpersector_path = (output_location + "/XSecPerSector/").c_str();
  if (!std::filesystem::exists(xsecpersector_path))
    std::filesystem::create_directory(xsecpersector_path);

  if (store_root)
  {
    TFile root_file((output_location + "/XSecPerSector/" + output_name + ".root").c_str(), "recreate");
    for (unsigned int i = 0; i < 6; ++i)
    {
      if (data_per_sector.size() != 0 && data_per_sector[i])
      {
        data_per_sector[i]->Write();
      }
    }
    c_sector->Write();
  }

  c_sector->SaveAs((output_location + "/XSecPerSector/" + output_name + ".pdf").c_str());
  delete c_sector;
}

void plotting::PlotLegend(std::vector<TH1D *> mc_hists, std::vector<TH1D *> breakdown, TH1D *data, std::string observable,
                          std::string data_name, std::vector<std::string> model, std::string output_location,
                          std::string output_file_name, bool store_root, bool log_scale)
{
  if( log_scale ) gPad->SetLogy();
  // Store legend in separate file
  TCanvas *c_leg = new TCanvas("c_leg", "c_leg");
  c_leg->cd();
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  double LegXmin = 0.1, LegYmin = 0.65, YSpread = 0.25;
  TLegend *leg = new TLegend(LegXmin, LegYmin, LegXmin + 0.9, LegYmin + YSpread);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.08);
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  std::string model_def = "MC";
  if (model.size() > 0)
    model_def = model[0];
  if (mc_hists[0])
  {
    leg->AddEntry(mc_hists[0], ("GENIE " + model_def).c_str(), "l");
    leg->AddEntry(breakdown[0], (model_def + " EMQEL").c_str(), "l");
    leg->AddEntry(breakdown[1], (model_def + " EMRES P33(1232)").c_str(), "l");
    leg->AddEntry(breakdown[2], (model_def + " EMRES Others").c_str(), "l");
    leg->AddEntry(breakdown[3], (model_def + " EMSIS").c_str(), "l");
    leg->AddEntry(breakdown[4], (model_def + " EMMEC").c_str(), "l");
    leg->AddEntry(breakdown[5], (model_def + " EMDIS").c_str(), "l");
  }

  if (mc_hists.size() > 1)
  {
    for (unsigned int id = 1; id < mc_hists.size(); ++id)
    {
      if (!mc_hists[id])
        continue;
      std::string model_id = "Model " + std::to_string(id);
      if (model.size() == mc_hists.size())
        model_id = model[id];
      leg->AddEntry(mc_hists[id], ("GENIE " + model_id).c_str(), "l");
    }
  }

  std::string data_def = "CLAS6 data";
  if (data_name != "")
    data_def = data_name;
  if (data)
    leg->AddEntry(data, data_def.c_str(), "lp");
  leg->Draw();
  std::string output_name = output_file_name;
  if (store_root)
    c_leg->SaveAs((output_location + "/" + output_name + "_legend.root").c_str());
  c_leg->SaveAs((output_location + "/" + output_name + "_legend.pdf").c_str());

  delete c_leg;
}

void plotting::PlotPerSector(std::vector<TH1D *> mc_per_sector, std::vector<TH1D *> data_per_sector, std::string observable,
                             std::string title, std::string data_name, std::vector<std::string> model,
                             std::string input_MC_location, std::string input_data_location, std::string output_location,
                             std::string output_file_name, std::map<string, double> systematic_map,
                             std::string analysis_id, bool store_root, bool log_scale)
{
  if( log_scale ) gPad->SetLogy();
  // Format plots
  if (data_per_sector.size() != 0)
  {
    if (data_per_sector[0])
      StandardFormat(data_per_sector[0], title + " Sector  0", kOrange + 1, 8, observable, log_scale );
    if (data_per_sector[1])
      StandardFormat(data_per_sector[1], title + " Sector  1", kPink + 4, 8, observable, log_scale );
    if (data_per_sector[2])
      StandardFormat(data_per_sector[2], title + " Sector  2", kViolet + 5, 8, observable, log_scale );
    if (data_per_sector[3])
      StandardFormat(data_per_sector[3], title + " Sector  3", kAzure - 5, 8, observable, log_scale );
    if (data_per_sector[4])
      StandardFormat(data_per_sector[4], title + " Sector  4", kTeal - 7, 8, observable, log_scale );
    if (data_per_sector[5])
      StandardFormat(data_per_sector[5], title + " Sector  5", kGreen - 3, 8, observable, log_scale );
  }

  if (mc_per_sector.size() != 0)
  {
    StandardFormat(mc_per_sector[0], title + " Sector  0", kBlack, 1, observable, log_scale );
    StandardFormat(mc_per_sector[1], title + " Sector  1", kBlack, 1, observable, log_scale );
    StandardFormat(mc_per_sector[2], title + " Sector  2", kBlack, 1, observable, log_scale );
    StandardFormat(mc_per_sector[3], title + " Sector  3", kBlack, 1, observable, log_scale );
    StandardFormat(mc_per_sector[4], title + " Sector  4", kBlack, 1, observable, log_scale );
    StandardFormat(mc_per_sector[5], title + " Sector  5", kBlack, 1, observable, log_scale );
  }

  // Draw total xsec per sectors
  TCanvas *c_sector = new TCanvas("c_sector", "c_sector", 200, 10, 700, 500);
  c_sector->cd();
  TPad *pad_sector = new TPad("pad1", "", 0, 0, 1, 1);
  pad_sector->Draw();
  pad_sector->cd();
  pad_sector->SetBottomMargin(0.15);
  pad_sector->SetLeftMargin(0.15);
  pad_sector->Divide(3, 2);

  for (unsigned int i = 0; i < 6; ++i)
  {
    TPad *pad_sector_i = (TPad *)pad_sector->cd(i + 1);
    pad_sector_i->cd();
    pad_sector_i->SetBottomMargin(0.15);
    pad_sector_i->SetLeftMargin(0.15);
    mc_per_sector[i]->GetYaxis()->SetTitleOffset(1.2);
    mc_per_sector[i]->Draw("hist");

    if (data_per_sector.size() != 0 && data_per_sector[i])
    {
      data_per_sector[i]->SetMarkerSize(0.7);
      data_per_sector[i]->Draw(" err same ");
    }
  }

  std::string output_name = output_file_name + "_dxsec_d" + observable + "_persector";
  std::filesystem::path xsecpersector_path{(output_location + "/XSecPerSector/").c_str()};
  if (!std::filesystem::exists(xsecpersector_path))
    std::filesystem::create_directory(xsecpersector_path);
  if (store_root)
    c_sector->SaveAs((output_location + "/XSecPerSector/" + output_name + ".root").c_str());
  c_sector->SaveAs((output_location + "/XSecPerSector/" + output_name + ".pdf").c_str());

  for (unsigned int i = 0; i < 6; ++i)
  {
    TPad *pad_sector_i = (TPad *)pad_sector->cd(i + 1);
    pad_sector_i->cd();
    pad_sector_i->SetBottomMargin(0.15);
    pad_sector_i->SetLeftMargin(0.15);

    if (data_per_sector.size() != 0 && data_per_sector[i])
    {
      data_per_sector[i]->SetMarkerSize(0.7);
      data_per_sector[i]->Draw(" err ");
    }
  }

  output_name = output_file_name + "_dataonly_dxsec_d" + observable + "_persector";
  xsecpersector_path = (output_location + "/XSecPerSector/").c_str();
  if (!std::filesystem::exists(xsecpersector_path))
    std::filesystem::create_directory(xsecpersector_path);

  if (store_root)
  {
    TFile root_file((output_location + "/XSecPerSector/" + output_name + ".root").c_str(), "recreate");
    for (unsigned int i = 0; i < 6; ++i)
    {
      mc_per_sector[i]->Write();
      if (data_per_sector.size() != 0 && data_per_sector[i])
      {
        data_per_sector[i]->Write();
      }
    }
    c_sector->Write();
  }

  c_sector->SaveAs((output_location + "/XSecPerSector/" + output_name + ".pdf").c_str());
  delete c_sector;
}

void plotting::Plot2DXSec(std::vector<std::string> MC_files_name, std::string data_file_name,
                          std::string acceptance_file_name, std::string radcorr_file, std::string x_observable, std::string y_observable,
                          std::vector<double> & y_cuts, std::string title, std::string data_name, std::vector<std::string> model,
                          std::string input_MC_location, std::string input_data_location, std::string output_location,
                          std::string output_file_name, bool plot_data, std::map<string, double> systematic_map, std::map<std::string,std::vector<double>> cuts,
                          std::string analysis_id, bool store_root, bool log_scale)
{
  TCanvas *c = new TCanvas("c2", "c2", 800, 600);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  std::vector<TFile *> files_true_MC;
  for (unsigned int id = 0; id < MC_files_name.size(); ++id)
  {
    files_true_MC.push_back(new TFile((input_MC_location + MC_files_name[id] + "_true.root").c_str(), "ROOT"));
    if (!files_true_MC[id])
    {
      std::cout << "ERROR: the " << input_MC_location << MC_files_name[id] << "_true.root does not exist." << std::endl;
      return;
    }
  }
  TFile *file_data = nullptr;
  if (plot_data)
  {
    file_data = new TFile((input_data_location + data_file_name + ".root").c_str(), "READ");
    if (!file_data)
    {
      std::cout << "ERROR: file data doesn't exist" << std::endl;
      return;
    }
  }
  else
  {
    std::cout << "Not plotting data" << std::endl;
  }

  TFile *file_acceptance = new TFile((output_location + acceptance_file_name + ".root").c_str(), "READ");
  TFile *file_radcorr = nullptr;
  if (radcorr_file != "")
    file_radcorr = new TFile((output_location + radcorr_file + ".root").c_str(), "READ");
  if (!file_data && plot_data)
  {
    std::cout << "ERROR: the " << input_data_location << data_file_name << ".root does not exist." << std::endl;
    return;
  }
  if (!file_acceptance)
  {
    std::cout << "ERROR: the " << output_location << acceptance_file_name << ".root does not exist." << std::endl;
    return;
  }

  TH2D *h_acceptance = (TH2D *)file_acceptance->Get("Acceptance");
  TH2D *h_acceptance_0 = (TH2D *)file_acceptance->Get("Acceptance_0");
  TH2D *h_acceptance_1 = (TH2D *)file_acceptance->Get("Acceptance_1");
  TH2D *h_acceptance_2 = (TH2D *)file_acceptance->Get("Acceptance_2");
  TH2D *h_acceptance_3 = (TH2D *)file_acceptance->Get("Acceptance_3");
  TH2D *h_acceptance_4 = (TH2D *)file_acceptance->Get("Acceptance_4");
  TH2D *h_acceptance_5 = (TH2D *)file_acceptance->Get("Acceptance_5");

  TH2D *h_radcorr = nullptr;
  if (file_radcorr)
  {
    h_radcorr = (TH2D *)file_radcorr->Get("Acceptance");
  }

  // Get Tree for main model
  TTree *tree_true = (TTree *)files_true_MC[0]->Get("MCCLAS6Tree");
  // Get configured energy, used for plotting
  double BeamE;
  tree_true->SetBranchAddress("BeamE", &BeamE);
  tree_true->GetEntry(0);

  // For submodels only total prediction is plotted
  std::vector<TTree *> tree_submodels;
  for (unsigned int id = 1; id < MC_files_name.size(); ++id)
  {
    tree_submodels.push_back((TTree *)files_true_MC[id]->Get("MCCLAS6Tree"));
    if (!tree_submodels[id - 1])
    {
      std::cout << "ERROR: the threes do not exist." << std::endl;
      return;
    }
  }

  TTree *tree_data = nullptr;
  if (plot_data)
  {
    tree_data = (TTree *)file_data->Get("CLAS6Tree");
    if (!tree_data)
    {
      std::cout << "ERROR: tree data doesn't exist" << std::endl;
      return;
    }
  }

  if (!h_acceptance )
  {
    std::cout << "ERROR: Acceptance is not defined" << std::endl;
    return;
  }
  if (!tree_true || (!tree_data && plot_data))
  {
    std::cout << "ERROR: the threes do not exist." << std::endl;
    return;
  }

  // Create histogram for total and total xsec per sector
  TH2D *hist_true = (TH2D *)h_acceptance->Clone();
  hist_true->SetName("MC_True");
  hist_true->Reset();

  TH2D *hist_true_0 = (TH2D *)h_acceptance_0->Clone();
  hist_true_0->SetName("MC_True_Sector_0");
  hist_true_0->Reset();
  TH2D *hist_true_1 = (TH2D *)h_acceptance_1->Clone();
  hist_true_1->SetName("MC_True_Sector_1");
  hist_true_1->Reset();
  TH2D *hist_true_2 = (TH2D *)h_acceptance_2->Clone();
  hist_true_2->SetName("MC_True_Sector_2");
  hist_true_2->Reset();
  TH2D *hist_true_3 = (TH2D *)h_acceptance_3->Clone();
  hist_true_3->SetName("MC_True_Sector_3");
  hist_true_3->Reset();
  TH2D *hist_true_4 = (TH2D *)h_acceptance_4->Clone();
  hist_true_4->SetName("MC_True_Sector_4");
  hist_true_4->Reset();
  TH2D *hist_true_5 = (TH2D *)h_acceptance_5->Clone();
  hist_true_5->SetName("MC_True_Sector_5");
  hist_true_5->Reset();

  // Breakdown histograms for total (all sectors only):
  TH2D *hist_true_QEL = (TH2D *)h_acceptance->Clone();
  hist_true_QEL->SetName("MC_True_QEL");
  hist_true_QEL->Reset();
  TH2D *hist_true_RES_Delta = (TH2D *)h_acceptance->Clone();
  hist_true_RES_Delta->SetName("MC_True_RES_Delta");
  hist_true_RES_Delta->Reset();
  TH2D *hist_true_RES = (TH2D *)h_acceptance->Clone();
  hist_true_RES->SetName("MC_True_RES");
  hist_true_RES->Reset();
  TH2D *hist_true_SIS = (TH2D *)h_acceptance->Clone();
  hist_true_SIS->SetName("MC_True_SIS");
  hist_true_SIS->Reset();
  TH2D *hist_true_MEC = (TH2D *)h_acceptance->Clone();
  hist_true_MEC->SetName("MC_True_MEC");
  hist_true_MEC->Reset();
  TH2D *hist_true_DIS = (TH2D *)h_acceptance->Clone();
  hist_true_DIS->SetName("MC_True_DIS");
  hist_true_DIS->Reset();

  // Same per model - only total prediction
  std::vector<TH2D *> hists_true_submodel;
  for (unsigned int id = 1; id < MC_files_name.size(); ++id)
  {
    hists_true_submodel.push_back((TH2D *)h_acceptance->Clone());
    hists_true_submodel[id - 1]->SetName(("MC_True_Model_" + std::to_string(id)).c_str());
    hists_true_submodel[id - 1]->Reset();
    hists_true_submodel[id - 1]->SetLineWidth(3);
  }

  // total and per sector
  TH2D *hist_data = nullptr, *hist_data_0 = nullptr, *hist_data_1 = nullptr, *hist_data_2 = nullptr, *hist_data_3 = nullptr, *hist_data_4 = nullptr, *hist_data_5 = nullptr;
  if (plot_data)
  {
    hist_data = (TH2D *)h_acceptance->Clone();
    hist_data->SetName("Data");
    hist_data->Reset();

    hist_data_0 = (TH2D *)h_acceptance_0->Clone();
    hist_data_0->SetName("Data_Sector_0");
    hist_data_0->Reset();
    hist_data_1 = (TH2D *)h_acceptance_1->Clone();
    hist_data_1->SetName("Data_Sector_1");
    hist_data_1->Reset();
    hist_data_2 = (TH2D *)h_acceptance_2->Clone();
    hist_data_2->Reset();
    hist_data_2->SetName("Data_Sector_2");
    hist_data_3 = (TH2D *)h_acceptance_3->Clone();
    hist_data_3->SetName("Data_Sector_3");
    hist_data_3->Reset();
    hist_data_4 = (TH2D *)h_acceptance_4->Clone();
    hist_data_4->SetName("Data_Sector_4");
    hist_data_4->Reset();
    hist_data_5 = (TH2D *)h_acceptance_5->Clone();
    hist_data_5->SetName("Data_Sector_5");
    hist_data_5->Reset();
  }

  std::vector<TTree *> trees = {tree_true};
  if (plot_data)
    trees.push_back(tree_data);

  std::vector<TH2D *> hists = {hist_true, hist_data, hist_true_0, hist_data_0,
                               hist_true_1, hist_data_1, hist_true_2, hist_data_2,
                               hist_true_3, hist_data_3, hist_true_4, hist_data_4,
                               hist_true_5, hist_data_5};

  unsigned int size_primary_trees = trees.size();
  unsigned int size_primary_hists = hists.size();
  // Adding total predictions for alternative models
  for (unsigned int id = 1; id < MC_files_name.size(); ++id)
  {
    trees.push_back(tree_submodels[id - 1]);
    hists.push_back(hists_true_submodel[id - 1]);
  }

  // If data is plot, the position of its trees and histograms is 1.
  // Otherwise we set it to a big number so it is ignored
  unsigned int id_data = 9999;
  if (plot_data)
    id_data = 1;

  // OBSERVABLE DEFINITION specific for MC
  long NEntries;
  bool IsBkg;
  int ElectronSector;
  bool QEL, RES, DIS, MEC;
  double MCNormalization, DataNormalization;
  std::vector<double> mc_norm;
  int resid;

  for (unsigned int i = 0; i < trees.size(); ++i)
  {
    NEntries = trees[i]->GetEntries();
    if (!trees[i])
      continue;
    trees[i]->SetBranchAddress("TotWeight", &TotWeight);
    trees[i]->SetBranchAddress("IsBkg", &IsBkg);
    trees[i]->SetBranchAddress("ECal", &ECal);
    trees[i]->SetBranchAddress("pfl_theta", &pfl_theta);
    trees[i]->SetBranchAddress("pfl_phi", &pfl_phi);
    trees[i]->SetBranchAddress("pfl", &pfl);
    trees[i]->SetBranchAddress("proton_mom", &proton_mom);
    trees[i]->SetBranchAddress("proton_theta", &proton_theta);
    trees[i]->SetBranchAddress("proton_phi", &proton_phi);
    trees[i]->SetBranchAddress("pim_mom", &pim_mom);
    trees[i]->SetBranchAddress("pim_theta", &pim_theta);
    trees[i]->SetBranchAddress("pim_phi", &pim_phi);
    trees[i]->SetBranchAddress("pip_mom", &pip_mom);
    trees[i]->SetBranchAddress("pip_theta", &pip_theta);
    trees[i]->SetBranchAddress("pip_phi", &pip_phi);
    trees[i]->SetBranchAddress("RecoW", &RecoW);
    trees[i]->SetBranchAddress("Recoq3", &Recoq3);
    trees[i]->SetBranchAddress("RecoQELEnu", &RecoQELEnu);
    trees[i]->SetBranchAddress("RecoXBJK", &RecoXBJK);
    trees[i]->SetBranchAddress("RecoQ2", &RecoQ2);
    trees[i]->SetBranchAddress("RecoEnergyTransfer", &RecoEnergyTransfer);
    trees[i]->SetBranchAddress("AlphaT", &AlphaT);
    trees[i]->SetBranchAddress("HadAlphaT", &HadAlphaT);
    trees[i]->SetBranchAddress("DeltaPT", &DeltaPT);
    trees[i]->SetBranchAddress("HadDeltaPT", &HadDeltaPT);
    trees[i]->SetBranchAddress("HadDeltaPTx", &HadDeltaPTx);
    trees[i]->SetBranchAddress("HadDeltaPTy", &HadDeltaPTy);
    trees[i]->SetBranchAddress("DeltaPhiT", &DeltaPhiT);
    trees[i]->SetBranchAddress("HadDeltaPhiT", &HadDeltaPhiT);
    trees[i]->SetBranchAddress("ElectronSector", &ElectronSector);
    trees[i]->SetBranchAddress("HadSystemMass", &HadSystemMass);
    trees[i]->SetBranchAddress("MissingEnergy", &MissingEnergy);
    trees[i]->SetBranchAddress("MissingAngle", &MissingAngle);
    trees[i]->SetBranchAddress("MissingMomentum", &MissingMomentum);
    trees[i]->SetBranchAddress("InferedNucleonMom", &InferedNucleonMom);
    trees[i]->SetBranchAddress("HadronsAngle", &HadronsAngle);
    trees[i]->SetBranchAddress("AdlerAngleThetaP", &AdlerAngleThetaP);
    trees[i]->SetBranchAddress("AdlerAnglePhiP", &AdlerAnglePhiP);
    trees[i]->SetBranchAddress("AdlerAngleThetaPi", &AdlerAngleThetaPi);
    trees[i]->SetBranchAddress("AdlerAnglePhiPi", &AdlerAnglePhiPi);
    trees[i]->SetBranchAddress("Angleqvshad", &Angleqvshad);
    trees[i]->SetBranchAddress("RecoEvPion", &RecoEvPion);
    trees[i]->SetBranchAddress("RecoWPion", &RecoWPion);
    trees[i]->SetBranchAddress("ElectronPT", &ElectronPT);
    trees[i]->SetBranchAddress("PionPT", &PionPT);

    // Only fill true info for the first model:
    if (i == 0)
    {
      trees[i]->SetBranchAddress("QEL", &QEL);
      trees[i]->SetBranchAddress("RES", &RES);
      trees[i]->SetBranchAddress("MEC", &MEC);
      trees[i]->SetBranchAddress("DIS", &DIS);
      trees[i]->SetBranchAddress("resid", &resid);
    }
    // Only second tree corresponds to data
    if (i != id_data)
    {
      trees[i]->SetBranchAddress("MCNormalization", &MCNormalization);
    }
    else
    {
      trees[i]->SetBranchAddress("DataNormalization", &DataNormalization);
    }

    for (int j = 0; j < NEntries; ++j)
    {
      trees[i]->GetEntry(j);
      double content_x = GetObservable(x_observable);
      double content_y = GetObservable(y_observable);
      double w = TotWeight;
      if (i != id_data && j == 0)
        mc_norm.push_back(MCNormalization);

      // Check if passes cuts
      bool do_fill =true ;
      for (auto it = cuts.begin(); it != cuts.end(); it++)
      {
        double min = it->second[0] ;
        double max = it->second[1] ;
        if( GetObservable(it->first) < min || GetObservable(it->first) > max ) {
          do_fill = false;
          continue;
        }
      }
      if( !do_fill ) continue ;

      unsigned int id_hist = i;

      // Fill the per Sector  histogram. Only for primary model
      if (hists[size_primary_trees * (ElectronSector + 1) + i])
      {
        if (i < size_primary_trees)
        {
          hists[size_primary_trees * (ElectronSector + 1) + i]->Fill(content_x, content_y, w);
        }
      }
      if (i > size_primary_trees - 1)
        id_hist = size_primary_hists + (i - size_primary_trees);

      if (hists[id_hist])
      {
        hists[id_hist]->Fill(content_x, content_y, w);
        hists[id_hist]->SetLineWidth(3);
      }

      if (i == 0)
      {
        if (QEL)
          hist_true_QEL->Fill(content_x, content_y, w);
        if (RES)
        {
          if (resid == 0)
            hist_true_RES_Delta->Fill(content_x, content_y, w);
          else
            hist_true_RES->Fill(content_x, content_y, w);
        }
        if (DIS)
        {
          if (RecoW < 1.7)
            hist_true_SIS->Fill(content_x, content_y, w);
          else
            hist_true_DIS->Fill(content_x, content_y, w);
        }
        if (MEC)
          hist_true_MEC->Fill(content_x, content_y, w);
      }
    }
  }

  // Store uncorrected data
  TH2D *hist_data_uncorr = nullptr, *hist_data_uncorr_0 = nullptr, *hist_data_uncorr_1 = nullptr, *hist_data_uncorr_2 = nullptr, *hist_data_uncorr_3 = nullptr, *hist_data_uncorr_4 = nullptr, *hist_data_uncorr_5 = nullptr;
  // Store corrected for acceptance but not for radiation
  TH2D *hist_data_uncorrrad = nullptr;
  // Corr event rate
  TH2D *hist_data_correventrate = nullptr, *hist_data_correventrate_0 = nullptr, *hist_data_correventrate_1 = nullptr, *hist_data_correventrate_2 = nullptr, *hist_data_correventrate_3 = nullptr, *hist_data_correventrate_4 = nullptr, *hist_data_correventrate_5 = nullptr;
  // Event rate with Systematics
  TH2D *hist_data_correventrate_wsyst = nullptr, *hist_data_correventrate_wsyst_0 = nullptr, *hist_data_correventrate_wsyst_1 = nullptr, *hist_data_correventrate_wsyst_2 = nullptr, *hist_data_correventrate_wsyst_3 = nullptr, *hist_data_correventrate_wsyst_4 = nullptr, *hist_data_correventrate_wsyst_5 = nullptr;
  // Event rate with Systematics
  TH2D *hist_xsec_wsyst = nullptr, *hist_xsec_wsyst_0 = nullptr, *hist_xsec_wsyst_1 = nullptr, *hist_xsec_wsyst_2 = nullptr, *hist_xsec_wsyst_3 = nullptr, *hist_xsec_wsyst_4 = nullptr, *hist_xsec_wsyst_5 = nullptr;

  // Store data event rate before acceptance correction:
  if (plot_data && hist_data)
  {
    hist_data_uncorr = (TH2D *)hist_data->Clone();
    hist_data_uncorr->SetName("Uncorrected Data");
    hist_data_uncorr_0 = (TH2D *)hist_data_0->Clone();
    hist_data_uncorr_0->SetName("Uncorrected Data Sector  0");
    hist_data_uncorr_1 = (TH2D *)hist_data_1->Clone();
    hist_data_uncorr_1->SetName("Uncorrected Data Sector  1");
    hist_data_uncorr_2 = (TH2D *)hist_data_2->Clone();
    hist_data_uncorr_2->SetName("Uncorrected Data Sector  2");
    hist_data_uncorr_3 = (TH2D *)hist_data_3->Clone();
    hist_data_uncorr_3->SetName("Uncorrected Data Sector  3");
    hist_data_uncorr_4 = (TH2D *)hist_data_4->Clone();
    hist_data_uncorr_4->SetName("Uncorrected Data Sector  4");
    hist_data_uncorr_5 = (TH2D *)hist_data_5->Clone();
    hist_data_uncorr_5->SetName("Uncorrected Data Sector  5");

    // Normaize for bin size
    NormalizeHist(hist_data_uncorr, 1);
    NormalizeHist(hist_data_uncorr_0, 1);
    NormalizeHist(hist_data_uncorr_1, 1);
    NormalizeHist(hist_data_uncorr_2, 1);
    NormalizeHist(hist_data_uncorr_3, 1);
    NormalizeHist(hist_data_uncorr_4, 1);
    NormalizeHist(hist_data_uncorr_5, 1);

    // Correct data for detector acceptance :
    // Notice this step already propagates the error from the acceptance to the corrected event rate
    // (Err_Corr_eventrate)^2 = (Err_Raw_EventRate)^2 * (Acc)^2 + (Raw_EventRate)^2 * (Err_Acc)^2
    std::cout << " Acceptance correcting data " << std::endl;
    CorrectData(hist_data, h_acceptance);
    CorrectData(hist_data_0, h_acceptance_0);
    CorrectData(hist_data_1, h_acceptance_1);
    CorrectData(hist_data_2, h_acceptance_2);
    CorrectData(hist_data_3, h_acceptance_3);
    CorrectData(hist_data_4, h_acceptance_4);
    CorrectData(hist_data_5, h_acceptance_5);

    hist_data_uncorrrad = (TH2D *)hist_data->Clone();
    hist_data_uncorrrad->SetName("Corrected for acceptance before rad corr data");
    NormalizeHist(hist_data_uncorrrad, 1);

    // Apply radiative correction
    if (h_radcorr)
    {
      std::cout << " Radiative correcting data " << std::endl;
      CorrectData(hist_data, h_radcorr);
      CorrectData(hist_data_0, h_radcorr);
      CorrectData(hist_data_1, h_radcorr);
      CorrectData(hist_data_2, h_radcorr);
      CorrectData(hist_data_3, h_radcorr);
      CorrectData(hist_data_4, h_radcorr);
      CorrectData(hist_data_5, h_radcorr);
    }

    hist_data_correventrate = (TH2D *)hist_data->Clone();
    hist_data_correventrate->SetName("Corrected_Event_Rate_Data");
    hist_data_correventrate_0 = (TH2D *)hist_data_0->Clone();
    hist_data_correventrate_0->SetName("Corrected_Event_Rate_Data_Sector_0");
    hist_data_correventrate_1 = (TH2D *)hist_data_1->Clone();
    hist_data_correventrate_1->SetName("Corrected_Event_Rate_Data_Sector_1");
    hist_data_correventrate_2 = (TH2D *)hist_data_2->Clone();
    hist_data_correventrate_2->SetName("Corrected_Event_Rate_Data_Sector_2");
    hist_data_correventrate_3 = (TH2D *)hist_data_3->Clone();
    hist_data_correventrate_3->SetName("Corrected_Event_Rate_Data_Sector_3");
    hist_data_correventrate_4 = (TH2D *)hist_data_4->Clone();
    hist_data_correventrate_4->SetName("Corrected_Event_Rate_Data_Sector_4");
    hist_data_correventrate_5 = (TH2D *)hist_data_5->Clone();
    hist_data_correventrate_5->SetName("Corrected_Event_Rate_Data_Sector_5");

    // Store event rate with systematics
    hist_data_correventrate_wsyst = (TH2D *)hist_data->Clone();
    hist_data_correventrate_wsyst->SetName("Corrected Event Rate with Systematics Data");
    hist_data_correventrate_wsyst_0 = (TH2D *)hist_data_0->Clone();
    hist_data_correventrate_wsyst_0->SetName("Corrected Event Rate with Systematics Sector 0");
    hist_data_correventrate_wsyst_1 = (TH2D *)hist_data_1->Clone();
    hist_data_correventrate_wsyst_1->SetName("Corrected Event Rate with Systematics Sector  1");
    hist_data_correventrate_wsyst_2 = (TH2D *)hist_data_2->Clone();
    hist_data_correventrate_wsyst_2->SetName("Corrected Event Rate with Systematics Sector  2");
    hist_data_correventrate_wsyst_3 = (TH2D *)hist_data_3->Clone();
    hist_data_correventrate_wsyst_3->SetName("Corrected Event Rate with Systematics Sector  3");
    hist_data_correventrate_wsyst_4 = (TH2D *)hist_data_4->Clone();
    hist_data_correventrate_wsyst_4->SetName("Corrected Event Rate with Systematics Sector  4");
    hist_data_correventrate_wsyst_5 = (TH2D *)hist_data_5->Clone();
    hist_data_correventrate_wsyst_5->SetName("Corrected Event Rate with Systematics Sector  5");
    std::cout << " Normalizing data "<<std::endl;
    // Normaize by bin width
    NormalizeHist(hist_data_correventrate_wsyst, 1);
    NormalizeHist(hist_data_correventrate_wsyst_0, 1);
    NormalizeHist(hist_data_correventrate_wsyst_1, 1);
    NormalizeHist(hist_data_correventrate_wsyst_2, 1);
    NormalizeHist(hist_data_correventrate_wsyst_3, 1);
    NormalizeHist(hist_data_correventrate_wsyst_4, 1);
    NormalizeHist(hist_data_correventrate_wsyst_5, 1);

    // Normaize by bin width
    NormalizeHist(hist_data_correventrate, 1);
    NormalizeHist(hist_data_correventrate_0, 1);
    NormalizeHist(hist_data_correventrate_1, 1);
    NormalizeHist(hist_data_correventrate_2, 1);
    NormalizeHist(hist_data_correventrate_3, 1);
    NormalizeHist(hist_data_correventrate_4, 1);
    NormalizeHist(hist_data_correventrate_5, 1);

    // Normalize to cross-section
    NormalizeHist(hist_data, DataNormalization);
    NormalizeHist(hist_data_0, DataNormalization);
    NormalizeHist(hist_data_1, DataNormalization);
    NormalizeHist(hist_data_2, DataNormalization);
    NormalizeHist(hist_data_3, DataNormalization);
    NormalizeHist(hist_data_4, DataNormalization);
    NormalizeHist(hist_data_5, DataNormalization);

    // Add Systematics
    // 1 - Acceptance model dependence (already included)
    // 2 - Sector Sector Variation
    // 3 - Relative uncertanties from configuration
    std::cout << " Adding systematics " << std::endl;
    // Adding Acceptance correction systematics from model dependence
    TH2D *hist_syst_acc = systematics::AddSystematic(*hist_data, *h_acceptance);

    TCanvas *cacc = new TCanvas("cacc", "cacc", 800, 600);
    hist_syst_acc->Draw("hist");
    cacc->SaveAs((output_location + "/XSecPerSector/" + output_file_name + "_syst_accmodel_" + x_observable + "_vs_" + y_observable + ".root").c_str());
    delete cacc;

    // Add sector variation ERROR. Store relative error in histogram
    // We use the bkg substracted, eff corrected distributions for the calculation
    TH2D *hist_syst_sector = systematics::SectorVariationError(*hist_data, {hist_data_0, hist_data_1, hist_data_2, hist_data_3, hist_data_4, hist_data_5});
    systematics::SectorVariationError(*hist_data_correventrate_wsyst, {hist_data_0, hist_data_1, hist_data_2, hist_data_3, hist_data_4, hist_data_5});

    TCanvas *csect = new TCanvas("csect", "csect", 800, 600);
    hist_syst_sector->Draw("hist");
    csect->SaveAs((output_location + "/XSecPerSector/" + output_file_name + "_syst_persector_" + x_observable + "_vs_" + y_observable + ".root").c_str());
    delete csect;

    // adding systematics from systematic map. Relative systematic added to all bins
    for (auto it = systematic_map.begin(); it != systematic_map.end(); ++it)
    {
      std::cout << " Adding " << it->second << " % systematic on " << it->first << std::endl;
      systematics::AddSystematic(*hist_data, it->second, it->first);
      systematics::AddSystematic(*hist_data_correventrate_wsyst, it->second, it->first);
    }

    // Hard coding some well known systematics
    systematics::AddSystematic(*hist_data, 5, "Radiative");
    systematics::AddSystematic(*hist_data, 3, "Normalization");
    systematics::AddSystematic(*hist_data, 1, "AnglDependence");
    systematics::AddSystematic(*hist_data, 1, "MaxMultiplicity");

    systematics::AddSystematic(*hist_data_correventrate_wsyst, 5, "Radiative");
    systematics::AddSystematic(*hist_data_correventrate_wsyst, 3, "Normalization");
    systematics::AddSystematic(*hist_data_correventrate_wsyst, 1, "AnglDependence");
    systematics::AddSystematic(*hist_data_correventrate_wsyst, 1, "MaxMultiplicity");

  } // end if data

  // Normalize MC to cross-section
  for (unsigned int id = 0; id < hists_true_submodel.size(); ++id)
  {
    NormalizeHist(hists_true_submodel[id], mc_norm[id + 1]);
    StandardFormat(hists_true_submodel[id], title, kBlack, 2 + id, x_observable, y_observable);
  }

  NormalizeHist(hist_true, mc_norm[0]);
  NormalizeHist(hist_true_0, mc_norm[0]);
  NormalizeHist(hist_true_1, mc_norm[0]);
  NormalizeHist(hist_true_2, mc_norm[0]);
  NormalizeHist(hist_true_3, mc_norm[0]);
  NormalizeHist(hist_true_4, mc_norm[0]);
  NormalizeHist(hist_true_5, mc_norm[0]);
  NormalizeHist(hist_true_QEL, mc_norm[0]);
  NormalizeHist(hist_true_RES, mc_norm[0]);
  NormalizeHist(hist_true_RES_Delta, mc_norm[0]);
  NormalizeHist(hist_true_SIS, mc_norm[0]);
  NormalizeHist(hist_true_MEC, mc_norm[0]);
  NormalizeHist(hist_true_DIS, mc_norm[0]);

  // Store histograms for plotting
  std::vector<TH2D> mc_hists = {*hist_true};
  std::vector<TH2D *> mc_hists_xsec = {hist_true};
  for (unsigned int id = 0; id < hists_true_submodel.size(); ++id)
  {
    mc_hists.push_back(*hists_true_submodel[id]);
    mc_hists_xsec.push_back(hists_true_submodel[id]);
  }

  std::vector<TH2D> breakdown = {*hist_true_QEL, *hist_true_RES_Delta, *hist_true_RES, *hist_true_SIS, *hist_true_MEC, *hist_true_DIS};
  std::vector<TH2D *> breakdown_xsec = {hist_true_QEL, hist_true_RES_Delta, hist_true_RES, hist_true_SIS, hist_true_MEC, hist_true_DIS};

  plotting::PlotTotal2DXSec(mc_hists_xsec, breakdown_xsec, hist_data, x_observable, y_observable, y_cuts, title, data_name, model, input_MC_location,
                          input_data_location, output_location, output_file_name + "_with_breakdown", systematic_map, true, analysis_id, store_root, log_scale);


}

void plotting::CreateCanvasWithPads(TPad *& pad, std::vector<TPad*>& topPad, std::vector<TPad*>& bottomPad, const std::string& canvasName) {
    // We are using the canvas previously created and sliced into N pads. For each pad we create a main plot and a ratio plot.
    pad->cd();

    topPad.push_back(new TPad("topPad", "Top Pad", 0, 0.35, 1, 1));
    topPad.back()->SetBottomMargin(0.03);
    topPad.back()->SetLeftMargin(0.2);
    topPad.back()->SetRightMargin(0.02);
    topPad.back()->Draw();

    bottomPad.push_back( new TPad("bottomPad", "Bottom Pad", 0, 0, 1, 0.35));
    bottomPad.back()->SetTopMargin(0.03);
    bottomPad.back()->SetBottomMargin(0.25);
    bottomPad.back()->SetLeftMargin(0.2);
    bottomPad.back()->SetRightMargin(0.02);
    bottomPad.back()->Draw();
}

void plotting::PlotProjectionWithRatio( const std::vector<TH2D*>& mcHists, const std::vector<TH2D*>& breakdown, const TH2D* data, const std::string& xobservable, const std::string& yobservable,
    TPad* topPad, TPad* bottomPad, bool logScale, const std::string& axis, double y_cut_min, double y_cut_max )
{
    if( y_cut_min > y_cut_max ) throw std::invalid_argument("Requested cuts are not valid.");
    topPad->cd();
    if (logScale) topPad->SetLogy();

    // Change min and max for limits if outside of range
    double max_range = mcHists[0]->GetYaxis()->GetBinLowEdge(mcHists[0]->GetYaxis()->GetNbins())+ mcHists[0]->GetYaxis()->GetBinWidth(mcHists[0]->GetYaxis()->GetNbins());;
    double min_range = mcHists[0]->GetYaxis()->GetBinLowEdge(1);
    if( y_cut_min < min_range ) y_cut_min = min_range;
    if( y_cut_max > max_range ) y_cut_max = max_range;

    // Projection and formatting
    std::vector<TH1D*> mcProjections;
    auto breakdownStack = new THStack(("projection_components_" + axis).c_str(), "");

    // Store input for projection details
    std::string projection_name = axis == "X" ? "_px_Range_"+std::to_string(y_cut_min)+"_"+std::to_string(y_cut_max) : "_py";
    double first_bin = 0 ;
    double last_bin = 999;
    const double bin_max_cut_id = GetClosestBin( mcHists[0],y_cut_max, "Y" );
    const double bin_min_cut_id = GetClosestBin( mcHists[0],y_cut_min, "Y" );

    // Compute slice for main predictions
    for ( unsigned int i = 0 ; i < mcHists.size(); ++i ) {
        double first_bin_axis = mcHists[i]->GetNbinsX();
        last_bin = axis == "X" ? mcHists[i]->GetNbinsX() : mcHists[i]->GetNbinsY();
        if( axis == "X" && last_bin > bin_max_cut_id ) last_bin = bin_max_cut_id ;
        if( axis == "X" && first_bin < bin_min_cut_id ) first_bin = bin_min_cut_id ;
        mcProjections.push_back(axis == "X" ? mcHists[i]->ProjectionX((mcHists[i]->GetName()+projection_name).c_str(),first_bin,last_bin) : mcHists[i]->ProjectionY((mcHists[i]->GetName()+projection_name).c_str(),first_bin,last_bin));
        StandardFormat(mcProjections[i], "", kBlack, i+1, xobservable, logScale);
    }

    double total_min = 0 ;
    if( logScale) total_min = 1E-4;
    // Store maximum of all TH2D histograms to set the same range to all
    double total_max = plotting::GetMaximum(mcProjections);

    // Store breakdown for default prediction (first)
    for (unsigned int i = 0 ; i < breakdown.size(); ++i ) {
        last_bin = axis == "X" ? breakdown[i]->GetNbinsX() : breakdown[i]->GetNbinsY();
        if( axis == "X" && last_bin > bin_max_cut_id ) last_bin = bin_max_cut_id ;
        if( axis == "X" && first_bin < bin_min_cut_id ) first_bin = bin_min_cut_id ;
        breakdownStack->Add(axis == "X" ? breakdown[i]->ProjectionX((breakdown[i]->GetName()+projection_name).c_str(),first_bin,last_bin) : breakdown[i]->ProjectionY((breakdown[i]->GetName()+projection_name).c_str(),first_bin,last_bin));
    }

    // Determine maximum y-axis range
    mcProjections[0]->GetYaxis()->SetTitleOffset(1.4);
    mcProjections[0]->GetXaxis()->SetLabelSize(0);
    mcProjections[0]->SetMarkerSize(0);
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << y_cut_min << "<" << plotting::GetAxisLabel(yobservable, 0) << "<" << y_cut_max;
    mcProjections[0]->SetTitle(oss.str().c_str());
    mcProjections[0]->SetTitleSize(0.1);
    mcProjections[0]->GetYaxis()->SetTitleOffset(1.2);
    mcProjections[0]->Draw("hist err");
    breakdownStack->Draw("hist err same");
    for (size_t i = 0; i < mcProjections.size(); ++i) {
        mcProjections[i]->GetYaxis()->SetRangeUser(total_min,total_max); // TO AUTOMATIZE!
        mcProjections[i]->GetXaxis()->SetTitleSize(0.05);
        mcProjections[i]->GetYaxis()->SetTitleSize(0.09);
        mcProjections[i]->GetYaxis()->SetTitleOffset(1);
        mcProjections[i]->SetMarkerSize(0);
        mcProjections[i]->SetLineWidth(1.2);
        mcProjections[i]->Draw("hist err same");
    }

    // Data projection
    last_bin = axis == "X" ? data->GetNbinsX() : data->GetNbinsY();
    if( axis == "X" && last_bin > bin_max_cut_id ) last_bin = bin_max_cut_id ;
    if( axis == "X" && first_bin < bin_min_cut_id ) first_bin = bin_min_cut_id ;
    TH1D* dataProjection = axis == "X" ? data->ProjectionX((data->GetName()+projection_name).c_str(),first_bin,last_bin) : data->ProjectionY((data->GetName()+projection_name).c_str(),first_bin,last_bin);
    StandardFormat(dataProjection, "", kBlack, 8, xobservable, logScale);
    dataProjection->SetLineStyle(1);
    dataProjection->SetMarkerSize(0.6);
    dataProjection->SetLineWidth(1);
    dataProjection->Draw("err same");

    // Ratio calculation
    bottomPad->cd();
    if( logScale ) bottomPad->SetLogy();
    TH1D* dataRatio = (TH1D*)dataProjection->Clone();
    std::vector<TH1D*> mcRatios;

    for (unsigned int i = 0 ; i < mcProjections.size(); ++i){
        auto ratio = (TH1D*)mcProjections[i]->Clone();
        ratio->Divide(dataProjection);
        StandardFormat(ratio, "", kBlack, i + 1, xobservable, logScale);
        ratio->SetMarkerSize(0);
        mcRatios.push_back(ratio);
    }

    dataRatio->Divide(dataProjection);
    StandardFormat(dataRatio, "", kBlack, 8, xobservable, logScale);
    dataRatio->GetXaxis()->SetLabelSize(0.1);
    dataRatio->GetXaxis()->SetTitleSize(0.1);
    dataRatio->GetYaxis()->SetLabelSize(0.1);
    dataRatio->GetYaxis()->SetTitleSize(0.1);
    dataRatio->SetLineStyle(1);
    dataRatio->SetMinimum(0);
    dataRatio->GetYaxis()->SetMaxDigits(5);
    dataRatio->GetYaxis()->SetTitleOffset(0.7);
    dataRatio->GetYaxis()->SetTitle("Pred / Data ");

    double maxRatio = GetMaximum(mcRatios)*(1-0.35);
    double minRatio = GetMinimum(mcRatios)*(1-0.15);
    dataRatio->GetYaxis()->SetRangeUser(minRatio, maxRatio);
    dataRatio->SetLineWidth(1);
    dataRatio->GetXaxis()->SetLabelSize(0.08);
    dataRatio->GetYaxis()->SetTitleSize(0.1);
    dataRatio->GetYaxis()->SetLabelSize(0.07);
    dataRatio->GetYaxis()->SetTitleSize(0.1);
    dataRatio->GetYaxis()->SetTitleOffset(0.77);
    dataRatio->Draw("err");
    for (const auto& ratio : mcRatios) {
        ratio->SetLineWidth(1.2);
        ratio->Draw("hist err same");
    }
}

void plotting::PlotSlicesGeneralized(const std::vector<TH2D*>& mcHists, const std::vector<TH2D*>& breakdown, const TH2D* data,
    const std::string& xObservable, const std::string& yObservable, std::vector<double> & y_cuts, const std::string& title, const std::string& outputLocation,
    const std::string& outputName, bool storeRoot, bool logScale ){

    // We want to slice the histogram on the Y axis - which is the alternative axis
    // The input is the cuts - ignoring the histogram edges.
    // The number of canvas will be y_cuts.size()+1
    int v_slice = 1, h_slice = y_cuts.size() + 1 ;
    if( y_cuts.size() > 2 ) {v_slice = 2; h_slice = 3; }
    else if( y_cuts.size() > 5 ) {v_slice = 3; h_slice = 3; }

    TCanvas *c_slices = new TCanvas("c_slices", "c_slices", 200, 10, 700, 500);
    c_slices->cd();
    TPad *pad_slice = new TPad("pad1", "", 0, 0, 1, 1);
    pad_slice->Draw();
    pad_slice->cd();
    pad_slice->SetBottomMargin(0.15);
    pad_slice->SetLeftMargin(0.15);
    pad_slice->SetRightMargin(0.15);
    pad_slice->Divide(h_slice, v_slice);

    // We add two edges for the cuts so we include, effectively, the min and max range of a histogram
    y_cuts.insert(y_cuts.begin(), -99999);
    y_cuts.push_back(99999);

    std::vector<TPad *> pad_slices ;
    std::vector<TPad *> top_pads, bottom_pads ;
    for(unsigned int i = 0 ; i < y_cuts.size() - 1; ++i ){
      pad_slices.push_back((TPad *)pad_slice->cd(i+1));
      pad_slices[i]->cd();
      CreateCanvasWithPads(pad_slices[i], top_pads, bottom_pads, "CanvasX");
      PlotProjectionWithRatio(mcHists, breakdown, data, xObservable, yObservable, top_pads[i], bottom_pads[i], logScale, "X", y_cuts[i], y_cuts[i+1]);
    }

    c_slices->SaveAs((outputLocation + "/TotalXSec/" + outputName + "_X.pdf").c_str());
    c_slices->SaveAs((outputLocation + "/TotalXSec/" + outputName + "_X.root").c_str());

    if (storeRoot) {
        TFile rootFile((outputLocation + "/TotalXSec/" + outputName + ".root").c_str(), "recreate");
        for (auto& proj : mcHists) {
            proj->Write();
        }
        data->Write();
        c_slices->Write();
        rootFile.Close();
    }
}

int plotting::GetClosestBin(TH2D* hist, double cut_value, std::string axis) {
    // Check if the histogram is valid
    if (!hist) {
        throw std::invalid_argument("Histogram pointer is null.");
    }

    // Determine the axis to iterate over
    TAxis* target_axis = nullptr;
    if (axis == "x" || axis == "X" ) {
        target_axis = hist->GetXaxis();
    } else if (axis == "y" || axis == "Y") {
        target_axis = hist->GetYaxis();
    } else {
        throw std::invalid_argument("Invalid axis. Use 'x' or 'y'.");
    }

    // Find the bin closest to the cut value
    int n_bins = target_axis->GetNbins();

    if( cut_value > target_axis->GetBinCenter(n_bins) ) return n_bins;
    if( cut_value < target_axis->GetBinCenter(1) ) return 0;

    double min_diff = std::numeric_limits<double>::max(); // Initialize to a large value
    int closest_bin = -1;

    for (int bin = 1; bin <= n_bins; ++bin) { // Bins start at 1 in ROOT
        double bin_center = target_axis->GetBinCenter(bin);
        double diff = std::abs(bin_center - cut_value);
        if (diff < min_diff) {
            min_diff = diff;
            closest_bin = bin;
        }
    }

    return closest_bin; // Return the bin number
}
