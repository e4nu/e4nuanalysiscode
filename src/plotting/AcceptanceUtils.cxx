#include <iomanip>
#include <filesystem>
#include "plotting/AcceptanceUtils.h"

using namespace e4nu;
using namespace e4nu::plotting;

std::string plotting::Compute1DAcceptance(std::vector<std::string> mc_files, std::string observable, std::string title,
                                        std::string input_MC_location, std::string output_location, std::string output_file_name,
                                        std::map<std::string,std::vector<double>> cuts,
                                        std::string analysis_id, bool store_root )
{

  // Define trees
  std::vector<TFile *> files_mcrecoacc, files_mctrueacc;
  std::vector<TTree *> trees_mcrecoacc, trees_mctrueacc;

  // Define Hists
  // The _# correspond to histograms for each sector
  std::vector<TH1D *> hists_recoacc, hists_trueacc, hists_recoacc_0, hists_trueacc_0, hists_recoacc_1, hists_trueacc_1,
      hists_recoacc_2, hists_trueacc_2, hists_recoacc_3, hists_trueacc_3,
      hists_recoacc_4, hists_trueacc_4, hists_recoacc_5, hists_trueacc_5;
  std::vector<std::vector<TH1D *>> hists_recoacc_slices, hists_trueacc_slices;
  std::vector<TTree *> trees;
  std::vector<TH1D *> hists, ratios, ratios_0, ratios_1, ratios_2, ratios_3, ratios_4, ratios_5;
  std::vector<std::vector<TH1D *>> ratios_slices;
  std::vector<double> binning;
  // Get energy from tree to define range
  double BeamE;

  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    files_mcrecoacc.push_back(new TFile((input_MC_location + mc_files[i] + "_truereco.root").c_str(), "ROOT"));
    files_mctrueacc.push_back(new TFile((input_MC_location + mc_files[i] + "_true.root").c_str(), "ROOT"));
    if (!files_mcrecoacc[i])
    {
      std::cout << "ERROR: the " << mc_files[i] << "_truereco.root does not exist." << std::endl;
      return "";
    }
    if (!files_mctrueacc[i])
    {
      std::cout << "ERROR: the " << mc_files[i] << "_true.root  does not exist." << std::endl;
      return "";
    }
    trees_mcrecoacc.push_back((TTree *)files_mcrecoacc[i]->Get("MCCLAS6Tree"));
    trees_mctrueacc.push_back((TTree *)files_mctrueacc[i]->Get("MCCLAS6Tree"));
    if (!trees_mctrueacc[i] || !trees_mcrecoacc[i])
    {
      std::cout << "ERROR: the threes do not exist." << std::endl;
      return "";
    }

    trees_mctrueacc[0]->SetBranchAddress("BeamE", &BeamE);
    trees_mctrueacc[0]->GetEntry(0);
    binning = plotting::GetBinning(observable, BeamE, analysis_id);
    if( binning.size() == 0 ){
      std::cout << " ERROR: Binning is zero! Exiting... "<<std::endl;
      exit(0);
    }

    hists_recoacc.push_back(new TH1D(("Reco MC ACC Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_trueacc.push_back(new TH1D(("True MC ACC Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_recoacc_0.push_back(new TH1D(("Reco MC ACC Sector  0 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_trueacc_0.push_back(new TH1D(("True MC ACC Sector  0 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_recoacc_1.push_back(new TH1D(("Reco MC ACC Sector  1 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_trueacc_1.push_back(new TH1D(("True MC ACC Sector  1 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_recoacc_2.push_back(new TH1D(("Reco MC ACC Sector  2 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_trueacc_2.push_back(new TH1D(("True MC ACC Sector  2 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_recoacc_3.push_back(new TH1D(("Reco MC ACC Sector  3 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_trueacc_3.push_back(new TH1D(("True MC ACC Sector  3 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_recoacc_4.push_back(new TH1D(("Reco MC ACC Sector  4 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_trueacc_4.push_back(new TH1D(("True MC ACC Sector  4 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_recoacc_5.push_back(new TH1D(("Reco MC ACC Sector  5 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_trueacc_5.push_back(new TH1D(("True MC ACC Sector  5 Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));

    unsigned int initial_size_trees = trees.size();
    unsigned int initial_size_hists = hists.size();
    trees.push_back(trees_mcrecoacc[i]);
    trees.push_back(trees_mctrueacc[i]);
    hists.push_back(hists_recoacc[i]);   // 0
    hists.push_back(hists_trueacc[i]);   // 1
    hists.push_back(hists_recoacc_0[i]); // 2
    hists.push_back(hists_trueacc_0[i]); // 3
    hists.push_back(hists_recoacc_1[i]); // 4
    hists.push_back(hists_trueacc_1[i]); // 5
    hists.push_back(hists_recoacc_2[i]); // 6
    hists.push_back(hists_trueacc_2[i]); // 7
    hists.push_back(hists_recoacc_3[i]); // 8
    hists.push_back(hists_trueacc_3[i]); // 9
    hists.push_back(hists_recoacc_4[i]); // 10
    hists.push_back(hists_trueacc_4[i]); // 11
    hists.push_back(hists_recoacc_5[i]); // 12
    hists.push_back(hists_trueacc_5[i]); // 13

    //Set condition for new hists
    for (unsigned int id = initial_size_hists; id < hists.size(); id++)
    {
      hists[id]->Sumw2();
    }

    // Observables definition in Plotting Utils
    for (unsigned int j = initial_size_trees; j < trees.size(); ++j)
    {
      plotting::SetAnalysisBranch( trees[j] ) ;

      for (int k = 0; k < NEntries; ++k)
      {
        trees[j]->GetEntry(k);
        double content = 0;
        // Weight is the total weight devided by the number of entries.
        // This ensures that we get the same results even if we run less radiated events.
        double w = TotWeight / (double) InitialNEvents;
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

        // Fill the per Sector  histogram
        hists[2 * (ElectronSector + 1) + (j - initial_size_trees) + initial_size_hists]->Fill(content, w);
        hists[2 * (ElectronSector + 1) + (j - initial_size_trees) + initial_size_hists]->SetLineWidth(3);

        hists[j + initial_size_hists - initial_size_trees]->Fill(content, w);
        hists[j + initial_size_hists - initial_size_trees]->SetLineWidth(3);
      }
    }

    ratios.push_back((TH1D *)hists_trueacc[i]->Clone());
    ratios[i]->Divide(hists_recoacc[i]);
    ratios[i]->SetName(("Acceptance_model_" + std::to_string(i)).c_str());
    StandardFormat(ratios[i], title, kBlack + i + 1, 2 + i, observable);
    ratios[i]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
    ratios[i]->GetYaxis()->SetTitle("Acceptance correction");

    ratios_0.push_back((TH1D *)hists_trueacc_0[i]->Clone());
    ratios_0[i]->Divide(hists_recoacc_0[i]);
    ratios_0[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_0").c_str());
    StandardFormat(ratios_0[i], title, kOrange + 1 + i, 2 + i, observable);
    ratios_0[i]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
    ratios_0[i]->GetYaxis()->SetTitle("Acceptance correction e-Sector  0");

    ratios_1.push_back((TH1D *)hists_trueacc_1[i]->Clone());
    ratios_1[i]->Divide(hists_recoacc_1[i]);
    ratios_1[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_1").c_str());
    StandardFormat(ratios_1[i], title, kPink + 4 - i, 2 + i, observable);
    ratios_1[i]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
    ratios_1[i]->GetYaxis()->SetTitle("Acceptance correction e-Sector 1");

    ratios_2.push_back((TH1D *)hists_trueacc_2[i]->Clone());
    ratios_2[i]->Divide(hists_recoacc_2[i]);
    ratios_2[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_2").c_str());
    StandardFormat(ratios_2[i], title, kViolet + 5 - i, 2 + i, observable);
    ratios_2[i]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
    ratios_2[i]->GetYaxis()->SetTitle("Acceptance correction e-Sector 2");

    ratios_3.push_back((TH1D *)hists_trueacc_3[i]->Clone());
    ratios_3[i]->Divide(hists_recoacc_3[i]);
    ratios_3[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_3").c_str());
    StandardFormat(ratios_3[i], title, kAzure - 5 + i, 2 + i, observable);
    ratios_3[i]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
    ratios_3[i]->GetYaxis()->SetTitle("Acceptance correction e-Sector 3");

    ratios_4.push_back((TH1D *)hists_trueacc_4[i]->Clone());
    ratios_4[i]->Divide(hists_recoacc_4[i]);
    ratios_4[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_4").c_str());
    StandardFormat(ratios_4[i], title, kTeal - 7 - i, 2 + i, observable);
    ratios_4[i]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
    ratios_4[i]->GetYaxis()->SetTitle("Acceptance correction e-Sector 4");

    ratios_5.push_back((TH1D *)hists_trueacc_5[i]->Clone());
    ratios_5[i]->Divide(hists_recoacc_5[i]);
    ratios_5[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_5").c_str());
    StandardFormat(ratios_5[i], title, kGreen - 3 - i, 2 + i, observable);
    ratios_5[i]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
    ratios_5[i]->GetYaxis()->SetTitle("Acceptance correction e-Sector 5");

    std::vector<TH1D *> temp_ratios_slices;
    if (hists_trueacc_slices.size() != 0)
    {
      for (unsigned int l = 0; l < hists_trueacc_slices[i].size(); ++l)
      {
        temp_ratios_slices.push_back((TH1D *)hists_trueacc_slices[i][l]->Clone());
        temp_ratios_slices[l]->Divide(hists_recoacc_slices[i][l]);
        StandardFormat(temp_ratios_slices[l], title, kGreen - 3 - i, 2 + i, observable);
        std::string name = "Acceptance for slice " + std::to_string(l);
        temp_ratios_slices[l]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
        temp_ratios_slices[l]->GetYaxis()->SetTitle(name.c_str());
      }
    }
    ratios_slices.push_back(temp_ratios_slices);
  }

  TH1D *ratio = (TH1D *)ratios[0]->Clone();
  ratio->SetName("Acceptance");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio->Add(ratios[i]);
  }
  ratio->Scale(1. / mc_files.size());

  // We want to store the ratio before smoothing to account for this as an error
  TH1D *ratio_aSmooth = (TH1D *)ratio->Clone();
  ratio_aSmooth->Smooth(1);

  // Compute Acceptance error from model variation
  for (unsigned int i = 1; i < ratio->GetNbinsX() + 1; ++i)
  {
    double bin_cont_max = 0;
    double bin_cont_min = 999;

    double error_smoothing_2 = pow(ratio->GetBinContent(i) - ratio_aSmooth->GetBinContent(i), 2) / 2.;
    double error_stat_2 = pow(ratio->GetBinError(i), 2);
    double error_model_2 = 0;
    for (unsigned int j = 0; j < mc_files.size(); ++j)
    {
      if (ratios[j]->GetBinContent(i) > bin_cont_max)
      {
        bin_cont_max = ratios[j]->GetBinContent(i);
      }
      if (ratios[j]->GetBinContent(i) < bin_cont_min)
      {
        bin_cont_min = ratios[j]->GetBinContent(i);
      }
    }
    // Compute the error assuming a uniform distribution
    error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;

    // Adding all errors together
    ratio->SetBinError(i, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
  }
  StandardFormat(ratio, title, kBlack, 1, observable);
  ratio->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
  ratio->GetYaxis()->SetTitle("Acceptance correction");


  TH1D *ratio_0 = (TH1D *)ratios_0[0]->Clone();
  ratio_0->SetName("Acceptance_0");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_0->Add(ratios_0[i]);
  }
  ratio_0->Scale(1. / mc_files.size());
  TH1D *ratio_aSmooth_0 = (TH1D *)ratio_0->Clone();
  ratio_aSmooth_0->Smooth(3);

  // Compute Acceptance error from model variation
  for (unsigned int i = 1; i < ratio_0->GetNbinsX()+1; ++i)
  {
    double bin_cont_max = 0;
    double bin_cont_min = 999;
    double error_smoothing_2 = pow(ratio_0->GetBinContent(i) - ratio_aSmooth_0->GetBinContent(i), 2) / 12.;
    double error_stat_2 = pow(ratio_0->GetBinError(i), 2);
    double error_model_2 = 0;
    for (unsigned int j = 0; j < mc_files.size(); ++j)
    {
      if (ratios_0[j]->GetBinContent(i) > bin_cont_max)
        bin_cont_max = ratios_0[j]->GetBinContent(i);
      if (ratios_0[j]->GetBinContent(i) < bin_cont_min)
        bin_cont_min = ratios_0[j]->GetBinContent(i);
    }
    error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
    ratio_0->SetBinError(i, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
  }
  StandardFormat(ratio_0, title, kOrange + 1, 1, observable);
  ratio_0->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
  ratio_0->GetYaxis()->SetTitle("Acceptance correction e-Sector 0");

  TH1D *ratio_1 = (TH1D *)ratios_1[0]->Clone();
  ratio_1->SetName("Acceptance_1");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_1->Add(ratios_1[i]);
  }
  ratio_1->Scale(1. / mc_files.size());
  TH1D *ratio_aSmooth_1 = (TH1D *)ratio_1->Clone();
  ratio_aSmooth_1->Smooth(3);

  // Compute Acceptance error from model variation
  for (unsigned int i = 1; i < ratio_1->GetNbinsX()+1; ++i)
  {
    double bin_cont_max = 0;
    double bin_cont_min = 999;
    double error_smoothing_2 = pow(ratio_1->GetBinContent(i) - ratio_aSmooth_1->GetBinContent(i), 2) / 12.;
    double error_stat_2 = pow(ratio_1->GetBinError(i), 2);
    double error_model_2 = 0;
    for (unsigned int j = 0; j < mc_files.size(); ++j)
    {
      if (ratios_1[j]->GetBinContent(i) > bin_cont_max)
        bin_cont_max = ratios_1[j]->GetBinContent(i);
      if (ratios_1[j]->GetBinContent(i) < bin_cont_min)
        bin_cont_min = ratios_1[j]->GetBinContent(i);
    }
    error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
    ratio_1->SetBinError(i, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
  }
  StandardFormat(ratio_1, title, kPink + 4, 1, observable);
  ratio_1->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
  ratio_1->GetYaxis()->SetTitle("Acceptance correction e-Sector 1");

  TH1D *ratio_2 = (TH1D *)ratios_2[0]->Clone();
  ratio_2->SetName("Acceptance_2");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_2->Add(ratios_2[i]);
  }
  ratio_2->Scale(1. / mc_files.size());
  TH1D *ratio_aSmooth_2 = (TH1D *)ratio_2->Clone();
  ratio_aSmooth_2->Smooth(3);

  // Compute Acceptance error from model variation
  for (unsigned int i = 1; i < ratio_2->GetNbinsX()+1; ++i)
  {
    double bin_cont_max = 0;
    double bin_cont_min = 999;
    double error_smoothing_2 = pow(ratio_2->GetBinContent(i) - ratio_aSmooth_2->GetBinContent(i), 2) / 12.;
    double error_stat_2 = pow(ratio_2->GetBinError(i), 2);
    double error_model_2 = 0;
    for (unsigned int j = 0; j < mc_files.size(); ++j)
    {
      if (ratios_2[j]->GetBinContent(i) > bin_cont_max)
        bin_cont_max = ratios_2[j]->GetBinContent(i);
      if (ratios_2[j]->GetBinContent(i) < bin_cont_min)
        bin_cont_min = ratios_2[j]->GetBinContent(i);
    }
    error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12. ;
    ratio_2->SetBinError(i, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
  }
  StandardFormat(ratio_2, title, kViolet + 5, 1, observable);
  ratio_2->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
  ratio_2->GetYaxis()->SetTitle("Acceptance correction e-Sector 2");

  TH1D *ratio_3 = (TH1D *)ratios_3[0]->Clone();
  ratio_3->SetName("Acceptance_3");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_3->Add(ratios_3[i]);
  }
  ratio_3->Scale(1. / mc_files.size());
  TH1D *ratio_aSmooth_3 = (TH1D *)ratio_3->Clone();
  ratio_aSmooth_3->Smooth(3);

  // Compute Acceptance error from model variation
  for (unsigned int i = 1; i < ratio_3->GetNbinsX()+1; ++i)
  {
    double bin_cont_max = 0;
    double bin_cont_min = 999;
    double error_smoothing_2 = pow(ratio_3->GetBinContent(i) - ratio_aSmooth_3->GetBinContent(i), 2) / 12.;
    double error_stat_2 = pow(ratio_3->GetBinError(i), 2);
    double error_model_2 = 0;
    for (unsigned int j = 0; j < mc_files.size(); ++j)
    {
      if (ratios_3[j]->GetBinContent(i) > bin_cont_max)
        bin_cont_max = ratios_3[j]->GetBinContent(i);
      if (ratios_3[j]->GetBinContent(i) < bin_cont_min)
        bin_cont_min = ratios_3[j]->GetBinContent(i);
    }
    error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
    ratio_3->SetBinError(i, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
  }
  StandardFormat(ratio_3, title, kAzure - 5, 1, observable);
  ratio_3->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
  ratio_3->GetYaxis()->SetTitle("Acceptance correction e-Sector 3");

  TH1D *ratio_4 = (TH1D *)ratios_4[0]->Clone();
  ratio_4->SetName("Acceptance_4");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_4->Add(ratios_4[i]);
  }
  ratio_4->Scale(1. / mc_files.size());
  TH1D *ratio_aSmooth_4 = (TH1D *)ratio_4->Clone();
  ratio_aSmooth_4->Smooth(3);

  // Compute Acceptance error from model variation
  for (unsigned int i = 1; i < ratio_4->GetNbinsX()+1; ++i)
  {
    double bin_cont_max = 0;
    double bin_cont_min = 999;
    double error_smoothing_2 = pow(ratio_4->GetBinContent(i) - ratio_aSmooth_4->GetBinContent(i), 2) / 12.;
    double error_stat_2 = pow(ratio_4->GetBinError(i), 2);
    double error_model_2 = 0;
    for (unsigned int j = 0; j < mc_files.size(); ++j)
    {
      if (ratios_4[j]->GetBinContent(i) > bin_cont_max)
        bin_cont_max = ratios_4[j]->GetBinContent(i);
      if (ratios_4[j]->GetBinContent(i) < bin_cont_min)
        bin_cont_min = ratios_4[j]->GetBinContent(i);
    }
    error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
    ratio_4->SetBinError(i, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
  }
  StandardFormat(ratio_4, title, kTeal - 7, 1, observable);
  ratio_4->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
  ratio_4->GetYaxis()->SetTitle("Acceptance correction e-Sector 4");

  TH1D *ratio_5 = (TH1D *)ratios_5[0]->Clone();
  ratio_5->SetName("Acceptance_5");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_5->Add(ratios_5[i]);
  }
  ratio_5->Scale(1. / mc_files.size());
  TH1D *ratio_aSmooth_5 = (TH1D *)ratio_5->Clone();
  ratio_aSmooth_5->Smooth(3);
  // Compute Acceptance error from model variation
  for (unsigned int i = 1; i < ratio_5->GetNbinsX()+1; ++i)
  {
    double bin_cont_max = 0;
    double bin_cont_min = 999;
    double error_smoothing_2 = pow(ratio_5->GetBinContent(i) - ratio_aSmooth_5->GetBinContent(i), 2) / 12.;
    double error_stat_2 = pow(ratio_5->GetBinError(i), 2);
    double error_model_2 = 0;
    for (unsigned int j = 0; j < mc_files.size(); ++j)
    {
      if (ratios_5[j]->GetBinContent(i) > bin_cont_max)
        bin_cont_max = ratios_5[j]->GetBinContent(i);
      if (ratios_5[j]->GetBinContent(i) < bin_cont_min)
        bin_cont_min = ratios_5[j]->GetBinContent(i);
    }
    error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
    ratio_5->SetBinError(i, sqrt(error_stat_2 + error_model_2 + error_smoothing_2));
  }
  StandardFormat(ratio_5, title, kGreen - 3, 1, observable);
  ratio_5->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
  ratio_5->GetYaxis()->SetTitle("Acceptance correction e-Sector 5");

  std::vector<TH1D *> ratio_slices;
  if (ratios_slices.size() != 0)
  {
    for (unsigned l = 0; l < ratios_slices[0].size(); ++l)
    {
      std::string name = "MC Acceptance for ";

      TH1D *temp_slice_ratio = (TH1D *)ratios_slices[0][l]->Clone();
      temp_slice_ratio->SetName(("Acceptance_Slice_" + std::to_string(l)).c_str());
      for (unsigned int i = 1; i < mc_files.size(); ++i)
      {
        temp_slice_ratio->Add(ratios_slices[i][l]);
      }
      temp_slice_ratio->Scale(1. / mc_files.size());
      StandardFormat(temp_slice_ratio, title, kGreen - 3, 1, observable);
      temp_slice_ratio->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
      temp_slice_ratio->GetYaxis()->SetTitle("Acceptance");
      temp_slice_ratio->SetTitle((name + ". Slice " + std::to_string(l)).c_str());
      ratio_slices.push_back(temp_slice_ratio);
    }
  }
  std::string output_name = output_file_name + "_acceptance_correction_" + observable;
  std::string acc_file = "/AcceptanceFiles/" + output_name;

  std::filesystem::path acceptance_path{(output_location + "/AcceptanceFiles").c_str()};
  if (!std::filesystem::exists(acceptance_path))
    std::filesystem::create_directory(acceptance_path);

  TFile outputFile((output_location + acc_file + ".root").c_str(), "RECREATE");

  TCanvas *c_1 = new TCanvas("c_1", "c_1", 200, 10, 700, 500);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Store total contribution (averaged)
  ratio->Write();
  ratio_0->Write();
  ratio_1->Write();
  ratio_2->Write();
  ratio_3->Write();
  ratio_4->Write();
  ratio_5->Write();

  // Store per model
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios_0[i]->Write();
    ratios_1[i]->Write();
    ratios_2[i]->Write();
    ratios_3[i]->Write();
    ratios_4[i]->Write();
    ratios_5[i]->Write();
  }

  for (unsigned int i = 0; i < ratio_slices.size(); ++i)
  {
    ratio_slices[i]->Write();
  }

  // Plot it
  ratio->Draw("hist err P");
  ratio->SetMarkerStyle(8);
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios[i]->Draw("hist err same");
    ratios[i]->Write();
  }
  ratio->Draw("hist err P same");
  // teff->Draw("AP");

  if (store_root)
    c_1->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_total.root").c_str());
  c_1->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_total.pdf").c_str());
  delete c_1;

  // Draw total xsec per sectors
  TCanvas *c_sector_2 = new TCanvas("c_sector_2", "c_sector_2", 200, 10, 700, 500);
  c_sector_2->cd();
  TPad *pad_sector = new TPad("pad1", "", 0, 0, 1, 1);
  pad_sector->Draw();
  pad_sector->cd();
  pad_sector->SetBottomMargin(0.15);
  pad_sector->SetLeftMargin(0.15);
  pad_sector->Divide(3, 2);

  TPad *pad_sector_0 = (TPad *)pad_sector->cd(1);
  pad_sector_0->cd();
  pad_sector_0->SetBottomMargin(0.15);
  pad_sector_0->SetLeftMargin(0.15);
  ratio_0->GetYaxis()->SetTitleOffset(1.2);
  ratio_0->SetMarkerStyle(8);
  ratio_0->Draw("hist err P ");
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios_0[i]->Draw("hist err P same");
  }

  TPad *pad_sector_1 = (TPad *)pad_sector->cd(2);
  pad_sector_1->cd();
  pad_sector_1->SetBottomMargin(0.15);
  pad_sector_1->SetLeftMargin(0.15);
  ratio_1->GetYaxis()->SetTitleOffset(1.2);
  ratio_1->SetMarkerStyle(8);
  ratio_1->Draw("hist err P ");
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios_1[i]->Draw("hist err P same");
  }

  TPad *pad_sector_2 = (TPad *)pad_sector->cd(3);
  pad_sector_2->cd();
  pad_sector_2->SetBottomMargin(0.15);
  pad_sector_2->SetLeftMargin(0.15);
  ratio_2->GetYaxis()->SetTitleOffset(1.2);
  ratio_2->SetMarkerStyle(8);
  ratio_2->Draw("hist err P");
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios_2[i]->Draw("hist err P same");
  }

  TPad *pad_sector_3 = (TPad *)pad_sector->cd(4);
  pad_sector_3->cd();
  pad_sector_3->SetBottomMargin(0.15);
  pad_sector_3->SetLeftMargin(0.15);
  ratio_3->GetYaxis()->SetTitleOffset(1.2);
  ratio_3->SetMarkerStyle(8);
  ratio_3->Draw("hist err P ");
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios_3[i]->Draw("hist err P same");
  }

  TPad *pad_sector_4 = (TPad *)pad_sector->cd(5);
  pad_sector_4->cd();
  pad_sector_4->SetBottomMargin(0.15);
  pad_sector_4->SetLeftMargin(0.15);
  ratio_4->GetYaxis()->SetTitleOffset(1.2);
  ratio_4->SetMarkerStyle(8);
  ratio_4->Draw("hist err P");
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios_4[i]->Draw("hist err P same");
  }

  TPad *pad_sector_5 = (TPad *)pad_sector->cd(6);
  pad_sector_5->cd();
  pad_sector_5->SetBottomMargin(0.15);
  pad_sector_5->SetLeftMargin(0.15);
  ratio_5->GetYaxis()->SetTitleOffset(1.2);
  ratio_5->SetMarkerStyle(8);
  ratio_5->Draw("hist err P");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratios_5[i]->Draw("hist err P same");
  }

  if (store_root)
    c_sector_2->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_persector.root").c_str());
  c_sector_2->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_persector.pdf").c_str());
  delete c_sector_2;

  for (size_t i = 0; i < files_mcrecoacc.size(); i++)
  {
    delete files_mcrecoacc[i];
    delete files_mctrueacc[i];
  }
  outputFile.Close();
  return acc_file;
}


std::string plotting::Compute2DAcceptance(std::vector<std::string> mc_files, std::string x_observable, std::string y_observable, std::string title,
                                        std::string input_MC_location, std::string output_location, std::string output_file_name,
                                        std::map<std::string,std::vector<double>> cuts,
                                        std::string analysis_id, bool store_root )
{

  // Define trees
  std::vector<TFile *> files_mcrecoacc, files_mctrueacc;
  std::vector<TTree *> trees_mcrecoacc, trees_mctrueacc;

  // Define Hists
  // The _# correspond to histograms for each sector
  std::vector<TH2D *> hists_recoacc, hists_trueacc, hists_recoacc_0, hists_trueacc_0, hists_recoacc_1, hists_trueacc_1,
      hists_recoacc_2, hists_trueacc_2, hists_recoacc_3, hists_trueacc_3,
      hists_recoacc_4, hists_trueacc_4, hists_recoacc_5, hists_trueacc_5;
  std::vector<std::vector<TH2D *>> hists_recoacc_slices, hists_trueacc_slices;
  std::vector<TTree *> trees;
  std::vector<TH2D *> hists, ratios, ratios_0, ratios_1, ratios_2, ratios_3, ratios_4, ratios_5;
  std::vector<double> x_binning;
  std::vector<double> y_binning;
  // Get energy from tree to define range
  double BeamE;

  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    files_mcrecoacc.push_back(new TFile((input_MC_location + mc_files[i] + "_truereco.root").c_str(), "ROOT"));
    files_mctrueacc.push_back(new TFile((input_MC_location + mc_files[i] + "_true.root").c_str(), "ROOT"));
    if (!files_mcrecoacc[i])
    {
      std::cout << "ERROR: the " << mc_files[i] << "_truereco.root does not exist." << std::endl;
      return "";
    }
    if (!files_mctrueacc[i])
    {
      std::cout << "ERROR: the " << mc_files[i] << "_true.root  does not exist." << std::endl;
      return "";
    }
    trees_mcrecoacc.push_back((TTree *)files_mcrecoacc[i]->Get("MCCLAS6Tree"));
    trees_mctrueacc.push_back((TTree *)files_mctrueacc[i]->Get("MCCLAS6Tree"));
    if (!trees_mctrueacc[i] || !trees_mcrecoacc[i])
    {
      std::cout << "ERROR: the threes do not exist." << std::endl;
      return "";
    }

    trees_mctrueacc[0]->SetBranchAddress("BeamE", &BeamE);
    trees_mctrueacc[0]->GetEntry(0);
    x_binning = plotting::GetBinning(x_observable, BeamE, analysis_id);
    if( x_binning.size() == 0 ){
      std::cout << " ERROR: Binning is zero! Exiting... "<<std::endl;
      exit(0);
    }

    y_binning = plotting::GetBinning(y_observable, BeamE, analysis_id);
    if( y_binning.size() == 0 ){
      std::cout << " ERROR: Binning is zero! Exiting... "<<std::endl;
      exit(0);
    }

    hists_recoacc.push_back(new TH2D(("Reco MC ACC Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_trueacc.push_back(new TH2D(("True MC ACC Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_recoacc_0.push_back(new TH2D(("Reco MC ACC Sector  0 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_trueacc_0.push_back(new TH2D(("True MC ACC Sector  0 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_recoacc_1.push_back(new TH2D(("Reco MC ACC Sector  1 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_trueacc_1.push_back(new TH2D(("True MC ACC Sector  1 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_recoacc_2.push_back(new TH2D(("Reco MC ACC Sector  2 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_trueacc_2.push_back(new TH2D(("True MC ACC Sector  2 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_recoacc_3.push_back(new TH2D(("Reco MC ACC Sector  3 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_trueacc_3.push_back(new TH2D(("True MC ACC Sector  3 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_recoacc_4.push_back(new TH2D(("Reco MC ACC Sector  4 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_trueacc_4.push_back(new TH2D(("True MC ACC Sector  4 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_recoacc_5.push_back(new TH2D(("Reco MC ACC Sector  5 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));
    hists_trueacc_5.push_back(new TH2D(("True MC ACC Sector  5 Model " + std::to_string(i)).c_str(), "", x_binning.size() - 1, &x_binning[0], y_binning.size() - 1, &y_binning[0]));

    unsigned int initial_size_trees = trees.size();
    unsigned int initial_size_hists = hists.size();
    trees.push_back(trees_mcrecoacc[i]);
    trees.push_back(trees_mctrueacc[i]);
    hists.push_back(hists_recoacc[i]);   // 0
    hists.push_back(hists_trueacc[i]);   // 1
    hists.push_back(hists_recoacc_0[i]); // 2
    hists.push_back(hists_trueacc_0[i]); // 3
    hists.push_back(hists_recoacc_1[i]); // 4
    hists.push_back(hists_trueacc_1[i]); // 5
    hists.push_back(hists_recoacc_2[i]); // 6
    hists.push_back(hists_trueacc_2[i]); // 7
    hists.push_back(hists_recoacc_3[i]); // 8
    hists.push_back(hists_trueacc_3[i]); // 9
    hists.push_back(hists_recoacc_4[i]); // 10
    hists.push_back(hists_trueacc_4[i]); // 11
    hists.push_back(hists_recoacc_5[i]); // 12
    hists.push_back(hists_trueacc_5[i]); // 13

    // Set condition for new hists
    for (unsigned int id = initial_size_hists; id < hists.size(); id++)
    {
      hists[id]->Sumw2();
    }

    // Observables definition in Plotting Utils
    for (unsigned int j = initial_size_trees; j < trees.size(); ++j)
    {
      plotting::SetAnalysisBranch( trees[j] ) ;

      for (int k = 0; k < NEntries; ++k)
      {
        trees[j]->GetEntry(k);
        double content_x = GetObservable(x_observable);
        double content_y = GetObservable(y_observable);
        // Weight is the total weight devided by the number of entries.
        // This ensures that we get the same results even if we run less radiated events.
        double w = TotWeight / (double) InitialNEvents;

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

        // Fill the per Sector  histogram
        hists[2 * (ElectronSector + 1) + (j - initial_size_trees) + initial_size_hists]->Fill(content_x, content_y, w);
        hists[2 * (ElectronSector + 1) + (j - initial_size_trees) + initial_size_hists]->SetLineWidth(3);

        hists[j + initial_size_hists - initial_size_trees]->Fill(content_x, content_y, w);
        hists[j + initial_size_hists - initial_size_trees]->SetLineWidth(3);
      }
    }

    ratios.push_back((TH2D *)hists_trueacc[i]->Clone());
    ratios[i]->Divide(hists_recoacc[i]);
    ratios[i]->SetName(("Acceptance_model_" + std::to_string(i)).c_str());
    StandardFormat(ratios[i], title, kBlack + i + 1, 2 + i, x_observable, y_observable);
    ratios[i]->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios[i]->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
    ratios[i]->GetZaxis()->SetTitle("Acceptance correction");

    ratios_0.push_back((TH2D *)hists_trueacc_0[i]->Clone());
    ratios_0[i]->Divide(hists_recoacc_0[i]);
    ratios_0[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_0").c_str());
    StandardFormat(ratios_0[i], title, kOrange + 1 + i, 2 + i, x_observable, y_observable);
    ratios_0[i]->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios_0[i]->GetYaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios_0[i]->GetZaxis()->SetTitle("Acceptance correction e-Sector  0");

    ratios_1.push_back((TH2D *)hists_trueacc_1[i]->Clone());
    ratios_1[i]->Divide(hists_recoacc_1[i]);
    ratios_1[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_1").c_str());
    StandardFormat(ratios_1[i], title, kPink + 4 - i, 2 + i, x_observable, y_observable);
    ratios_1[i]->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios_1[i]->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
    ratios_1[i]->GetZaxis()->SetTitle("Acceptance correction e-Sector 1");

    ratios_2.push_back((TH2D *)hists_trueacc_2[i]->Clone());
    ratios_2[i]->Divide(hists_recoacc_2[i]);
    ratios_2[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_2").c_str());
    StandardFormat(ratios_2[i], title, kViolet + 5 - i, 2 + i, x_observable, y_observable);
    ratios_2[i]->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios_2[i]->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
    ratios_2[i]->GetZaxis()->SetTitle("Acceptance correction e-Sector 2");

    ratios_3.push_back((TH2D *)hists_trueacc_3[i]->Clone());
    ratios_3[i]->Divide(hists_recoacc_3[i]);
    ratios_3[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_3").c_str());
    StandardFormat(ratios_3[i], title, kAzure - 5 + i, 2 + i, x_observable, y_observable);
    ratios_3[i]->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios_3[i]->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
    ratios_3[i]->GetZaxis()->SetTitle("Acceptance correction e-Sector 3");

    ratios_4.push_back((TH2D *)hists_trueacc_4[i]->Clone());
    ratios_4[i]->Divide(hists_recoacc_4[i]);
    ratios_4[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_4").c_str());
    StandardFormat(ratios_4[i], title, kTeal - 7 - i, 2 + i, x_observable, y_observable);
    ratios_4[i]->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios_4[i]->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
    ratios_4[i]->GetZaxis()->SetTitle("Acceptance correction e-Sector 4");

    ratios_5.push_back((TH2D *)hists_trueacc_5[i]->Clone());
    ratios_5[i]->Divide(hists_recoacc_5[i]);
    ratios_5[i]->SetName(("Acceptance_model_" + std::to_string(i) + "_sector_5").c_str());
    StandardFormat(ratios_5[i], title, kGreen - 3 - i, 2 + i, x_observable, y_observable);
    ratios_5[i]->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios_5[i]->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
    ratios_5[i]->GetZaxis()->SetTitle("Acceptance correction e-Sector 5");
  }

  TH2D *ratio = (TH2D *)ratios[0]->Clone();
  ratio->SetName("Acceptance");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio->Add(ratios[i]);
  }
  ratio->Scale(1. / mc_files.size());

  // We want to store the ratio before smoothing to account for this as an error
  TH2D *ratio_aSmooth = (TH2D *)ratio->Clone();
  ratio_aSmooth->Smooth(1);

  // Compute Acceptance error from model variation
  for (int i = 1; i <= ratio->GetNbinsX(); ++i) {
    for (int j = 1; j <= ratio->GetNbinsY(); ++j) {
        double bin_cont_max = 0;
        double bin_cont_min = 999;
        double error_smoothing_2 = pow(ratio->GetBinContent(i, j) - ratio_aSmooth->GetBinContent(i, j), 2) / 12.;
        double error_stat_2 = pow(ratio->GetBinError(i, j), 2);
        double error_model_2 = 0;

        for (unsigned int k = 0; k < mc_files.size(); ++k) {
            if (ratios[k]->GetBinContent(i, j) > bin_cont_max)
                bin_cont_max = ratios[k]->GetBinContent(i, j);
            if (ratios[k]->GetBinContent(i, j) < bin_cont_min)
                bin_cont_min = ratios[k]->GetBinContent(i, j);
        }

        error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
        ratio->SetBinError(i, j, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
    }
  }
  // StandardFormat(ratio, title, kBlack, 1, observable);
  ratio->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
  ratio->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
  ratio->GetZaxis()->SetTitle("Acceptance correction");

  TH2D *ratio_0 = (TH2D *)ratios_0[0]->Clone();
  ratio_0->SetName("Acceptance_0");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_0->Add(ratios_0[i]);
  }
  ratio_0->Scale(1. / mc_files.size());
  TH2D *ratio_aSmooth_0 = (TH2D *)ratio_0->Clone();
  ratio_aSmooth_0->Smooth(3);

  // Compute Acceptance error from model variation
  for (int i = 1; i <= ratio_0->GetNbinsX(); ++i) {
    for (int j = 1; j <= ratio_0->GetNbinsY(); ++j) {
        double bin_cont_max = 0;
        double bin_cont_min = 999;
        double error_smoothing_2 = pow(ratio_0->GetBinContent(i, j) - ratio_aSmooth_0->GetBinContent(i, j), 2) / 12.;
        double error_stat_2 = pow(ratio_0->GetBinError(i, j), 2);
        double error_model_2 = 0;

        for (unsigned int k = 0; k < mc_files.size(); ++k) {
            if (ratios_0[k]->GetBinContent(i, j) > bin_cont_max)
                bin_cont_max = ratios_0[k]->GetBinContent(i, j);
            if (ratios_0[k]->GetBinContent(i, j) < bin_cont_min)
                bin_cont_min = ratios_0[k]->GetBinContent(i, j);
        }

        error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
        ratio_0->SetBinError(i, j, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
    }
  }
  StandardFormat(ratio_0, title, kOrange + 1, 1, x_observable, y_observable);
  ratio_0->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
  ratio_0->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
  ratio_0->GetZaxis()->SetTitle("Acceptance correction e-Sector 0");

  TH2D *ratio_1 = (TH2D *)ratios_1[0]->Clone();
  ratio_1->SetName("Acceptance_1");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_1->Add(ratios_1[i]);
  }
  ratio_1->Scale(1. / mc_files.size());
  TH2D *ratio_aSmooth_1 = (TH2D *)ratio_1->Clone();
  ratio_aSmooth_1->Smooth(3);

  // Compute Acceptance error from model variation
  for (int i = 1; i <= ratio_1->GetNbinsX(); ++i) {
    for (int j = 1; j <= ratio_1->GetNbinsY(); ++j) {
        double bin_cont_max = 0;
        double bin_cont_min = 999;
        double error_smoothing_2 = pow(ratio_1->GetBinContent(i, j) - ratio_aSmooth_1->GetBinContent(i, j), 2) / 12.;
        double error_stat_2 = pow(ratio_1->GetBinError(i, j), 2);
        double error_model_2 = 0;

        for (unsigned int k = 0; k < mc_files.size(); ++k) {
            if (ratios_1[k]->GetBinContent(i, j) > bin_cont_max)
                bin_cont_max = ratios_1[k]->GetBinContent(i, j);
            if (ratios_1[k]->GetBinContent(i, j) < bin_cont_min)
                bin_cont_min = ratios_1[k]->GetBinContent(i, j);
        }

        error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
        ratio_1->SetBinError(i, j, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
    }
  }
  StandardFormat(ratio_1, title, kPink + 4, 1, x_observable, y_observable);
  ratio_1->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
  ratio_1->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
  ratio_1->GetZaxis()->SetTitle("Acceptance correction e-Sector 1");

  TH2D *ratio_2 = (TH2D *)ratios_2[0]->Clone();
  ratio_2->SetName("Acceptance_2");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_2->Add(ratios_2[i]);
  }
  ratio_2->Scale(1. / mc_files.size());
  TH2D *ratio_aSmooth_2 = (TH2D *)ratio_2->Clone();
  ratio_aSmooth_2->Smooth(3);

  // Compute Acceptance error from model variation
  for (int i = 1; i <= ratio_2->GetNbinsX(); ++i) {
    for (int j = 1; j <= ratio_2->GetNbinsY(); ++j) {
        double bin_cont_max = 0;
        double bin_cont_min = 999;
        double error_smoothing_2 = pow(ratio_2->GetBinContent(i, j) - ratio_aSmooth_2->GetBinContent(i, j), 2) / 12.;
        double error_stat_2 = pow(ratio_2->GetBinError(i, j), 2);
        double error_model_2 = 0;

        for (unsigned int k = 0; k < mc_files.size(); ++k) {
            if (ratios_2[k]->GetBinContent(i, j) > bin_cont_max)
                bin_cont_max = ratios_2[k]->GetBinContent(i, j);
            if (ratios_2[k]->GetBinContent(i, j) < bin_cont_min)
                bin_cont_min = ratios_2[k]->GetBinContent(i, j);
        }

        error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
        ratio_2->SetBinError(i, j, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
    }
  }
  StandardFormat(ratio_2, title, kViolet + 5, 1, x_observable, y_observable);
  ratio_2->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
  ratio_2->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
  ratio_2->GetZaxis()->SetTitle("Acceptance correction e-Sector 2");

  TH2D *ratio_3 = (TH2D *)ratios_3[0]->Clone();
  ratio_3->SetName("Acceptance_3");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_3->Add(ratios_3[i]);
  }
  ratio_3->Scale(1. / mc_files.size());
  TH2D *ratio_aSmooth_3 = (TH2D *)ratio_3->Clone();
  ratio_aSmooth_3->Smooth(3);

  // Compute Acceptance error from model variation
  for (int i = 1; i <= ratio_3->GetNbinsX(); ++i) {
    for (int j = 1; j <= ratio_3->GetNbinsY(); ++j) {
        double bin_cont_max = 0;
        double bin_cont_min = 999;
        double error_smoothing_2 = pow(ratio_3->GetBinContent(i, j) - ratio_aSmooth_3->GetBinContent(i, j), 2) / 12.;
        double error_stat_2 = pow(ratio_3->GetBinError(i, j), 2);
        double error_model_2 = 0;

        for (unsigned int k = 0; k < mc_files.size(); ++k) {
            if (ratios_3[k]->GetBinContent(i, j) > bin_cont_max)
                bin_cont_max = ratios_3[k]->GetBinContent(i, j);
            if (ratios_3[k]->GetBinContent(i, j) < bin_cont_min)
                bin_cont_min = ratios_3[k]->GetBinContent(i, j);
        }

        error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
        ratio_3->SetBinError(i, j, sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
    }
  }
  StandardFormat(ratio_3, title, kAzure - 5, 1, x_observable, y_observable);
  ratio_3->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
  ratio_3->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
  ratio_3->GetZaxis()->SetTitle("Acceptance correction e-Sector 3");

  TH2D *ratio_4 = (TH2D *)ratios_4[0]->Clone();
  ratio_4->SetName("Acceptance_4");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_4->Add(ratios_4[i]);
  }
  ratio_4->Scale(1. / mc_files.size());
  TH2D *ratio_aSmooth_4 = (TH2D *)ratio_4->Clone();
  ratio_aSmooth_4->Smooth(3);

  // Compute Acceptance error from model variation
  for (int i = 1; i <= ratio_4->GetNbinsX(); ++i) {
    for (int j = 1; j <= ratio_4->GetNbinsY(); ++j) {
        double bin_cont_max = 0;
        double bin_cont_min = 999;
        double error_smoothing_2 = pow(ratio_4->GetBinContent(i, j) - ratio_aSmooth_4->GetBinContent(i, j), 2) / 12.;
        double error_stat_2 = pow(ratio_4->GetBinError(i, j), 2);
        double error_model_2 = 0;

        for (unsigned int k = 0; k < mc_files.size(); ++k) {
            if (ratios_4[k]->GetBinContent(i, j) > bin_cont_max)
                bin_cont_max = ratios_4[k]->GetBinContent(i, j);
            if (ratios_4[k]->GetBinContent(i, j) < bin_cont_min)
                bin_cont_min = ratios_4[k]->GetBinContent(i, j);
        }

        error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
        ratio_4->SetBinError(i, j, 1E30); // sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
    }
  }
  StandardFormat(ratio_4, title, kTeal - 7, 1, x_observable,y_observable);
  ratio_4->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
  ratio_4->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
  ratio_4->GetZaxis()->SetTitle("Acceptance correction e-Sector 4");

  TH2D *ratio_5 = (TH2D *)ratios_5[0]->Clone();
  ratio_5->SetName("Acceptance_5");
  for (unsigned int i = 1; i < mc_files.size(); ++i)
  {
    ratio_5->Add(ratios_5[i]);
  }
  ratio_5->Scale(1. / mc_files.size());
  TH2D *ratio_aSmooth_5 = (TH2D *)ratio_5->Clone();
  ratio_aSmooth_5->Smooth(3);
  // Compute Acceptance error from model variation
  for (int i = 1; i <= ratio_5->GetNbinsX(); ++i) {
    for (int j = 1; j <= ratio_5->GetNbinsY(); ++j) {
        double bin_cont_max = 0;
        double bin_cont_min = 999;
        double error_smoothing_2 = pow(ratio_5->GetBinContent(i, j) - ratio_aSmooth_5->GetBinContent(i, j), 2) / 12.;
        double error_stat_2 = pow(ratio_0->GetBinError(i, j), 2);
        double error_model_2 = 0;

        for (unsigned int k = 0; k < mc_files.size(); ++k) {
            if (ratios_5[k]->GetBinContent(i, j) > bin_cont_max)
                bin_cont_max = ratios_5[k]->GetBinContent(i, j);
            if (ratios_5[k]->GetBinContent(i, j) < bin_cont_min)
                bin_cont_min = ratios_5[k]->GetBinContent(i, j);
        }

        error_model_2 = pow(bin_cont_max - bin_cont_min, 2) / 12.;
        ratio_5->SetBinError(i, j, 1E30); // sqrt(error_stat_2 + error_smoothing_2 + error_model_2));
    }
  }
  StandardFormat(ratio_5, title, kGreen - 3, 1, x_observable,y_observable);
  ratio_5->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
  ratio_5->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
  ratio_5->GetZaxis()->SetTitle("Acceptance correction e-Sector 5");
  std::string output_name = output_file_name + "_acceptance_correction_" + x_observable + "_vs_" + y_observable;
  std::string acc_file = "/AcceptanceFiles/" + output_name;

  std::filesystem::path acceptance_path{(output_location + "/AcceptanceFiles").c_str()};
  if (!std::filesystem::exists(acceptance_path))
    std::filesystem::create_directory(acceptance_path);

  TFile outputFile((output_location + acc_file + ".root").c_str(), "RECREATE");

  TCanvas *c_1 = new TCanvas("c_1", "c_1", 200, 10, 700, 500);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Store total contribution (averaged)
  ratio->Write();
  ratio_0->Write();
  ratio_1->Write();
  ratio_2->Write();
  ratio_3->Write();
  ratio_4->Write();
  ratio_5->Write();

  // Store per model
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios_0[i]->Write();
    ratios_1[i]->Write();
    ratios_2[i]->Write();
    ratios_3[i]->Write();
    ratios_4[i]->Write();
    ratios_5[i]->Write();
  }


  // Plot it
  ratio->Draw("COLZ");
  // ratio->SetMarkerStyle(8);
  for (unsigned int i = 0; i < mc_files.size(); ++i)
  {
    ratios[i]->Write();
  }

  if (store_root)
    c_1->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_total.root").c_str());
  c_1->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_total.pdf").c_str());
  delete c_1;

  // Draw total xsec per sectors
  TCanvas *c_sector_2 = new TCanvas("c_sector_2", "c_sector_2", 200, 10, 700, 500);
  c_sector_2->cd();
  TPad *pad_sector = new TPad("pad1", "", 0, 0, 1, 1);
  pad_sector->Draw();
  pad_sector->cd();
  pad_sector->SetBottomMargin(0.15);
  pad_sector->SetLeftMargin(0.15);
  pad_sector->Divide(3, 2);

  TPad *pad_sector_0 = (TPad *)pad_sector->cd(1);
  pad_sector_0->cd();
  pad_sector_0->SetBottomMargin(0.15);
  pad_sector_0->SetLeftMargin(0.15);
  ratio_0->GetYaxis()->SetTitleOffset(1.2);
  ratio_0->SetMarkerStyle(8);
  ratio_0->Draw("COLZ");

  TPad *pad_sector_1 = (TPad *)pad_sector->cd(2);
  pad_sector_1->cd();
  pad_sector_1->SetBottomMargin(0.15);
  pad_sector_1->SetLeftMargin(0.15);
  ratio_1->GetYaxis()->SetTitleOffset(1.2);
  ratio_1->SetMarkerStyle(8);
  ratio_1->Draw("COLZ");

  TPad *pad_sector_2 = (TPad *)pad_sector->cd(3);
  pad_sector_2->cd();
  pad_sector_2->SetBottomMargin(0.15);
  pad_sector_2->SetLeftMargin(0.15);
  ratio_2->GetYaxis()->SetTitleOffset(1.2);
  ratio_2->SetMarkerStyle(8);
  ratio_2->Draw("COLZ");

  TPad *pad_sector_3 = (TPad *)pad_sector->cd(4);
  pad_sector_3->cd();
  pad_sector_3->SetBottomMargin(0.15);
  pad_sector_3->SetLeftMargin(0.15);
  ratio_3->GetYaxis()->SetTitleOffset(1.2);
  ratio_3->SetMarkerStyle(8);
  ratio_3->Draw("COLZ");

  TPad *pad_sector_4 = (TPad *)pad_sector->cd(5);
  pad_sector_4->cd();
  pad_sector_4->SetBottomMargin(0.15);
  pad_sector_4->SetLeftMargin(0.15);
  ratio_4->GetYaxis()->SetTitleOffset(1.2);
  ratio_4->SetMarkerStyle(8);
  ratio_4->Draw("COLZ");

  TPad *pad_sector_5 = (TPad *)pad_sector->cd(6);
  pad_sector_5->cd();
  pad_sector_5->SetBottomMargin(0.15);
  pad_sector_5->SetLeftMargin(0.15);
  ratio_5->GetYaxis()->SetTitleOffset(1.2);
  ratio_5->SetMarkerStyle(8);
  ratio_5->Draw("COLZ");

  if (store_root)
    c_sector_2->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_2D_persector.root").c_str());
  c_sector_2->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_2D_persector.pdf").c_str());
  delete c_sector_2;

  for (size_t i = 0; i < files_mcrecoacc.size(); i++)
  {
    delete files_mcrecoacc[i];
    delete files_mctrueacc[i];
  }
  outputFile.Close();
  return acc_file;
}

std::string plotting::Compute1DRadCorr(std::vector<std::string> mc_files, std::vector<std::string> rad_files,
                                       std::string observable, std::string title,
                                       std::string input_MC_location, std::string output_location, std::string output_file_name,
                                       std::map<std::string,std::vector<double>> cuts,
                                       std::string analysis_id, bool store_root )
{

  // Define trees
  std::vector<TFile *> files_mctrueacc, files_mcradcorr;
  std::vector<TTree *> trees_mctrueacc, trees_mcradcorr;

  // Define Hists
  // The _# correspond to histograms for each sector
  std::vector<TH1D *> ratios, hists_true, hists_radcorr;
  std::vector<double> binning;
  // Get energy from tree to define range
  double BeamE;

  for (unsigned int i = 0; i < rad_files.size(); ++i)
  {
    files_mctrueacc.push_back(new TFile((input_MC_location + mc_files[i] + "_true.root").c_str(), "ROOT"));
    files_mcradcorr.push_back(new TFile((input_MC_location + rad_files[i] + "_true_radcorr.root").c_str(), "ROOT"));
    if (!files_mctrueacc[i])
    {
      std::cout << "ERROR: the " << mc_files[i] << "_true.root  does not exist." << std::endl;
      return "";
    }
    if (!files_mcradcorr[i])
    {
      std::cout << "ERROR: the " << mc_files[i] << "_true_radcorr.root does not exist." << std::endl;
      return "";
    }
    trees_mctrueacc.push_back((TTree *)files_mctrueacc[i]->Get("MCCLAS6Tree"));
    trees_mcradcorr.push_back((TTree *)files_mcradcorr[i]->Get("MCCLAS6Tree"));
    if (!trees_mctrueacc[i] || !trees_mcradcorr[i])
    {
      std::cout << "ERROR: the threes do not exist." << std::endl;
      return "";
    }

    trees_mctrueacc[0]->SetBranchAddress("BeamE", &BeamE);
    trees_mctrueacc[0]->GetEntry(0);
    binning = plotting::GetBinning(observable, BeamE, analysis_id);

    hists_radcorr.push_back(new TH1D(("Rad Corr MC ACC Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));
    hists_true.push_back(new TH1D(("True MC ACC Model " + std::to_string(i)).c_str(), "", binning.size() - 1, &binning[0]));

    hists_radcorr[i]->Sumw2();
    hists_true[i]->Sumw2();

    std::vector<TTree *> trees = {trees_mctrueacc[i], trees_mcradcorr[i]};
    for (unsigned int j = 0; j < trees.size(); ++j)
    {
      plotting::SetAnalysisBranch( trees[j] ) ;

      for (int k = 0; k < NEntries; ++k)
      {
        trees[j]->GetEntry(k);
        double content = 0;
        // Weight is the total weight devided by the number of entries.
        // This ensures that we get the same results even if we run less radiated events.
        double w = TotWeight / (double) InitialNEvents ;

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

        // Fill the per Sector  histogram
        if (j == 0)
          hists_true[i]->Fill(content, w);
        else if (j == 1) {
          hists_radcorr[i]->Fill(content, w);
        }

      }
    }

    ratios.push_back((TH1D *)hists_true[i]->Clone());
    ratios[i]->Divide(hists_radcorr[i]);
    ratios[i]->SetName(("RadCorrModel_" + std::to_string(i)).c_str());
    StandardFormat(ratios[i], title, kBlack + i + 1, 2 + i, observable);
    ratios[i]->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
    ratios[i]->GetYaxis()->SetTitle("Radiation correction");
  }

  // Calculate average correction
  TH1D *ratio = (TH1D *)ratios[0]->Clone();
  ratio->SetName("Acceptance");
  for (unsigned int i = 1; i < rad_files.size(); ++i)
  {
    ratio->Add(ratios[i]);
  }

  //ratio->Scale(1. / mc_files.size());
  StandardFormat(ratio, title, kBlack, 1, observable);
  ratio->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());
  ratio->GetYaxis()->SetTitle("Radiation Correction");

  std::string output_name = output_file_name + "_acceptance_correction_rad_" + observable;
  std::string acc_file = "/AcceptanceFiles/" + output_name;

  std::filesystem::path acceptance_path{(output_location + "/AcceptanceFiles").c_str()};
  if (!std::filesystem::exists(acceptance_path))
    std::filesystem::create_directory(acceptance_path);

  TFile outputFile((output_location + acc_file + ".root").c_str(), "RECREATE");

  TCanvas *c_1 = new TCanvas("c_1", "c_1", 200, 10, 700, 500);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Store total contribution (averaged)
  ratio->Write();

  // Plot it
  ratio->SetMarkerStyle(8);
  ratio->Draw("hist err P ");
  for (unsigned int i = 0; i < rad_files.size(); ++i)
  {
    ratios[i]->Draw("hist err same");
    ratios[i]->Write();
  }
  ratio->Draw("hist err P same");
  // teff->Draw("AP");

  if (store_root)
    c_1->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_total.root").c_str());
  c_1->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_total.pdf").c_str());

  delete c_1;
  return acc_file;
}


std::string plotting::Compute2DRadCorr(std::vector<std::string> mc_files, std::vector<std::string> rad_files,
                                        std::string x_observable, std::string y_observable, std::string title,
                                        std::string input_MC_location, std::string output_location, std::string output_file_name,
                                        std::map<std::string,std::vector<double>> cuts,
                                        std::string analysis_id, bool store_root )
{

  // Define trees
  std::vector<TFile *> files_mctrueacc, files_mcradcorr;
  std::vector<TTree *> trees_mctrueacc, trees_mcradcorr;

  // Define Hists
  // The _# correspond to histograms for each sector
  std::vector<TH2D *> ratios, hists_true, hists_radcorr;
  std::vector<double> binning_x, binning_y;
  // Get energy from tree to define range
  double BeamE;

  for (unsigned int i = 0; i < rad_files.size(); ++i)
  {
    files_mctrueacc.push_back(new TFile((input_MC_location + mc_files[i] + "_true.root").c_str(), "ROOT"));
    files_mcradcorr.push_back(new TFile((input_MC_location + rad_files[i] + "_true_radcorr.root").c_str(), "ROOT"));
    if (!files_mctrueacc[i])
    {
      std::cout << "ERROR: the " << mc_files[i] << "_true.root  does not exist." << std::endl;
      return "";
    }
    if (!files_mcradcorr[i])
    {
      std::cout << "ERROR: the " << mc_files[i] << "_true_radcorr.root does not exist." << std::endl;
      return "";
    }
    trees_mctrueacc.push_back((TTree *)files_mctrueacc[i]->Get("MCCLAS6Tree"));
    trees_mcradcorr.push_back((TTree *)files_mcradcorr[i]->Get("MCCLAS6Tree"));
    if (!trees_mctrueacc[i] || !trees_mcradcorr[i])
    {
      std::cout << "ERROR: the threes do not exist." << std::endl;
      return "";
    }

    trees_mctrueacc[0]->SetBranchAddress("BeamE", &BeamE);
    trees_mctrueacc[0]->GetEntry(0);
    binning_x = plotting::GetBinning(x_observable, BeamE, analysis_id);
    binning_y = plotting::GetBinning(y_observable, BeamE, analysis_id);

    hists_radcorr.push_back(new TH2D(("Rad Corr MC ACC Model " + std::to_string(i)).c_str(), "", binning_x.size() - 1, &binning_x[0], binning_y.size() - 1, &binning_y[0]));
    hists_true.push_back(new TH2D(("True MC ACC Model " + std::to_string(i)).c_str(), "", binning_x.size() - 1, &binning_x[0], binning_y.size() - 1, &binning_y[0]));

    hists_radcorr[i]->Sumw2();
    hists_true[i]->Sumw2();

    std::vector<TTree *> trees = {trees_mctrueacc[i], trees_mcradcorr[i]};

    for (unsigned int j = 0; j < trees.size(); ++j)
    {
      plotting::SetAnalysisBranch( trees[j] ) ;
      for (int k = 0; k < NEntries; ++k)
      {
        trees[j]->GetEntry(k);
        double content_x = GetObservable(x_observable);
        double content_y = GetObservable(y_observable);
        // Weight is the total weight devided by the number of entries.
        // This ensures that we get the same results even if we run less radiated events.
        double w = TotWeight / (double) InitialNEvents;

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

        // Fill the per Sector  histogram
        if (j == 0)
          hists_true[i]->Fill(content_x, content_y, w);
        else if (j == 1)
          hists_radcorr[i]->Fill(content_x, content_y, w);
      }
    }

    ratios.push_back((TH2D *)hists_true[i]->Clone());
    ratios[i]->Divide(hists_radcorr[i]);
    ratios[i]->SetName(("RadCorrModel_" + std::to_string(i)).c_str());
    StandardFormat(ratios[i], title, kBlack + i + 1, 2 + i, x_observable, y_observable);
    ratios[i]->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
    ratios[i]->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
    ratios[i]->GetZaxis()->SetTitle("Radiation correction");
  }

  // Calculate average correction
  TH2D *ratio = (TH2D *)ratios[0]->Clone();
  ratio->SetName("Acceptance");
  for (unsigned int i = 1; i < rad_files.size(); ++i)
  {
    ratio->Add(ratios[i]);
  }
  ratio->Scale(1. / mc_files.size());
  // ratio->Smooth(3);
  StandardFormat(ratio, title, kBlack, 1, x_observable, y_observable);
  ratio->GetXaxis()->SetTitle(GetAxisLabel(x_observable, 0).c_str());
  ratio->GetYaxis()->SetTitle(GetAxisLabel(y_observable, 0).c_str());
  ratio->GetZaxis()->SetTitle("Radiation Correction");

  std::string output_name = output_file_name + "_acceptance_correction_rad_" + x_observable + "_vs_" + y_observable;
  std::string acc_file = "/AcceptanceFiles/" + output_name;

  std::filesystem::path acceptance_path{(output_location + "/AcceptanceFiles").c_str()};
  if (!std::filesystem::exists(acceptance_path))
    std::filesystem::create_directory(acceptance_path);

  TFile outputFile((output_location + acc_file + ".root").c_str(), "RECREATE");

  TCanvas *c_1 = new TCanvas("c_1", "c_1", 200, 10, 700, 500);
  TPad *pad1 = new TPad("pad1", "", 0, 0, 1, 1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Store total contribution (averaged)
  ratio->Write();

  // Plot it
  ratio->SetMarkerStyle(8);
  ratio->Draw("COLZ");

  if (store_root)
    c_1->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_total.root").c_str());
  c_1->SaveAs((output_location + "/AcceptanceFiles/" + output_name + "_total.pdf").c_str());
  delete c_1;
  return acc_file;
}
