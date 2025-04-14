#include "plotting/PlottingUtils.h"

using namespace e4nu;
using namespace e4nu::plotting;

namespace e4nu
{
  namespace plotting {
    // Defining variables to be read from root file here.
    double TotWeight = -9999, ECal= -9999, BeamE = -9999, Recoq3= -9999, RecoW= -9999;
    double EventWght = 1, AccWght = 1, MottXSecScale = 1;
    double Efl= -9999, pfl= -9999, pfl_theta= -9999, pfl_phi= -9999;
    double proton_mom= -9999, proton_phi= -9999, proton_theta= -9999;
    double pim_mom= -9999, pim_theta= -9999, pim_phi= -9999;
    double pip_mom= -9999, pip_theta= -9999, pip_phi= -9999;
    double HadAlphaT= -9999, HadDeltaPT= -9999, HadDeltaPTx= -9999, HadDeltaPTy= -9999, HadDeltaPhiT= -9999;
    double AlphaT= -9999, DeltaPT= -9999, DeltaPhiT= -9999;
    double RecoXBJK= -9999, RecoEnergyTransfer= -9999, RecoQ2= -9999, HadSystemMass= -9999, RecoQELEnu= -9999;
    double MissingEnergy= -9999, MissingAngle= -9999, MissingMomentum= -9999, MissingTransMomentum= -9999, pfl_T= -9999;
    double InferedNucleonMom= -9999, HadronsAngle= -9999, Angleqvshad= -9999;
    double AdlerAngleThetaP= -9999, AdlerAnglePhiP= -9999, AdlerAngleThetaPi= -9999, AdlerAnglePhiPi= -9999;
    double RecoEvPion= -9999, RecoWPion= -9999, ElectronPT= -9999, PionPT= -9999;
    bool IsBkg= 0;
    int ElectronSector= -9999, resid= -9999, InitialNEvents= 1;
    bool QEL= 0, RES= 0, DIS= 0, MEC= 0;
    double MCNormalization= 0, DataNormalization= 0;
    long NEntries = 0;
    TGraph2D* graph_oscillations = nullptr, * graph_oscillations_1 = nullptr, * graph_oscillations_2 = nullptr, * graph_oscillations_3 = nullptr;
  }
}

void plotting::SetAnalysisBranch( TTree * tree ) {
  NEntries = tree->GetEntries();
  if( NEntries == 0 ){
    std::cout << " ROOT tree is empty. Exiting...!" << std::endl;
    exit(0);
  }

  if(tree->GetBranch("TotWeight")) tree->SetBranchAddress("TotWeight", &TotWeight);
  else {
    // This variable should exist. Use as validation
    std::cout << " ROOT file corrupted. TotWeight Branch missing...!, Exiting"<< std::endl;
    exit(0);
  }
  if(tree->GetBranch("InitialNEvents")) tree->SetBranchAddress("InitialNEvents",&InitialNEvents);
  if(tree->GetBranch("EventWght")) tree->SetBranchAddress("EventWght",&EventWght);
  if(tree->GetBranch("AccWght")) tree->SetBranchAddress("AccWght",&AccWght);
  if(tree->GetBranch("MottXSecScale")) tree->SetBranchAddress("MottXSecScale",&MottXSecScale);
  if(tree->GetBranch("IsBkg")) tree->SetBranchAddress("IsBkg", &IsBkg);
  if(tree->GetBranch("ECal")) tree->SetBranchAddress("ECal", &ECal);
  if(tree->GetBranch("BeamE")) tree->SetBranchAddress("BeamE", &BeamE);
  if(tree->GetBranch("pfl_theta")) tree->SetBranchAddress("pfl_theta", &pfl_theta);
  if(tree->GetBranch("pfl_phi")) tree->SetBranchAddress("pfl_phi", &pfl_phi);
  if(tree->GetBranch("pfl")) tree->SetBranchAddress("pfl", &pfl);
  if(tree->GetBranch("Efl")) tree->SetBranchAddress("Efl", &Efl);
  if(tree->GetBranch("pfl_T")) tree->SetBranchAddress("pfl_T", &pfl_T);
  if(tree->GetBranch("proton_mom")) tree->SetBranchAddress("proton_mom", &proton_mom);
  if(tree->GetBranch("proton_theta")) tree->SetBranchAddress("proton_theta", &proton_theta);
  if(tree->GetBranch("proton_phi")) tree->SetBranchAddress("proton_phi", &proton_phi);
  if(tree->GetBranch("pim_mom")) tree->SetBranchAddress("pim_mom", &pim_mom);
  if(tree->GetBranch("pim_theta")) tree->SetBranchAddress("pim_theta", &pim_theta);
  if(tree->GetBranch("pim_phi")) tree->SetBranchAddress("pim_phi", &pim_phi);
  if(tree->GetBranch("pip_mom")) tree->SetBranchAddress("pip_mom", &pip_mom);
  if(tree->GetBranch("pip_theta")) tree->SetBranchAddress("pip_theta", &pip_theta);
  if(tree->GetBranch("pip_phi")) tree->SetBranchAddress("pip_phi", &pip_phi);
  if(tree->GetBranch("RecoW")) tree->SetBranchAddress("RecoW", &RecoW);
  if(tree->GetBranch("Recoq3")) tree->SetBranchAddress("Recoq3", &Recoq3);
  if(tree->GetBranch("RecoQELEnu")) tree->SetBranchAddress("RecoQELEnu", &RecoQELEnu);
  if(tree->GetBranch("RecoXBJK")) tree->SetBranchAddress("RecoXBJK", &RecoXBJK);
  if(tree->GetBranch("RecoQ2")) tree->SetBranchAddress("RecoQ2", &RecoQ2);
  if(tree->GetBranch("RecoEnergyTransfer")) tree->SetBranchAddress("RecoEnergyTransfer", &RecoEnergyTransfer);
  if(tree->GetBranch("AlphaT")) tree->SetBranchAddress("AlphaT", &AlphaT);
  if(tree->GetBranch("HadAlphaT")) tree->SetBranchAddress("HadAlphaT", &HadAlphaT);
  if(tree->GetBranch("DeltaPT")) tree->SetBranchAddress("DeltaPT", &DeltaPT);
  if(tree->GetBranch("HadDeltaPT")) tree->SetBranchAddress("HadDeltaPT", &HadDeltaPT);
  if(tree->GetBranch("HadDeltaPTx")) tree->SetBranchAddress("HadDeltaPTx", &HadDeltaPTx);
  if(tree->GetBranch("HadDeltaPTy")) tree->SetBranchAddress("HadDeltaPTy", &HadDeltaPTy);
  if(tree->GetBranch("DeltaPhiT")) tree->SetBranchAddress("DeltaPhiT", &DeltaPhiT);
  if(tree->GetBranch("HadDeltaPhiT")) tree->SetBranchAddress("HadDeltaPhiT", &HadDeltaPhiT);
  if(tree->GetBranch("ElectronSector")) tree->SetBranchAddress("ElectronSector", &ElectronSector);
  if(tree->GetBranch("HadSystemMass")) tree->SetBranchAddress("HadSystemMass", &HadSystemMass);
  if(tree->GetBranch("MissingEnergy")) tree->SetBranchAddress("MissingEnergy", &MissingEnergy);
  if(tree->GetBranch("MissingTransMomentum")) tree->SetBranchAddress("MissingTransMomentum", &MissingTransMomentum);
  if(tree->GetBranch("MissingAngle")) tree->SetBranchAddress("MissingAngle", &MissingAngle);
  if(tree->GetBranch("MissingMomentum")) tree->SetBranchAddress("MissingMomentum", &MissingMomentum);
  if(tree->GetBranch("InferedNucleonMom")) tree->SetBranchAddress("InferedNucleonMom", &InferedNucleonMom);
  if(tree->GetBranch("HadronsAngle")) tree->SetBranchAddress("HadronsAngle", &HadronsAngle);
  if(tree->GetBranch("AdlerAngleThetaP")) tree->SetBranchAddress("AdlerAngleThetaP", &AdlerAngleThetaP);
  if(tree->GetBranch("AdlerAnglePhiP")) tree->SetBranchAddress("AdlerAnglePhiP", &AdlerAnglePhiP);
  if(tree->GetBranch("AdlerAngleThetaPi")) tree->SetBranchAddress("AdlerAngleThetaPi", &AdlerAngleThetaPi);
  if(tree->GetBranch("AdlerAnglePhiPi")) tree->SetBranchAddress("AdlerAnglePhiPi", &AdlerAnglePhiPi);
  if(tree->GetBranch("Angleqvshad")) tree->SetBranchAddress("Angleqvshad", &Angleqvshad);
  if(tree->GetBranch("RecoEvPion")) tree->SetBranchAddress("RecoEvPion", &RecoEvPion);
  if(tree->GetBranch("RecoWPion")) tree->SetBranchAddress("RecoWPion", &RecoWPion);
  if(tree->GetBranch("ElectronPT")) tree->SetBranchAddress("ElectronPT", &ElectronPT);
  if(tree->GetBranch("PionPT")) tree->SetBranchAddress("PionPT", &PionPT);
  if(tree->GetBranch("QEL")) tree->SetBranchAddress("QEL", &QEL);
  if(tree->GetBranch("RES")) tree->SetBranchAddress("RES", &RES);
  if(tree->GetBranch("MEC")) tree->SetBranchAddress("MEC", &MEC);
  if(tree->GetBranch("DIS")) tree->SetBranchAddress("DIS", &DIS);
  if(tree->GetBranch("resid")) tree->SetBranchAddress("resid", &resid);
  if(tree->GetBranch("MCNormalization")) tree->SetBranchAddress("MCNormalization", &MCNormalization);
  if(tree->GetBranch("DataNormalization")) tree->SetBranchAddress("DataNormalization", &DataNormalization);
}


double plotting::GetObservable(const std::string observable)
{

  double content = 0;
  if (observable == "ECal")
    content = ECal;
  else if (observable == "BeamE")
      content = BeamE;
  else if (observable == "Efl")
      content = Efl;
  else if (observable == "pfl")
    content = pfl;
  else if (observable == "pfl_theta")
    content = pfl_theta;
  else if (observable == "pfl_phi")
    content = pfl_phi;
  else if (observable == "pfl_T")
    content = pfl_T;
  else if (observable == "proton_mom")
    content = proton_mom;
  else if (observable == "proton_theta")
    content = proton_theta;
  else if (observable == "proton_phi")
    content = proton_phi;
  else if (observable == "pim_mom")
    content = pim_mom;
  else if (observable == "pim_theta")
    content = pim_theta;
  else if (observable == "pim_phi")
    content = pim_phi;
  else if (observable == "pip_mom")
    content = pip_mom;
  else if (observable == "pip_theta")
    content = pip_theta;
  else if (observable == "pip_phi")
    content = pip_phi;
  else if (observable == "RecoW")
    content = RecoW;
  else if (observable == "Recoq3")
    content = Recoq3;
  else if (observable == "RecoQELEnu")
    content = RecoQELEnu;
  else if (observable == "RecoXBJK")
    content = RecoXBJK;
  else if (observable == "RecoQ2")
    content = RecoQ2;
  else if (observable == "RecoEnergyTransfer")
    content = RecoEnergyTransfer;
  else if (observable == "AlphaT")
    content = AlphaT;
  else if (observable == "HadAlphaT")
    content = HadAlphaT;
  else if (observable == "DeltaPT")
    content = DeltaPT;
  else if (observable == "HadDeltaPT")
    content = HadDeltaPT;
  else if (observable == "HadDeltaPTx")
    content = HadDeltaPTx;
  else if (observable == "HadDeltaPTy")
    content = HadDeltaPTy;
  else if (observable == "DeltaPhiT")
    content = DeltaPhiT;
  else if (observable == "HadDeltaPhiT")
    content = HadDeltaPhiT;
  else if (observable == "HadSystemMass")
    content = HadSystemMass;
  else if (observable == "MissingEnergy")
    content = MissingEnergy;
  else if (observable == "MissingTransMomentum")
    content = MissingTransMomentum;
  else if (observable == "CorrMissingEnergy")
    content = ComputeMissingEnergy( Efl, HadSystemMass );
  else if (observable == "CorrMissingEnergy1")
    content = ComputeMissingEnergy( Efl, HadSystemMass, 1 );
  else if (observable == "CorrMissingEnergy2")
    content = ComputeMissingEnergy( Efl, HadSystemMass, 2 );
  else if (observable == "CorrMissingEnergy3")
    content = ComputeMissingEnergy( Efl, HadSystemMass, 3 );
  else if (observable == "MissingAngle")
    content = MissingAngle;
  else if (observable == "MissingMomentum")
    content = MissingMomentum;
  else if (observable == "InferedNucleonMom")
    content = InferedNucleonMom;
  else if (observable == "HadronsAngle")
    content = HadronsAngle;
  else if (observable == "AdlerAngleThetaP")
    content = AdlerAngleThetaP;
  else if (observable == "AdlerAnglePhiP")
    content = AdlerAnglePhiP;
  else if (observable == "AdlerAngleThetaPi")
    content = AdlerAngleThetaPi;
  else if (observable == "AdlerAnglePhiPi")
    content = AdlerAnglePhiPi;
  else if (observable == "Angleqvshad")
    content = Angleqvshad;
  else if (observable == "RecoEvPion")
    content = RecoEvPion;
  else if (observable == "RecoWPion")
    content = RecoWPion;
  else if (observable == "ElectronPT")
    content = ElectronPT;
  else if (observable == "PionPT")
    content = PionPT;

  return content;
}

void plotting::NormalizeHist(TH1D *h, double normalization_factor)
{
  // Data normalization
  h->Scale(normalization_factor);
  // h->Sumw2(kFALSE);
  double NBins = h->GetNbinsX();

  for (int i = 1; i <= NBins; i++)
    {
      double content = h->GetBinContent(i);
      double error = h->GetBinError(i);
      double width = h->GetBinWidth(i);
      double newcontent = content / width;
      double newerror = error / width;
      h->SetBinContent(i, newcontent);
      h->SetBinError(i, newerror);
    }
}

void plotting::NormalizeHist(TH2D *h, double normalization_factor)
{
  // Data normalization
  h->Scale(normalization_factor);

  int NbinsX = h->GetNbinsX();
  int NbinsY = h->GetNbinsY();

  for (int i = 1; i <= NbinsX; i++) // Loop over X-axis bins
  {
    for (int j = 1; j <= NbinsY; j++) // Loop over Y-axis bins
    {
      double content = h->GetBinContent(i, j);
      double error = h->GetBinError(i, j);
      double width_x = h->GetXaxis()->GetBinWidth(i);
      double width_y = h->GetYaxis()->GetBinWidth(j);

      double newcontent = content / width_x ;/// width_y;
      double newerror = error / width_x ;/// width_y;

      h->SetBinContent(i, j, newcontent);
      h->SetBinError(i, j, newerror);
    }
  }
}

void plotting::CorrectData(TH1D *h, TH1D *acc)
{
    int NBins = h->GetNbinsX();  // Get the number of bins along the X-axis
    // Loop through each bin and multiply 'h' by the corresponding bin value from 'acc'
    for (int i = 1; i <= NBins; ++i) {
        // Get the value and error from histogram 'h' and 'acc' for bin i
        double h_value = h->GetBinContent(i);
        double h_error = h->GetBinError(i);
        double acc_value = acc->GetBinContent(i);
        double acc_error = acc->GetBinError(i);

        // Multiply the bin content of 'h' by the bin content of 'acc'
        double new_value = h_value * acc_value;
        h->SetBinContent(i, new_value);

        // Calculate the propagated error
        double new_error = std::sqrt( pow(acc_value * h_error,2) + pow(h_value * acc_error,2) );
        h->SetBinError(i, new_error);
    }
}

void plotting::CorrectData(TH2D *h, TH2D *acc)
{
    int NBinsX = h->GetNbinsX(); // Get number of bins along X-axis
    int NBinsY = h->GetNbinsY(); // Get number of bins along Y-axis

    // Loop through each bin in both X and Y directions
    for (int i = 1; i <= NBinsX; ++i) {
        for (int j = 1; j <= NBinsY; ++j) {
            // Get values and errors from histogram 'h' and 'acc' for bin (i, j)
            double h_value = h->GetBinContent(i, j);
            double h_error = h->GetBinError(i, j);
            double acc_value = acc->GetBinContent(i, j);
            double acc_error = acc->GetBinError(i, j);
            // Multiply the bin content of 'h' by the bin content of 'acc'
            double new_value = h_value * acc_value;
            h->SetBinContent(i, j, new_value);

            // Calculate the propagated error
            double new_error = std::sqrt(pow(acc_value * h_error, 2) + pow(h_value * acc_error, 2));
            h->SetBinError(i, j, new_error);
        }
    }
}

std::string plotting::GetAxisLabel(std::string observable, unsigned int id_axis)
{
  std::string x_axis, y_axis;
  if (observable == "ECal")
    {
      x_axis = "E_{Cal} [GeV]";
      y_axis = "d#sigma/dE_{Cal} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "pfl_theta")
    {
      x_axis = "#theta_{e'} [deg]";
      y_axis = "d#sigma/d#theta_{e'} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "pfl_phi")
    {
      x_axis = "#phi_{e'} [deg]";
      y_axis = "d#sigma/d#phi_{e'} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "pfl")
    {
      x_axis = "p_{e'} [GeV/c]";
      y_axis = "d#sigma/dp_{e'} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "proton_mom")
    {
      x_axis = "p_{p} [GeV/c]";
      y_axis = "d#sigma/dp_{p} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "proton_theta")
    {
      x_axis = "#theta_{p} [deg]";
      y_axis = "d#sigma/d#theta_{p} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "proton_phi")
    {
      x_axis = "E_{Cal} [GeV]";
      y_axis = "d#sigma/dE_{Cal} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "pim_mom")
    {
      x_axis = "p_{#pi^{-}} [GeV/c]";
      y_axis = "d#sigma/dp_{#pi^{-}} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "pim_theta")
    {
      x_axis = "#theta_{#pi^{-}} [deg]";
      y_axis = "d#sigma/d#theta_{#pi^{-}} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "pip_mom")
    {
      x_axis = "p_{#pi^{+}} [GeV/c]";
      y_axis = "d#sigma/dp_{#pi^{+}} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "pip_theta")
    {
      x_axis = "#theta_{#pi^{+}} [deg]";
      y_axis = "d#sigma/d#theta_{#pi^{+}} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "RecoW")
    {
      x_axis = "W [GeV]";
      y_axis = "d#sigma/dW #left[#mub GeV^{-1}#right#right]";
    }
  else if (observable == "RecoQELEnu")
    {
      x_axis = "E^{QE} [GeV]";
      y_axis = "d#sigma/dE^{QE} #left[#mub GeV^{-1}#right#right]";
    }
  else if (observable == "RecoXBJK")
    {
      x_axis = "x_{BJK} [GeV]";
      y_axis = "d#sigma/dx_{BJK} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "RecoQ2")
    {
      x_axis = "Q^{2} [GeV^{2}]";
      y_axis = "d#sigma/dQ^{2} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "Recoq3")
    {
      x_axis = "q_{3} [GeV]";
      y_axis = "d#sigma/dq_{3} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "DeltaPT")
    {
      x_axis = "#deltap_{T} [GeV]";
      y_axis = "d#sigma/d#deltap_{T} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "HadDeltaPT")
    {
      x_axis = "#deltap_{T} [GeV]";
      y_axis = "d#sigma/d#deltap_{T} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "HadDeltaPTx")
    {
      x_axis = "#deltap_{Tx} [GeV]";
      y_axis = "d#sigma/d#deltap_{Tx} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "HadDeltaPTy")
    {
      x_axis = "#deltap_{Ty} [GeV]";
      y_axis = "d#sigma/d#deltap_{Ty} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "InferedNucleonMom")
    {
      x_axis = "p_{N,proxy} [GeV]";
      y_axis = "d#sigma/dp_{N,proxy} #left[#mub #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "DeltaPhiT")
    {
      x_axis = "#delta#phi_{T} [deg]";
      y_axis = "d#sigma/d#delta#phi_{T} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "HadDeltaPhiT")
    {
      x_axis = "#delta#phi_{T}^{had} [deg]";
      y_axis = "d#sigma/d#delta#phi_{T}^{had} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "AlphaT")
    {
      x_axis = "#alpha_{T} [deg]";
      y_axis = "d#sigma/d#alpha_{T} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "HadAlphaT")
    {
      x_axis = "#alpha_{T}^{had} [deg]";
      y_axis = "d#sigma/d#alpha_{T}^{had} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "RecoEnergyTransfer")
    {
      x_axis = "#omega [GeV]";
      y_axis = "d#sigma/d#omega #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "HadSystemMass")
    {
      x_axis = "M_{R}[GeV]";
      y_axis = "d#sigma/dM_{R} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "MissingEnergy")
    {
      x_axis = "E_{miss}[GeV]";
      y_axis = "d#sigma/dE_{miss} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "MissingTransMomentum")
    {
      x_axis = "E_{miss}[GeV]";
      y_axis = "d#sigma/dp_{miss}^{T} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "CorrMissingEnergy")
    {
      x_axis = "E_{miss}[GeV]";
      y_axis = "d#sigma/dE_{miss}^{corr} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "CorrMissingEnergy1")
    {
      x_axis = "E_{miss}[GeV]";
      y_axis = "d#sigma/dE_{miss}^{corr,1} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "CorrMissingEnergy2")
    {
      x_axis = "E_{miss}[GeV]";
      y_axis = "d#sigma/dE_{miss}^{corr,2} #left[#mub GeV^{-1}#right]";
    }
  else if (observable == "CorrMissingEnergy3")
    {
      x_axis = "E_{miss}[GeV]";
    y_axis = "d#sigma/dE_{miss}^{corr,3} #left[#mub GeV^{-1}#right]";
  }
  else if (observable == "MissingAngle")
    {
      x_axis = "#theta_{miss}[deg]";
      y_axis = "d#sigma/d#theta_{miss} #left[#mub deg^{-1}#right]";
    }
  else if (observable == "MissingMomentum")
    {
      x_axis = "p_{miss}[GeV/c]";
      y_axis = "d#sigma/dp_{miss} #left[#mub (GeV/c)^{-1}#right]";
    }
  else if (observable == "HadronsAngle")
    {
      x_axis = "#theta_{had}[deg]";
      y_axis = "d#sigma/d#theta_{had} #left[#mub (deg)^{-1}#right]";
    }
  else if (observable == "AdlerAngleThetaP")
    {
      x_axis = "#theta_{p}^{*}[deg]";
      y_axis = "d#sigma/d#theta_{p}^{*} #left[#mub (deg)^{-1}#right]";
    }
  else if (observable == "AdlerAnglePhiP")
    {
      x_axis = "#phi_{p}^{*}[deg]";
      y_axis = "d#sigma/d#phi_{p}^{*} #left[#mub (deg)^{-1}#right]";
    }
  else if (observable == "AdlerAngleThetaPi")
    {
      x_axis = "#theta_{#pi}^{*}[deg]";
      y_axis = "d#sigma/d#theta_{#pi}^{*} #left[#mub (deg)^{-1}#right]";
    }
  else if (observable == "AdlerAnglePhiPi")
    {
      x_axis = "#phi_{#pi}^{*}[deg]";
      y_axis = "d#sigma/d#phi_{#pi}^{*} #left[#mub (deg)^{-1}#right]";
    }
  else if (observable == "Angleqvshad")
    {
      x_axis = "#theta_{#vec{q}#dot#vec{p}_{had}}[deg]";
      y_axis = "d#sigma/d#theta_{#vec{q}#dot#vec{p}_{had}} #left[#mub (deg)^{-1}#right]";
    }
  else if (observable == "RecoEvPion")
    {
      x_axis = "E_{rec} [GeV]";
      y_axis = "d#sigma/dE_{rec} #left[#mub GeV^{-1}#right#right]";
    }
  else if (observable == "RecoWPion")
    {
      x_axis = "W_{rec} [GeV]";
      y_axis = "d#sigma/dW_{rec} #left[#mub GeV^{-1}#right#right]";
    }
  else if (observable == "ElectronPT")
    {
      x_axis = "p_{e',T} [GeV]";
      y_axis = "d#sigma/dp_{e'T} #left[#mub GeV^{-1}#right#right]";
    }
  else if (observable == "PionPT")
    {
      x_axis = "p_{#pi,T} [GeV]";
      y_axis = "d#sigma/dp_{#pi,T} #left[#mub GeV^{-1}#right#right]";
    }

  if (id_axis == 0)
    return x_axis;
  return y_axis;
}

std::string plotting::GetAxisLabel(std::string observable_x, std::string observable_y, unsigned int id_axis) {

  return "Cross Section";
}


std::vector<double> plotting::GetUniformBinning(unsigned int nbins, double min, double max)
{

  std::vector<double> binning;
  double step = (max - min) / nbins;
  for (unsigned int i = 0; i < nbins + 1; ++i)
    {
      binning.push_back(min + i * step);
    }

  return binning;
}

std::vector<double> plotting::GetECalBinning(unsigned int nbins_tail, unsigned int nbins_peak, double min, double max, double EBeam)
{
  std::vector<double> binning;
  double temp_min = min;
  double temp_max = max;
  if (EBeam < max)
    temp_max = EBeam * (1 - 0.05);
  double step = (temp_max - temp_min) / nbins_tail;
  for (unsigned int i = 0; i < nbins_tail + 1; ++i)
    {
      binning.push_back(temp_min + i * step);
    }

  step = (max - temp_max) / nbins_peak;
  for (unsigned int i = 1; i < nbins_peak + 1; ++i)
    {
      binning.push_back(temp_max + i * step);
    }
  return binning;
}

std::vector<double> plotting::GetBinning(std::string observable, double EBeam, std::string analysis_key)
{
  std::vector<double> binning;

  if (observable == "ECal")
    {
      if (EBeam == 1.161)
	binning = plotting::GetECalBinning(13, 15, 0.7, EBeam + 0.15, EBeam);
      else if (EBeam == 2.261)
	binning = plotting::GetECalBinning(13, 15, 1, EBeam + 0.15, EBeam);
      else if (EBeam == 4.461)
	binning = plotting::GetECalBinning(13, 15, 1.5, EBeam + 0.15, EBeam);
    }
  else if (observable == "DiffECal")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(25, -0.6, 0.2);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(25, -0.6, 0.2);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(25, -0.6, 0.2);
    }
  else if (observable == "pfl_theta")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(25, 24, 45);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(25, 24, 45);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(25, 15, 45);
    }
  else if (observable == "pfl_phi")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 180);
    }
  else if (observable == "pfl")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0.35, 0.9);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0.5, 1.7);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 1.2, 3.8);
    }
  else if (observable == "proton_mom")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0.3, 1.1);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0.3, 2);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0.3, 3);
    }
  else if (observable == "proton_theta")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 12, 110);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 12, 110);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 12, 110);
    }
  else if (observable == "proton_phi")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(35, 0, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(35, 0, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(35, 0, 180);
    }
  else if (observable == "pim_mom")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(25, 0.15, 0.6);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(25, 0.15, 1.6);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(25, 0.15, 1.6);
    }
  else if (observable == "pim_theta")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(35, 12, 140);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(35, 12, 140);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(35, 12, 140);
    }
  else if (observable == "pip_mom")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(35, 0.15, 3);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(35, 0.15, 3);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(35, 0.15, 3);
    }
  else if (observable == "pip_theta")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 12, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 12, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 12, 180);
    }
  else if (observable == "RecoW")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 1, 1.5);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 1, 2);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 1, 2.4);
    }
  else if (observable == "RecoXBJK")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 0.9);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 0.9);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 1);
    }
  else if (observable == "RecoQ2")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0.15, 0.45);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0.3, 1.5);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0.9, 3);
    }
  else if (observable == "RecoQELEnu")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0., 1.3);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0.3, EBeam + 0.2);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0.9, EBeam + 0.2);
    }
  else if (observable == "Recoq3")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, EBeam + 0.2);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, EBeam + 0.2);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, EBeam + 0.2);
    }
  else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 1);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 1);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 1);
    }
  else if (observable == "HadDeltaPTx")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, -0.6, 0.6);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, -1, 1);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, -1, 1);
    }
  else if (observable == "HadDeltaPTy")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, -0.6, 0.6);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, -1, 1);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, -1, 1);
    }
  else if (observable == "HadDeltaPhiT" || observable == "DeltaPhiT")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 80);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 80);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 80);
    }
  else if (observable == "AlphaT")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 180);
    }
  else if (observable == "HadAlphaT")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 180);
    }
  else if (observable == "RecoEnergyTransfer")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0.3, 0.8);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0.5, 2);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 4);
    }
  else if (observable == "HadSystemMass")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 1, 1.6);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 1, 2);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 2.7);
    }
  else if (observable == "MissingEnergy")
    {
      if (EBeam == 1.161)
	     binning = plotting::GetECalBinning(20, 15, 0.3, 1.1, 0.9);
      else if (EBeam == 2.261)
        binning = plotting::GetECalBinning(20, 15, -0.7, 1.2, 0.9);
      else if (EBeam == 4.461)
	     binning = plotting::GetECalBinning(20, 15, -2.5, 1.2, 0.9);
     }
  else if (observable == "MissingTransMomentum")
     {
       if (EBeam == 1.161)
	      binning = plotting::GetECalBinning(20, 15, 0.3, 1.1, 0.9);
       else if (EBeam == 2.261)
         binning = plotting::GetECalBinning(20, 15, -0.7, 1.2, 0.9);
       else if (EBeam == 4.461)
   	    binning = plotting::GetECalBinning(20, 15, -2.5, 1.2, 0.9);
      }
  else if (observable == "CorrMissingEnergy" || observable == "CorrMissingEnergy1" || observable == "CorrMissingEnergy2" || observable == "CorrMissingEnergy3")
    {
      if (EBeam == 1.161)
        binning = plotting::GetECalBinning(20, 15, 0.3, 1.1, 0.9);
      else if (EBeam == 2.261)
        binning = plotting::GetECalBinning(20, 15, -0.7, 1.2, 0.9);
      else if (EBeam == 4.461)
        binning = plotting::GetECalBinning(20, 15, -2.5, 1.2, 0.9);
      }
  else if (observable == "MissingAngle")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 15, 180);
    }
  else if (observable == "MissingMomentum")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 1);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 2);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 1.5);
    }
  else if (observable == "InferedNucleonMom")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(30, 0, 0.8);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(30, 0, 1);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(30, 0, 1);
  }
  else if (observable == "HadronsAngle")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(30, 20, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(30, 20, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(30, 20, 180);
    }
  else if (observable == "AdlerAngleThetaP")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 180);
    }
  else if (observable == "AdlerAnglePhiP")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(20, 0, 180);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(20, 0, 180);
    } else if (observable == "AdlerAngleThetaPi")
{
  if (EBeam == 1.161)
    binning = plotting::GetUniformBinning(20, 0, 180);
  else if (EBeam == 2.261)
    binning = plotting::GetUniformBinning(20, 0, 180);
  else if (EBeam == 4.461)
    binning = plotting::GetUniformBinning(10, 0, 180);
} else if (observable == "AdlerAnglePhiPi")
  {
    if (EBeam == 1.161)
    binning = plotting::GetUniformBinning(30, 20, 180);
    else if (EBeam == 2.261)
    binning = plotting::GetUniformBinning(20, 20, 180);
    else if (EBeam == 4.461)
    binning = plotting::GetUniformBinning(20, 20, 180);
} else if (observable == "Angleqvshad")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(25, 0, 120);
      else if (EBeam == 2.261)
	binning = plotting::GetUniformBinning(25, 0, 120);
      else if (EBeam == 4.461)
	binning = plotting::GetUniformBinning(25, 0, 60);
    }
  else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam == 1.161)
	binning = plotting::GetUniformBinning(20, 0, 0.7);
    }
  else if (analysis_key == "1pip")
    {
      if (observable == "ECal")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 20, 180);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 20, 180);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 20, 180);
	}
      else if (observable == "AdlerAnglePhiPi")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(30, 20, 180);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 20, 180);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 20, 180);
	}
      else if (observable == "Angleqvshad")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0, 70);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 0, 50);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0, 30);
	}
      else if (observable == "RecoEvPion")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0, 2);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 0.5, 3.5);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 2, 6);
	}
      else if (observable == "RecoWPion")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0.5, 2);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 1, 2);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0.5, 3.5);
	}
      else if (observable == "ElectronPT")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0.2, 0.7);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 0, 1);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0, 1);
	}
      else if (observable == "PionPT")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0, 0.5);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 0, 1);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0, 1);
	}
      else if (observable == "Angleqvshad")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0, 120);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 0, 120);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0, 60);
	}
      else if (observable == "HadDeltaPT" || observable == "DeltaPT")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0, 0.7);
	}
    }

  if (analysis_key == "1p1pim")
    return binning;
  else if (analysis_key == "1p1pip")
    {
      if (observable == "ECal")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetECalBinning(15, 10, 0.8, EBeam + 0.2, EBeam);
	  else if (EBeam == 2.261)
	    binning = plotting::GetECalBinning(15, 10, 1.2, EBeam + 0.2, EBeam);
	  else if (EBeam == 4.461)
	    binning = plotting::GetECalBinning(15, 10, 2, EBeam + 0.2, EBeam);
	}
      else if (observable == "proton_mom")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(25, 0.3, 1);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(30, 0.3, 2);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(30, 0.3, 3);
	}
  else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam == 1.161)
  binning = plotting::GetUniformBinning(25, 0, 1);
      else if (EBeam == 2.261)
  binning = plotting::GetUniformBinning(25, 0, 1);
      else if (EBeam == 4.461)
  binning = plotting::GetUniformBinning(25, 0, 1);
    }
      else if (observable == "pip_mom")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(25, 0.2, 0.6);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(25, 0.3, 1.2);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(25, 0.3, 2.2);
	}
      else if (observable == "pip_theta")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 12, 120);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(25, 12, 130);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(30, 12, 130);
	}
      else if (observable == "MissingEnergy")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(25, 0.5, 1);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(25, 0, 1);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(25, 0, 1);
	}
  else if (observable == "CorrMissingEnergy"||observable == "CorrMissingEnergy1"||observable == "CorrMissingEnergy2"||observable == "CorrMissingEnergy3")
{
if (EBeam == 1.161)
  binning = plotting::GetUniformBinning(25, 0.5, 1);
else if (EBeam == 2.261)
  binning = plotting::GetUniformBinning(25, 0, 1);
else if (EBeam == 4.461)
  binning = plotting::GetUniformBinning(25, 0, 1);
}
      else if (observable == "MissingAngle")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 35, 180);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 35, 180);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(10, 35, 180);
	}
      else if (observable == "MissingMomentum")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(30, 0, 1);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(30, 0, 1.8);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(30, 0, 1.5);
	}
      else if (observable == "InferedNucleonMom")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(30, 0, 0.8);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(30, 0, 1);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0, 1);
	}
      else if (observable == "Angleqvshad")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(30, 0, 120);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(30, 0, 120);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0, 120);
	}
      else if (observable == "HadronsAngle")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(30, 0, 180);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(30, 0, 180);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0, 180);
	}
      else if (observable == "AdlerAngleThetaPi")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0, 180);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 0, 180);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(10, 0, 180);
	}
      else if (observable == "AdlerAnglePhiPi")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0, 180);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(15, 0, 180);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(10, 0, 180);
	}
      else if (observable == "Angleqvshad")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(25, 0, 120);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(25, 0, 120);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(25, 0, 60);
	}
      else if (observable == "RecoW")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 1.1, 1.45);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 1, 1.9);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 1.2, 2.3);
	}
    }
  else if (analysis_key == "1pim")
    {
      if (observable == "ECal")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetECalBinning(20, 10, 0.6, EBeam + 0.2, EBeam);
	  else if (EBeam == 2.261)
	    binning = plotting::GetECalBinning(20, 10, 0.6, EBeam + 0.2, EBeam);
	  else if (EBeam == 4.461)
	    binning = plotting::GetECalBinning(15, 10, 1.2, EBeam + 0.2, EBeam);
	}
      if (observable == "RecoEvPion")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(50, 0.5, 2);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(50, 0.5, 3.5);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(50, 1, 7);
	}
      if (observable == "RecoWPion")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0.9, 2);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 1, 2);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0, 4);
	}
      else if (observable == "ElectronPT")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(50, 0, 0.5);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(50, 0, 1);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(50, 0, 1);
	}
      else if (observable == "PionPT")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(50, 0, 0.5);
	}
      else if (observable == "RecoEvPion")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(30, 0.5, 2);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(30, 0.5, 3.5);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(30, 2, 6);
	}
      else if (observable == "RecoWPion")
	{
	  if (EBeam == 1.161)
	    binning = plotting::GetUniformBinning(20, 0.5, 2);
	  else if (EBeam == 2.261)
	    binning = plotting::GetUniformBinning(20, 1, 2);
	  else if (EBeam == 4.461)
	    binning = plotting::GetUniformBinning(20, 0.5, 3.5);
	}
    }
  else if (analysis_key == "1pip") {

    if (observable == "ECal")
      {
	if (EBeam == 1.161)
	  binning = plotting::GetECalBinning(15, 10, 0.6, EBeam + 0.2, EBeam);
	else if (EBeam == 2.261)
	  binning = plotting::GetECalBinning(15, 10, 1, EBeam + 0.2, EBeam);
	else if (EBeam == 4.461)
	  binning = plotting::GetECalBinning(15, 10, 1, EBeam + 0.2, EBeam);
      }
    else if (observable == "pip_mom")
      {
	if (EBeam == 1.161)
	  binning = plotting::GetUniformBinning(20, 0, 0.6);
	else if (EBeam == 2.261)
	  binning = plotting::GetUniformBinning(20, 0.3, 1.2);
	else if (EBeam == 4.461)
	  binning = plotting::GetUniformBinning(20, 0, 2.5);
      }
    else if (observable == "pip_theta")
      {
	if (EBeam == 1.161)
	  binning = plotting::GetUniformBinning(20, 0, 130);
	else if (EBeam == 2.261)
	  binning = plotting::GetUniformBinning(20, 0, 120);
	else if (EBeam == 4.461)
	  binning = plotting::GetUniformBinning(20, 0, 100);
      }
    if (observable == "RecoEvPion")
      {
	if (EBeam == 1.161)
	  binning = plotting::GetUniformBinning(50, 0.5, 2);
	else if (EBeam == 2.261)
	  binning = plotting::GetUniformBinning(50, 0.5, 3.5);
	else if (EBeam == 4.461)
	  binning = plotting::GetUniformBinning(50, 1, 7);
      }
    if (observable == "RecoWPion")
      {
	if (EBeam == 1.161)
	  binning = plotting::GetUniformBinning(50, 0.9, 2);
	else if (EBeam == 2.261)
	  binning = plotting::GetUniformBinning(50, 1, 2);
	else if (EBeam == 4.461)
	  binning = plotting::GetUniformBinning(50, 2, 4);
      }
    else if (observable == "ElectronPT")
      {
	if (EBeam == 1.161)
	  binning = plotting::GetUniformBinning(50, 0, 0.5);
	else if (EBeam == 2.261)
	  binning = plotting::GetUniformBinning(50, 0, 1);
	else if (EBeam == 4.461)
	  binning = plotting::GetUniformBinning(50, 0, 1);
      }
    else if (observable == "PionPT")
      {
	if (EBeam == 1.161)
	  binning = plotting::GetUniformBinning(50, 0, 0.5);
	else if (EBeam == 2.261)
	  binning = plotting::GetUniformBinning(50, 0, 1);
	else if (EBeam == 4.461)
	  binning = plotting::GetUniformBinning(50, 0, 1);
      }
  }


  if (binning.size() == 0)
    {
      std::cout << " ERROR: Binning for " << observable << " is nul" << std::endl;
    }
  return binning;
}

std::vector<double> plotting::GetAdditionalBinning(std::string second_observable, double EBeam, std::string analysis_id)
{
  // In some cases, we might want to add additional plots.
  // In particular, we might want to break the plot into additonal plots as a function of a second observable
  // This function returns the binning for this second observable.
  // The "binning" corresponds to the ranges of interest
  std::vector<double> binning;
  std::vector<double> original_binning = plotting::GetBinning(second_observable, EBeam, analysis_id);
  if (second_observable == "ECal")
    {
      binning.push_back(original_binning[0]);
      binning.push_back(EBeam * (1 - 0.05));
      binning.push_back(original_binning[original_binning.size() - 1]);
    }
  else if (second_observable == "HadDeltaPT" || second_observable == "DeltaPT")
    {
      binning.push_back(original_binning[0]);
      binning.push_back(0.2);
      // binning.push_back(0.4);
      binning.push_back(original_binning[original_binning.size() - 1]);
    }
  else if (second_observable == "HadAlphaT" || second_observable == "AlphaT")
    {
      binning.push_back(original_binning[0]);
      binning.push_back(45);
      binning.push_back(original_binning[original_binning.size() - 1]);
    }

  return binning;
}

std::string plotting::GetAlternativeObs(std::string observable)
{
  if (observable == "ECal")
    return "HadDeltaPT";
  else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    return "ECal";
  else if (observable == "HadAlphaT" || observable == "AlphaT")
    return "ECal";
  return "";
}

std::string plotting::GetObsName(std::string observable)
{
  if (observable == "ECal")
    return "E_{Cal}";
  else if (observable == "HadDeltaPT")
    return "#deltap^{had}_{T}";
  else if (observable == "HadAlphaT")
    return "#alpha^{had}_{T}";
  return "";
}

std::string plotting::GetUnit(std::string observable)
{
  if (observable == "ECal")
    return "[GeV]";
  else if (observable == "HadDeltaPT")
    return "[GeV/c]";
  else if (observable == "HadAlphaT")
    return "[deg]";
  return "";
}

double plotting::GetMaximum(std::vector<TH1D *> predictions)
{
  double max = 0;
  for (unsigned int i = 0; i < predictions.size(); ++i) {
      for (int j = 1; j <= predictions[i]->GetNbinsX(); ++j) { // Note the range: 1 to GetNbinsX()
          if (max < predictions[i]->GetBinContent(j)) max = predictions[i]->GetBinContent(j);
      }
  }
  return max * (1 + 0.05); // Add a 5% buffer to the maximum
}

double plotting::GetMaximum(std::vector<TH2D *> predictions)
{
    double max = 0;
    for (unsigned int i = 0; i < predictions.size(); ++i) {
        int nBinsX = predictions[i]->GetNbinsX();
        int nBinsY = predictions[i]->GetNbinsY();
        for (int xBin = 1; xBin <= nBinsX; ++xBin) { // Loop over X bins (1-based index)
            for (int yBin = 1; yBin <= nBinsY; ++yBin) { // Loop over Y bins (1-based index)
                double binContent = predictions[i]->GetBinContent(xBin, yBin);
                if (max < binContent) max = binContent;
            }
        }
    }
    return max * (1 + 0.05); // Add a 5% buffer to the maximum
}

double plotting::GetMinimum(std::vector<TH1D *> predictions)
{
  double min = 9999;
  for (unsigned int i = 0; i < predictions.size(); ++i)
    {
      for (unsigned int j = 0; j < predictions[i]->GetNbinsX(); ++j)
	{
	  if (min < predictions[i]->GetBinContent(j))
	    min = predictions[i]->GetBinContent(j);
	}
    }
  return min * (1 - 0.12);
}

void plotting::StandardFormat(TH1D *prediction, std::string title, int color, int style, std::string observable, bool is_log, double y_max, std::string y_axis_label )
{
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPaperSize(20, 26);
  gStyle->SetTitleFont(132, "pad");
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineStyleString(2, "[12 12]");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  prediction->SetLineColor(color);
  prediction->SetLineStyle(style);
  prediction->SetMarkerStyle(style);
  prediction->SetMarkerColor(color);
  prediction->SetLineWidth(2);

  prediction->SetTitle(title.c_str());
  // prediction -> SetTitleFont(13);
  prediction->GetXaxis()->SetTitle(GetAxisLabel(observable, 0).c_str());

  if (y_axis_label == "")
    prediction->GetYaxis()->SetTitle(GetAxisLabel(observable, 1).c_str());
  else
    prediction->GetYaxis()->SetTitle(y_axis_label.c_str());

  prediction->GetXaxis()->CenterTitle();
  prediction->GetYaxis()->CenterTitle();

  prediction->BufferEmpty(-1);
  if (y_max == 0)
    {
      double max = -999;
      for (unsigned int k = 0; k < prediction->GetNbinsX(); ++k)
	{
	  if (prediction->GetBinContent(k) > max)
	    max = prediction->GetBinContent(k);
	}
      // y_max = (prediction -> GetMaximum()) * ( 1+0.2 );
      y_max = max * (1 + 0.2);
  }

  if( y_max == 0 ) y_max = 100; // for empty plots.

  int FontStyle = 132;
  prediction->GetXaxis()->SetTitleOffset(1);
  prediction->GetXaxis()->SetLabelSize(0.05);
  prediction->GetXaxis()->SetTitleSize(0.08);
  prediction->GetXaxis()->SetNdivisions(6);
  prediction->GetXaxis()->SetLabelFont(FontStyle);
  prediction->GetXaxis()->SetTitleFont(FontStyle);

  prediction->GetYaxis()->SetNdivisions(8);
  prediction->GetYaxis()->SetTitleOffset(0.5);
  prediction->GetYaxis()->SetLabelSize(0.05);
  prediction->GetYaxis()->SetTitleSize(0.08);
  prediction->GetYaxis()->SetLabelFont(43);
  prediction->GetYaxis()->SetLabelFont(FontStyle);
  prediction->GetYaxis()->SetTitleFont(FontStyle);
  if( is_log ) prediction->GetYaxis()->SetRangeUser(1./1E4, y_max);
  else prediction->GetYaxis()->SetRangeUser(0, y_max);
  prediction->GetYaxis()->SetMaxDigits(3);
  prediction->SetTitleFont(FontStyle);

  return;
}

void plotting::StandardFormat(TH2D *prediction, std::string title, int color, int style, std::string observable_x, std::string observable_y, double z_max, std::string z_axis_label)
{
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPaperSize(20, 26);
  gStyle->SetTitleFont(132, "pad");
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineStyleString(2, "[12 12]");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  prediction->SetLineColor(color);
  prediction->SetLineStyle(style);
  prediction->SetMarkerStyle(style);
  prediction->SetMarkerColor(color);
  prediction->SetLineWidth(2);

  prediction->SetTitle(title.c_str());
  // prediction -> SetTitleFont(13);
  prediction->GetXaxis()->SetTitle(GetAxisLabel(observable_x, 0).c_str());
  prediction->GetYaxis()->SetTitle(GetAxisLabel(observable_y, 0).c_str());

  if (z_max == 0)
    {
      double max = -999;
      for (unsigned int k = 0; k < prediction->GetNbinsX(); ++k)
      {
        if (prediction->GetBinContent(k) > max)
        max = prediction->GetBinContent(k);
      }
      // y_max = (prediction -> GetMaximum()) * ( 1+0.2 );
      z_max = max * (1 + 0.2);
    }

  if (z_axis_label == "")
    prediction->GetZaxis()->SetTitle(GetAxisLabel(observable_x,observable_y, 1).c_str());
  else
    prediction->GetZaxis()->SetTitle(z_axis_label.c_str());

  prediction->GetXaxis()->CenterTitle();
  prediction->GetYaxis()->CenterTitle();
  prediction->GetZaxis()->CenterTitle();

  prediction->BufferEmpty(-1);

  int FontStyle = 132;
  prediction->GetXaxis()->SetTitleOffset(1);
  prediction->GetXaxis()->SetLabelSize(0.05);
  prediction->GetXaxis()->SetTitleSize(0.05);
  prediction->GetXaxis()->SetNdivisions(6);
  prediction->GetXaxis()->SetLabelFont(FontStyle);
  prediction->GetXaxis()->SetTitleFont(FontStyle);

  prediction->GetYaxis()->SetNdivisions(8);
  prediction->GetYaxis()->SetTitleOffset(0.9);
  prediction->GetYaxis()->SetLabelSize(0.05);
  prediction->GetYaxis()->SetTitleSize(0.05);
  prediction->GetYaxis()->SetLabelFont(43);
  prediction->GetYaxis()->SetLabelFont(FontStyle);
  prediction->GetYaxis()->SetTitleFont(FontStyle);
  // prediction->GetZaxis()->SetRangeUser(1./1E6, z_max);
  prediction->GetZaxis()->SetMaxDigits(3);
  prediction->SetTitleFont(FontStyle);

  return;
}

int plotting::ColorBlindPalette(int color_id){
  // using color blind palette developed by CMS
  // https://arxiv.org/abs/2107.02270
  const std::vector<std::string> sequence({"#3f90da","#ffa90e","#bd1f01","#94a4a2","#832db6","#a96b59","#e76300","#b9ac70","#717581","#92dadd"});
  return TColor::GetColor( sequence[ color_id%10 ].c_str() );
}


std::vector<std::string> plotting::SplitString(std::string s, char d)
{
  std::vector<std::string> strings;
  int startIndex = 0, endIndex = 0;
  for (unsigned int i = 0; i <= s.size(); i++)
    {

      // If we reached the end of the word or the end of the input.
      if (s[i] == d || i == s.size())
	{
	  endIndex = i;
	  std::string temp;
	  temp.append(s, startIndex, endIndex - startIndex);
	  strings.push_back(temp);
	  startIndex = endIndex + 1;
	}
    }
  return strings;
}

std::string plotting::GetArg(std::string op, int argc, char **argv)
{
  const int buf_size = 2048 * 128;
  char *argument = new char[buf_size];
  strcpy(argument, "");

  while (argc > 2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-')
	{

	  char op_cur[buf_size];
	  strcpy(op_cur, &argv[1][2]);

	  if (strcmp(op.c_str(), op_cur) == 0)
	    {
	      if (strlen(&argv[2][0]))
		{
		  strcpy(argument, &argv[2][0]);
		}
	    }
	}
      argc--;
      argv++;
    }

  std::string value = std::string(argument);
  delete[] argument;
  return value;
}

bool plotting::ExistArg(std::string op, int argc, char **argv)
{
  const int buf_size = 2048 * 128;
  char *argument = new char[buf_size];
  strcpy(argument, "");

  while (argc > 2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-')
	{

	  char op_cur[buf_size];
	  strcpy(op_cur, &argv[1][2]);

	  if (strcmp(op.c_str(), op_cur) == 0)
	    {
	      return true;
	    }
	}
      argc--;
      argv++;
    }
  delete[] argument;
  return false;
}

bool plotting::PlotZoomIn(std::string analysis_id)
{
  if (analysis_id == "1p1pim")
    return true;
  else if (analysis_id == "1p1pip" || analysis_id == "1pim" || analysis_id == "1pip")
    return false;
  return true;
}


void plotting::GetMissingEnergyGraph( const std::string mc_file ){
  // This function coputes a 2D graph for MC only with the following variables:
  // Efl, MissingEnergy and Ehad (HadSystemMass)
  // As done in NOVA, we can attempt to correct e-data with the MC calculation and directly estimate the bias

  // Open MC file
  TFile* in_root_file = new TFile((mc_file).c_str(),"ROOT") ;
  if( !in_root_file ) {
    std::cout << " ERROR: " << mc_file << " does not exist."<<std::endl;
    return ;
  }

  TTree* in_tree = (TTree*)in_root_file->Get("MCCLAS6Tree") ;
  if( !in_tree ) {
    std::cout << " ERROR: MCCLAS6Tree does not exist in "<< mc_file <<std::endl;
    return ;
  }

  plotting::SetAnalysisBranch(in_tree);
  // Retrieve the histogram contents (averages for each bin)
  std::vector<double> x_values, x_values_1, x_values_2, x_values_3;
  std::vector<double> y_values, y_values_1, y_values_2, y_values_3;
  std::vector<double> z_values, z_values_1, z_values_2, z_values_3;
  for( int j = 0 ; j < NEntries ; ++j ) {
    in_tree->GetEntry(j) ;
    double content_x = plotting::GetObservable("Efl");
    double content_y = plotting::GetObservable("HadSystemMass");
    double content_z = plotting::GetObservable("BeamE") - plotting::GetObservable("MissingEnergy");
    double content_pt = plotting::GetObservable("pfl_T");

    x_values.push_back(content_x);
    y_values.push_back(content_y);
    z_values.push_back(content_z);

    // Fill for slices of PT:
    if( content_pt < 0.4  && content_pt > 0 ) {
      x_values_1.push_back(content_x);
      y_values_1.push_back(content_y);
      z_values_1.push_back(content_z);
    } else if ( content_pt < 0.6  && content_pt > 0.4 ) {
      x_values_2.push_back(content_x);
      y_values_2.push_back(content_y);
      z_values_2.push_back(content_z);
    } else if ( content_pt > 0.6 ){
      x_values_3.push_back(content_x);
      y_values_3.push_back(content_y);
      z_values_3.push_back(content_z);
    }
  }

  // Create 2D graph:
  if( x_values.size() != 0) {
    graph_oscillations= new TGraph2D(x_values.size(), &x_values[0], &y_values[0], &z_values[0]);
    graph_oscillations->SetName("TGraph2D_0") ;
  }
  // For slices of PT:
  if( x_values_1.size() != 0) {
    graph_oscillations_1 = new TGraph2D(x_values_1.size(), &x_values_1[0], &y_values_1[0], &z_values_1[0]);
    graph_oscillations_1->SetName("TGraph2D_1") ;
  }
  if( x_values_2.size() != 0) {
    graph_oscillations_2 = new TGraph2D(x_values_2.size(), &x_values_2[0], &y_values_2[0], &z_values_2[0]);
    graph_oscillations_2->SetName("TGraph2D_2") ;
  }
  if( x_values_3.size() != 0) {
    graph_oscillations_3 = new TGraph2D(x_values_3.size(), &x_values_3[0], &y_values_3[0], &z_values_3[0]);
    graph_oscillations_3->SetName("TGraph2D_3") ;
  }

}

double plotting::ComputeMissingEnergy( const double event_efl, const double event_ehad, const unsigned int slice ){
  if( !graph_oscillations  ){
    std::cout << " ERROR: you did not compute graph oscillations from MC file. "<< std::endl;
    return 0;
  }
  if( slice ==0 ) graph_oscillations->Interpolate(event_efl,event_ehad);
  else if( slice == 1 && graph_oscillations_1 ) graph_oscillations_1->Interpolate(event_efl,event_ehad) ; // pt < 0.4
  else if( slice == 2 && graph_oscillations_2 ) graph_oscillations_2->Interpolate(event_efl,event_ehad) ; // 0.4 < pt < 0.6
  else if( slice == 3 && graph_oscillations_3 ) graph_oscillations_3->Interpolate(event_efl,event_ehad) ; // pt > 0.6
  return 0 ;
}
