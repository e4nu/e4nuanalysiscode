#include "plotting/PlottingUtils.h"

using namespace e4nu;

void plotting::NormalizeHist(TH1D *h, double normalization_factor)
{
  // Data normalization
  h->Sumw2(kFALSE);
  h->Scale(normalization_factor);
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

void plotting::CorrectData(TH1D *h, TH1D *acc)
{
  double NBins = h->GetNbinsX();
  h->Sumw2(kFALSE);
  h->Multiply(acc);
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
    x_axis = "Q^{2} [GeV]";
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
      binning = plotting::GetECalBinning(22, 10, 0.8, 1.2, EBeam);
    else if (EBeam == 2.261)
      binning = plotting::GetECalBinning(20, 10, 1, EBeam + 0.06, EBeam);
    else if (EBeam == 4.461)
      binning = plotting::GetECalBinning(20, 10, 1.5, EBeam + 0.1, EBeam);
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
      binning = plotting::GetUniformBinning(35, 20, 45);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(35, 20, 45);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 15, 45);
  }
  else if (observable == "pfl_phi")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(35, 0, 180);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(35, 0, 180);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 0, 180);
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
      binning = plotting::GetUniformBinning(30, 0.3, 1.1);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(30, 0.3, 2);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(30, 0.3, 3);
  }
  else if (observable == "proton_theta")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(35, 12, 110);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(35, 12, 110);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 12, 110);
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
      binning = plotting::GetUniformBinning(50, 0.15, 1.6);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 0.15, 1.6);
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
      binning = plotting::GetUniformBinning(35, 12, 180);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(35, 12, 180);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 12, 180);
  }
  else if (observable == "RecoW")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(35, 1, 1.5);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(35, 1, 1.9);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 1.2, 2.5);
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
      binning = plotting::GetUniformBinning(35, 0, 1);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(35, 0, 1);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 0, 1);
  }
  else if (observable == "HadDeltaPTx")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(20, -1, 1);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(20, -1, 1);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(20, -1, 1);
  }
  else if (observable == "HadDeltaPTy")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(20, -1, 1);
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
      binning = plotting::GetUniformBinning(35, 0, 180);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(35, 0, 180);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 0, 180);
  }
  else if (observable == "HadAlphaT")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(30, 0, 180);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(30, 0, 180);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(30, 0, 180);
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
      binning = plotting::GetUniformBinning(20, 0.5, 1);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(20, 0, 1);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(20, 0, 1);
  }
  else if (observable == "MissingAngle")
  {
    if (EBeam == -1.161)
      binning = plotting::GetUniformBinning(35, 0, 180);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(35, 0, 180);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(35, 15, 180);
  }
  else if (observable == "MissingMomentum")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(25, 0, 1);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(25, 0, 2);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(25, 0, 1.5);
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
      binning = plotting::GetUniformBinning(30, 20, 180);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(30, 20, 180);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(30, 20, 180);
  }
  else if (observable == "AdlerAnglePhiP")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(30, 20, 180);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(30, 20, 180);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(30, 20, 180);
  }
  else if (observable == "AdlerAngleThetaPi")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(30, 20, 180);
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
      binning = plotting::GetUniformBinning(30, 0, 70);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(30, 0, 50);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(30, 0, 30);
  }
  else if (observable == "RecoEvPion")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(50, 0, 2);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(50, 0.5, 3.5);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(50, 2, 6);
  }
  else if (observable == "RecoWPion")
  {
    if (EBeam == 1.161)
      binning = plotting::GetUniformBinning(50, 0.5, 2);
    else if (EBeam == 2.261)
      binning = plotting::GetUniformBinning(50, 0.5, 3);
    else if (EBeam == 4.461)
      binning = plotting::GetUniformBinning(50, 0.5, 3.5);
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

  if (analysis_key == "1p1pim")
    return binning;
  else if (analysis_key == "1p1pip")
  {
    if (observable == "ECal")
    {
      if (EBeam == 1.161)
        binning = plotting::GetECalBinning(15, 10, 0.8, 1.2, EBeam);
      else if (EBeam == 2.261)
        binning = plotting::GetECalBinning(15, 10, 1.2, EBeam + 0.06, EBeam);
      else if (EBeam == 4.461)
        binning = plotting::GetECalBinning(15, 10, 2, EBeam + 0.1, EBeam);
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
    else if (observable == "pip_mom")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(20, 0.15, 0.6);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(35, 0.3, 1.2);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(35, 0.3, 2.2);
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
    else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(20, 0, 0.8);
    }
    else if (observable == "MissingEnergy")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(20, 0.5, 1);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(20, 0, 1);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(10, 0, 1);
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
        binning = plotting::GetUniformBinning(25, 0, 1);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(25, 0, 1.8);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(25, 0, 1.5);
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
        binning = plotting::GetUniformBinning(30, 0, 120);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(30, 0, 120);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(20, 0, 60);
    }
    else if (observable == "RecoW")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(25, 1.1, 1.45);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(35, 1, 1.9);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(35, 1.2, 2.3);
    }
  }
  else if (analysis_key == "1pim")
  {
    if (observable == "ECal")
    {
      if (EBeam == 1.161)
        binning = plotting::GetECalBinning(20, 10, 0.6, 1.2, EBeam);
      else if (EBeam == 2.261)
        binning = plotting::GetECalBinning(20, 10, 0.6, EBeam + 0.06, EBeam);
      else if (EBeam == 4.461)
        binning = plotting::GetECalBinning(15, 10, 1.2, EBeam + 0.1, EBeam);
    }
    if (observable == "RecoEvPion")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(50, 0.8, 2);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(50, 0.5, 3.5);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(50, 0, 7);
    }
    if (observable == "RecoWPion")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(50, 0.9, 2);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(50, 0.5, 3);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(50, 0, 4);
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
    else if (observable == "Angleqvshad")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(30, 0, 120);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(30, 0, 120);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(20, 0, 60);
    }
    else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(20, 0, 0.7);
    }
  }
  else if (analysis_key == "1pip")
  {
    if (observable == "ECal")
    {
      if (EBeam == 1.161)
        binning = plotting::GetECalBinning(15, 10, 0.6, 1.2, EBeam);
      else if (EBeam == 2.261)
        binning = plotting::GetECalBinning(15, 10, 1, EBeam + 0.06, EBeam);
      else if (EBeam == 4.461)
        binning = plotting::GetECalBinning(15, 10, 1, EBeam + 0.1, EBeam);
    }
    else if (observable == "pip_mom")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(20, 0., 0.6);
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
        binning = plotting::GetUniformBinning(50, 0.8, 2);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(50, 0.5, 3.5);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(50, 0, 7);
    }
    if (observable == "RecoWPion")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(50, 0.9, 2);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(50, 0.5, 3);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(50, 0, 4);
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
    else if (observable == "Angleqvshad")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(30, 0, 120);
      else if (EBeam == 2.261)
        binning = plotting::GetUniformBinning(30, 0, 120);
      else if (EBeam == 4.461)
        binning = plotting::GetUniformBinning(20, 0, 60);
    }
    else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam == 1.161)
        binning = plotting::GetUniformBinning(20, 0, 0.7);
    }
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
  for (unsigned int i = 0; i < predictions.size(); ++i)
  {
    // if ( max < predictions[i] -> GetMaximum() ) max = predictions[i] -> GetMaximum();
    for (unsigned int j = 0; j < predictions[i]->GetNbinsX(); ++j)
    {
      if (max < predictions[i]->GetBinContent(j))
        max = predictions[i]->GetBinContent(j);
    }
  }
  return max * (1 + 0.12);
}

void plotting::StandardFormat(TH1D *prediction, std::string title, int color, int style, std::string observable, double y_max, std::string y_axis_label)
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

  int FontStyle = 132;
  prediction->GetXaxis()->SetTitleOffset(0.8);
  prediction->GetXaxis()->SetLabelSize(0.08);
  prediction->GetXaxis()->SetTitleSize(0.08);
  prediction->GetXaxis()->SetNdivisions(6);
  prediction->GetXaxis()->SetLabelFont(FontStyle);
  prediction->GetXaxis()->SetTitleFont(FontStyle);

  prediction->GetYaxis()->SetNdivisions(8);
  prediction->GetYaxis()->SetTitleOffset(0.8);
  prediction->GetYaxis()->SetLabelSize(0.08);
  prediction->GetYaxis()->SetTitleSize(0.08);
  prediction->GetYaxis()->SetLabelFont(43);
  prediction->GetYaxis()->SetLabelFont(FontStyle);
  prediction->GetYaxis()->SetTitleFont(FontStyle);
  prediction->GetYaxis()->SetRangeUser(0, y_max);
  prediction->GetYaxis()->SetMaxDigits(3);
  prediction->SetTitleFont(FontStyle);

  return;
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
