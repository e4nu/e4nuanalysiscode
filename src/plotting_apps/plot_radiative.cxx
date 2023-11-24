// __________________________________________________________________________
/* This app is used to plot quantities related to radiative corrections    */
// __________________________________________________________________________
#include <iostream>
#include <vector>
#include "utils/RadiativeCorrUtils.h"
#include "conf/ConstantsI.h"
#include "conf/ParticleI.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting e4nu analysis radiative corrections..." << std::endl;

  vector<double> EBeam = {1,2,4};
  double minQ2 = 0.01;
  double nbins = 100;
  int Z = 1 ; // To change
  vector<double> maxQ2,binwidth ;
  for(unsigned int i = 0 ; i < EBeam.size(); ++i ) {
    maxQ2.push_back(2*kNucleonMass*(EBeam[i]-kElectronMass));
    binwidth.push_back((maxQ2[i]-minQ2)/nbins);
  }

  vector<vector<double>> Q2bins, Epbins, QELVertex, QELVaumm, RadCorrQELRealRad, QELP1, QEL_Total ;
  for( unsigned int i = 0 ; i < EBeam.size() ; ++i ) {
    vector<double> Q2v, Epv, QELVertexv, QELVaummv, QELRealRadv, QELRealP1v, QEL_totalv;
    for( unsigned int j = 0 ; j < nbins ; ++j ) {
      double Q2 = minQ2 + j*binwidth[i] ;
      double Ep = EBeam[i]-Q2*0.5/kNucleonMass;
      Q2v.push_back(Q2);
      Epv.push_back(Ep);
      QELVertexv.push_back(utils::RadCorrQELVertex(Q2));
      QELVaummv.push_back(utils::RadCorrQELVacumm(Q2));
      QELRealRadv.push_back(utils::RadCorrQELRealRad(Q2,EBeam[i],EBeam[i]-0.1,0));
      QEL_totalv.push_back(1+QELVertexv[j]+QELVaummv[j]+QELRealRadv[i]);
    }
    Q2bins.push_back(Q2v);
    Epbins.push_back(Epv);
    QELVaumm.push_back(QELVaummv);
    QELVertex.push_back(QELVertexv);
    RadCorrQELRealRad.push_back(QELRealRadv);
  }

  TCanvas * c1 = new TCanvas("c1","c1",200,10,700,500);
  TGraph* gRadCorrQELTotal_1=new TGraph( nbins, &Q2bins[0][0], &QEL_Total[0][0]);
  TGraph* gRadCorrQELVertex_1=new TGraph( nbins, &Q2bins[0][0], &QELVertex[0][0]);
  TGraph* gRadCorrQELVacumm_1=new TGraph( nbins, &Q2bins[0][0], &QELVaumm[0][0]);
  TGraph* gRadCorrQELRealRad_1=new TGraph( nbins, &Q2bins[0][0], &RadCorrQELRealRad[0][0]);

  TGraph* gRadCorrQELTotal_Ep_1=new TGraph( nbins, &Epbins[0][0], &QEL_Total[0][0]);
  TGraph* gRadCorrQELVertex_Ep_1=new TGraph( nbins, &Epbins[0][0], &QELVertex[0][0]);
  TGraph* gRadCorrQELVacumm_Ep_1=new TGraph( nbins, &Epbins[0][0], &QELVaumm[0][0]);

  TGraph* gRadCorrQELRealRad_2=new TGraph( nbins, &Q2bins[1][0], &RadCorrQELRealRad[1][0]);
  TGraph* gRadCorrQELRealRad_4=new TGraph( nbins, &Q2bins[2][0], &RadCorrQELRealRad[2][0]);

  gRadCorrQELTotal_1->SetLineWidth(2);
  gRadCorrQELTotal_1->SetLineColor(kBlack);
  gRadCorrQELTotal_1->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  gRadCorrQELTotal_1->GetYaxis()->SetTitle("#delta");
  gRadCorrQELTotal_1->GetXaxis()->CenterTitle();
  gRadCorrQELTotal_1->GetYaxis()->CenterTitle();

  gRadCorrQELVertex_1->SetLineWidth(2);
  gRadCorrQELVertex_1->SetLineColor(216);
  gRadCorrQELVertex_1->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  gRadCorrQELVertex_1->GetYaxis()->SetTitle("#delta_{vertex}");
  gRadCorrQELVertex_1->GetXaxis()->CenterTitle();
  gRadCorrQELVertex_1->GetYaxis()->CenterTitle();
  gRadCorrQELVertex_1->Draw("AC");
  c1->SaveAs("RadCorrectionVertex.pdf");

  gRadCorrQELTotal_Ep_1->SetLineWidth(2);
  gRadCorrQELTotal_Ep_1->SetLineColor(kBlack);
  gRadCorrQELTotal_Ep_1->GetXaxis()->SetTitle("E'[GeV]");
  gRadCorrQELTotal_Ep_1->GetYaxis()->SetTitle("#delta");
  gRadCorrQELTotal_Ep_1->GetXaxis()->CenterTitle();
  gRadCorrQELTotal_Ep_1->GetYaxis()->CenterTitle();

  gRadCorrQELVertex_Ep_1->SetLineWidth(2);
  gRadCorrQELVertex_Ep_1->SetLineColor(216);
  gRadCorrQELVertex_Ep_1->GetXaxis()->SetTitle("E'[GeV]");
  gRadCorrQELVertex_Ep_1->GetYaxis()->SetTitle("#delta_{vertex}");
  gRadCorrQELVertex_Ep_1->GetXaxis()->CenterTitle();
  gRadCorrQELVertex_Ep_1->GetYaxis()->CenterTitle();
  gRadCorrQELVertex_Ep_1->Draw("AC");
  c1->SaveAs("RadCorrectionVertex_Ep.pdf");

  gRadCorrQELVacumm_1->SetLineWidth(2);
  gRadCorrQELVacumm_1->SetLineColor(209);
  gRadCorrQELVacumm_1->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  gRadCorrQELVacumm_1->GetYaxis()->SetTitle("#delta_{vac}");
  gRadCorrQELVacumm_1->GetXaxis()->CenterTitle();
  gRadCorrQELVacumm_1->GetYaxis()->CenterTitle();
  gRadCorrQELVacumm_1->Draw("AC");
  c1->SaveAs("RadCorrectionVacumm.pdf");

  gRadCorrQELVacumm_Ep_1->SetLineWidth(2);
  gRadCorrQELVacumm_Ep_1->SetLineColor(209);
  gRadCorrQELVacumm_Ep_1->GetXaxis()->SetTitle("E'[GeV]");
  gRadCorrQELVacumm_Ep_1->GetYaxis()->SetTitle("#delta_{vac}");
  gRadCorrQELVacumm_Ep_1->GetXaxis()->CenterTitle();
  gRadCorrQELVacumm_Ep_1->GetYaxis()->CenterTitle();
  gRadCorrQELVacumm_Ep_1->Draw("AC");
  c1->SaveAs("RadCorrectionVacumm_Ep.pdf");

  gRadCorrQELTotal_1->Draw("AC");
  c1->SaveAs("RadCorrectionTotal.pdf");

  gRadCorrQELTotal_Ep_1->Draw("AC");
  c1->SaveAs("RadCorrectionTotal_Ep.pdf");


  gRadCorrQELRealRad_1->SetLineWidth(2);
  gRadCorrQELRealRad_1->SetLineColor(kPink);
  gRadCorrQELRealRad_1->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  gRadCorrQELRealRad_1->GetYaxis()->SetTitle("#delta_{R}");
  gRadCorrQELRealRad_1->GetXaxis()->CenterTitle();
  gRadCorrQELRealRad_1->GetYaxis()->CenterTitle();
  gRadCorrQELRealRad_1->Draw("AC");

  gRadCorrQELRealRad_2->SetLineWidth(2);
  gRadCorrQELRealRad_2->SetLineColor(kGreen);
  gRadCorrQELRealRad_2->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  gRadCorrQELRealRad_2->GetYaxis()->SetTitle("#delta_{R}");
  gRadCorrQELRealRad_2->GetXaxis()->CenterTitle();
  gRadCorrQELRealRad_2->GetYaxis()->CenterTitle();
  gRadCorrQELRealRad_2->Draw("same");
  c1->SaveAs("gRadCorrQELRealRad.pdf");

  gRadCorrQELRealRad_4->SetLineWidth(2);
  gRadCorrQELRealRad_4->SetLineColor(kBlue);
  gRadCorrQELRealRad_4->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  gRadCorrQELRealRad_4->GetYaxis()->SetTitle("#delta_{R}");
  gRadCorrQELRealRad_4->GetXaxis()->CenterTitle();
  gRadCorrQELRealRad_4->GetYaxis()->CenterTitle();
  gRadCorrQELRealRad_4->Draw("same");
  c1->SaveAs("gRadCorrQELRealRad.pdf");

  return 0 ;
}
