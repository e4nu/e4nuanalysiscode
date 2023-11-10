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
  double W = kNucleonMass;
  double minQ2 = 0.01;
  double nbins = 100;
  vector<double> maxQ2,binwidth ;
  for(unsigned int i = 0 ; i < EBeam.size(); ++i ) {
      maxQ2.push_back(2*kNucleonMass*(EBeam[i]-kElectronMass));
      binwidth.push_back((maxQ2[i]-minQ2)/nbins);
  }

  vector<vector<double>> Q2bins, QELVertex, QELVaumm ;
  for( unsigned int i = 0 ; i < EBeam.size() ; ++i ) {
    vector<double> Q2v, QELVertexv, QELVaummv;
    for( unsigned int j = 0 ; j < nbins ; ++j ) {
      double Q2 = minQ2 + j*binwidth[i] ;
      Q2v.push_back(Q2);
      QELVertexv.push_back(utils::RadCorrQELVertex(Q2));
      QELVaummv.push_back(utils::RadCorrQELVacumm(Q2));
    }
    Q2bins.push_back(Q2v);
    QELVaumm.push_back(QELVaummv);
    QELVertex.push_back(QELVertexv);
  }

  TCanvas * c1 = new TCanvas("c1","c1",200,10,700,500);
  TGraph* gRadCorrQELVertex=new TGraph( nbins, &Q2bins[0][0], &QELVertex[0][0]);
  TGraph* gRadCorrQELVacumm=new TGraph( nbins, &Q2bins[0][0], &QELVaumm[0][0]);

  gRadCorrQELVertex->SetLineWidth(2);
  gRadCorrQELVertex->SetLineColor(216);
  gRadCorrQELVertex->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  gRadCorrQELVertex->GetYaxis()->SetTitle("#delta_{vertex}");
  gRadCorrQELVertex->GetXaxis()->CenterTitle();
  gRadCorrQELVertex->GetYaxis()->CenterTitle();
  gRadCorrQELVertex->Draw("AC");

  c1->SaveAs("RadCorrectionVertex.pdf");

  gRadCorrQELVacumm->SetLineWidth(2);
  gRadCorrQELVacumm->SetLineColor(209);
  gRadCorrQELVacumm->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  gRadCorrQELVacumm->GetYaxis()->SetTitle("#delta_{vac}");
  gRadCorrQELVacumm->GetXaxis()->CenterTitle();
  gRadCorrQELVacumm->GetYaxis()->CenterTitle();
  gRadCorrQELVacumm->Draw("AC");

  c1->SaveAs("RadCorrectionVacumm.pdf");
 return 0 ;
}
