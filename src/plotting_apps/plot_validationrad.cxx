// __________________________________________________________________________
/* This app is used to plot the validation plot from the e4nu paper figure 5*/
// __________________________________________________________________________
#include <iostream>
#include <vector>
#include "utils/RadiativeCorrUtils.h"
#include "conf/ConstantsI.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"
#include "utils/TargetUtils.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"
#include "analysis/MCCLAS6AnalysisI.h"
#include "TGraph.h"
#include "TH1D.h"
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

  std::cout << "Plotting e4nu analysis validation plot for radiative corrections..." << std::endl;

  vector<double> x_axis = {4.3005,4.3016,4.3028,4.3039,4.3051,4.3062,4.3074,4.3085,4.3096, 4.3108,4.3119,4.31301,4.31414, 4.3153,4.3164,4.3176,4.3187,4.3198,4.3210,4.3221,4.3232,4.3244,4.3255,4.3267,4.3278,4.3289,4.33};
  vector<double> y_axis = {142.86,142.86,178.57,214.29,214.29,178.57,214.29,214.29,250,250,285.71,321.43,357.14,392.86,428.57,500,571.43,714.29,821.43,1035.71,1428.57,2607.14,8392.86,6250,321.43,107.14,71.43};

  TCanvas * c1 = new TCanvas("c1","c1",200,10,700,500);
  TGraph * g_data = new TGraph(x_axis.size(), &x_axis[0], &y_axis[0] );

  g_data->SetLineWidth(2);
  g_data->SetLineColor(kBlack);
  g_data->GetXaxis()->SetTitle("^{1}H(e,e'p)E_{Cal}[GeV]");
  g_data->GetYaxis()->SetTitle("#events");
  g_data->GetXaxis()->CenterTitle();
  g_data->GetYaxis()->CenterTitle();
  g_data->SetLineWidth(4);
  g_data->SetMarkerColor(kBlack);
  g_data->SetMarkerSize(2);
  g_data->SetMarkerStyle(kFullCircle);

  TH1D * h_data = new TH1D( "hdata","hdata", 27, 4.3, 4.33 );
  auto nPoints = g_data->GetN(); // number of points in your TGraph
  for(int i=0; i < nPoints; ++i) {
    double x,y;
    g_data->GetPoint(i, x, y);
    h_data->Fill(x,y);
  }
  double data_integral = h_data->Integral();

  // Compute prediction
  string input_file = "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/MonoFlux_G18_10a_H_4320MeV.gst.root";
  int nevents = 1000;
  TH1D * h_sim1 = new TH1D( "hsim1","hsim1", 27,4.3, 4.33 );

  h_sim1->SetLineWidth(2);
  h_sim1->SetLineColor(kBlue);
  h_sim1->GetXaxis()->SetTitle("^{1}H(e,e'p)E_{Cal}[GeV]");
  h_sim1->GetYaxis()->SetTitle("#events");
  h_sim1->GetXaxis()->CenterTitle();
  h_sim1->GetYaxis()->CenterTitle();
  h_sim1->SetLineWidth(4);

  MCEventHolder * event_holder = new MCEventHolder( input_file, 0, nevents ) ; 
  if( !event_holder ) {
    std::cout << "Failed to instiantize event holder" << std::endl;
    return 0;
  }

  event_holder -> LoadBranch();
  for ( unsigned int i = 0 ; i < event_holder->GetNEvents() ; ++i ) { 
    e4nu::Event & event = *event_holder->GetEvent(i);
    std::map<int,std::vector<TLorentzVector>> particle_map = event.GetFinalParticles4Mom();
    std::map<int,std::vector<TLorentzVector>> new_particle_map;
    bool is_1p = true ;
    
    double ECal = event.GetOutLepton4Mom().E();
    for( auto it = particle_map.begin() ; it != particle_map.end() ; ++it ) {
      // Calculate ECal for visible particles
      
      if( it-> first == kPdgPiM && it->second.size() != 0 ) break ;
      if( it-> first == kPdgProton && it->second.size() ==1 ) {
	for( unsigned int j = 0 ; j < (it->second).size() ; ++j ) {
	  ECal += (it->second)[j].E() + utils::GetBindingEnergy( kPdgH ) - utils::GetParticleMass( it->first ) ; // Correct for proton binding energy
	}
      } else is_1p = false ;
    }

    if (!is_1p) continue ;
    h_sim1 -> Fill( ECal, event.GetEventWeight() );
  }
  
  double integral_sim1 = h_sim1->Integral();
  h_sim1->Scale(data_integral/integral_sim1);

  TH1D * h_sim2 = new TH1D( "hsim2","hsim2", 27,4.3, 4.33 );  
  h_sim2->SetLineWidth(2);
  h_sim2->SetLineColor(kRed);
  h_sim2->GetXaxis()->SetTitle("^{1}H(e,e'p)E_{Cal}[GeV]");
  h_sim2->GetYaxis()->SetTitle("#events");
  h_sim2->GetXaxis()->CenterTitle();
  h_sim2->GetYaxis()->CenterTitle();
  h_sim2->SetLineWidth(4);

  //string input_file_2 = "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/RadFlux_G18_10a_H_4320MeV.gst.root";
  string input_file_2 = "/genie/app/users/jtenavid/Software/e4v/E4NuAnalysis/Source/e4nuanalysiscode/radiated.gst.root";
  MCEventHolder * event_holder_2 = new MCEventHolder( input_file_2, 0, nevents ) ; 
  if( !event_holder_2 ) {
    std::cout << "Failed to instiantize event holder" << std::endl;
    return 0;
  }

  event_holder_2 -> LoadBranch();
  for ( unsigned int i = 0 ; i < event_holder_2->GetNEvents() ; ++i ) { 
    e4nu::Event & event = *event_holder_2->GetEvent(i);
    std::map<int,std::vector<TLorentzVector>> particle_map = event.GetFinalParticles4Mom();
    std::map<int,std::vector<TLorentzVector>> new_particle_map;
    bool is_1p = true ;
    
    double ECal = event.GetOutLepton4Mom().E();
    for( auto it = particle_map.begin() ; it != particle_map.end() ; ++it ) {
      // Calculate ECal for visible particles
      std::cout << it->first<< std::endl;
      if( it-> first == kPdgPiM && it->second.size() != 0 ) break ;
      if( it-> first == kPdgProton && it->second.size() ==1 ) {
	for( unsigned int j = 0 ; j < (it->second).size() ; ++j ) {
	  ECal += (it->second)[j].E() + utils::GetBindingEnergy( kPdgH ) - utils::GetParticleMass( it->first ) ; // Correct for proton binding energy
	  std::cout << "proton = " << (it->second)[j].E() << std::endl;
	}
      } else is_1p = false ;
    }
    
    ///    if (!is_1p) continue ;
    std::cout << ECal << " mom " << event.GetOutLepton4Mom().E() <<std::endl;
    h_sim2 -> Fill( ECal, event.GetEventWeight() );
  }

  double integral_sim2 = h_sim2->Integral();
  h_sim2->Scale(data_integral/integral_sim2);

  g_data->Draw("AP");
  //  //  h_data->Draw("hist same");
  h_sim1->Draw("hist same");
  h_sim2->Draw("hist same");
  c1->SaveAs("RadValidation_H_4320MeV_1p0pi.pdf");

  return 0;
}
