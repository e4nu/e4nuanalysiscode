// __________________________________________________________________________
/* This app is used to process a GENIE MC gst file with radiation weights  */
// __________________________________________________________________________

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include "TFile.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TLegend.h"
#include "analysis/MCCLAS6AnalysisI.h"
#include "plotting/PlottingUtils.h"
#include "utils/RadiativeCorrUtils.h"
#include "conf/RadConstants.h"
#include "utils/Utils.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf; 
using namespace e4nu::plotting;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////
// --input-gst-file)  InputFile in GENIE gst format            //
// --output-gst-file) OutputFile in GENIE gst format           //
// --true-EBeam)      True beam energy used for flux generation//
// --nevents          Number of events to process              //
// --first-event      Default 0                                //
// --target : target pdg                                       //
// --thickness : thickness of your experiment target           //
// --rad-model : model used for external radiation of in e-    //
//               simc or simple                                //
// --max-Ephoton : defaulted to 0.2 (of beam energy)           //
// --output-hist : name of output file with histograms         //
/////////////////////////////////////////////////////////////////
TTree * s_tree ;
void StoreToGstFormat( Event & event, string output_file ) ;

int main( int argc, char* argv[] ) {

  std::cout << "Processing GENIE gst file. Adding radiation weights and correcting electron kinematics..." << std::endl;

  string input_file, output_file = "radiated.gst.root";
  string output_hist = "" ; 
  double nevents = 10000000;
  double first_event = 0 ;
  double true_beam_energy = 1; 
  int tgt = 1000060120 ;
  double thickness = e4nu::conf::GetThickness(tgt); // Defaulted to CLAS6
  string rad_model = "simc";
  double MaxEPhoton = 0.2 ;

  if( argc == 0 ) { 
    cout << " Please, specify input gst file. Abort. " << endl;
    return 0;
  }

  if( argc > 1 ) { 
    if( ExistArg("input-gst-file",argc,argv)) {
      input_file = GetArg("input-gst-file",argc,argv) ;
    } else { 
      cout << " Please, specify input gst file. Abort. " << endl;
      return 0;
    }
    if( ExistArg("true-EBeam",argc,argv)){
      true_beam_energy = stod(GetArg("true-EBeam",argc,argv));
    } else { 
      cout << " Please specify true beam energy used for flux generation" << endl;
      return 0;
    }
    if( ExistArg("output-gst-file",argc,argv)) {
      output_file = GetArg("output-gst-file",argc,argv) ;
    } 
    if( ExistArg("nevents",argc,argv)) {
      nevents = stoi(GetArg("nevents",argc,argv));
    }
    if( ExistArg("first-event",argc,argv)) {
      first_event = stoi(GetArg("first-event",argc,argv));
    }
    if( ExistArg("output-hist",argc,argv)) {
      output_hist = GetArg("output-hist",argc,argv);
    }

    // Radiative specific options
    if( ExistArg("rad-model",argc,argv)) {
      rad_model = GetArg("rad-model",argc,argv); 
    }
    if( ExistArg("target",argc,argv)) {
      tgt = stoi(GetArg("target",argc,argv)); 
      thickness = e4nu::conf::GetThickness(tgt); // Defaulted to CLAS6
    }
    if( ExistArg("thickness",argc,argv)) {
      thickness = stoi(GetArg("thickness",argc,argv)); 
    }
    if( ExistArg("MaxEPhoton",argc,argv)) {
      MaxEPhoton = stod(GetArg("MaxEPhoton",argc,argv)); 
    }
    if( ExistArg("max-Ephoton",argc,argv)) {
      MaxEPhoton = stod(GetArg("max-Ephoton",argc,argv)); 
    }
    
  }

  TFile * out_file = TFile::Open(output_file.c_str(),"RECREATE");
  s_tree = new TTree("gst","GENIE Summary Event Tree");

  MCEventHolder * event_holder = new MCEventHolder( input_file, first_event, nevents ) ; 
  if( !event_holder ) {
    std::cout << "Failed to instiantize event holder" << std::endl;
    return 0;
  }

  // Initialize histograms for debugging
  TH1D * h_outemom = new TH1D( "h_outemom", "h_outemom", 100, 0, true_beam_energy) ; 
  TH1D * h_outcorremom = new TH1D( "h_outcorremom", "h_outcorremom", 100, 0, true_beam_energy) ; 
  TH1D * h_outgammamom = new TH1D( "h_outgammamom", "h_outgammamom", 100, 0, true_beam_energy) ; 
  
  event_holder -> LoadBranch();
  for ( unsigned int i = 0 ; i < event_holder->GetNEvents() ; ++i ) { 
    if( i==0 || i % 10000 == 0 ) utils::PrintProgressBar( i, event_holder->GetNEvents() ) ;
    e4nu::Event & event = *event_holder->GetEvent(i);
    
    // Store the in and out going electron kinematics at the interaction vertex
    // We refer to these as corrected electrons
    // Detected electrons correspond to the beam electron and final detected electron
    TLorentzVector CorrInElectron  = event.GetInLepton4Mom() ; // From GENIE file it corresponds to the interaction electron (already radiated)
    TLorentzVector CorrOutElectron = event.GetOutLepton4Mom() ; // From GENIE file it corresponds to the interaction electron (before radiation)

    event.SetInCorrLeptonKinematics( CorrInElectron.E(), CorrInElectron.Px(), CorrInElectron.Py(), CorrInElectron.Pz() ) ; 
    event.SetOutCorrLeptonKinematics( CorrOutElectron.E(), CorrOutElectron.Px(), CorrOutElectron.Py(), CorrOutElectron.Pz() ) ; 

    // Set true incoming electron kinematics from configuration
    event.SetInLeptonKinematics( true_beam_energy, 0, 0, true_beam_energy ) ; 

    // Compute energy of emited photon 
    TLorentzVector InGamma = event.GetInLepton4Mom() - CorrInElectron ;
    if( InGamma.E() < 0 ) InGamma.SetPxPyPzE(0,0,0,0);

    // Compute true detected outgoing electron kinematics with energy loss method
    TLorentzVector detected_electron; // will be set in RadOutElectron
    TLorentzVector OutGamma = RadOutElectron( CorrOutElectron, detected_electron,tgt, thickness, MaxEPhoton, "simc");//rad_model );
    h_outgammamom->Fill(OutGamma.E()); 
    h_outcorremom->Fill(CorrOutElectron.E());
    h_outemom->Fill(detected_electron.E());

    // Set true outcoming electron kinematics from configuration  
    event.SetOutLeptonKinematics( detected_electron ) ;

    // Compute correction weight
    double weight = SIMCRadCorrWeight( event, thickness, MaxEPhoton, rad_model ) ;
    
    event.SetEventWeight ( weight * event.GetEventWeight() ) ; 

    // Store emited photons in event record
    std::map<int,std::vector<TLorentzVector>> event_record = event.GetFinalParticles4Mom();
    if( event_record.find(kPdgPhoton) != event_record.end()){
      event_record[kPdgPhoton].push_back( InGamma ) ;
      event_record[kPdgPhoton].push_back( OutGamma ) ;
    } else { 
      event_record[kPdgPhoton] = { InGamma, OutGamma } ;
    }  
    event.SetFinalParticlesKinematics(event_record);

    // Store back to gst format
    StoreToGstFormat( event, output_file ) ;
  }

  s_tree ->Write();

  // Store debug info in histogram 
  if( output_hist != "" ) {
    TCanvas * c  = new TCanvas("","",800,800);
    c->SetFillColor(0);

    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetStatColor(0);
    gStyle->SetFillColor(0);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetPaperSize(20,26);
    gStyle->SetTitleFont(132,"pad");
    gStyle->SetMarkerStyle(20);
    gStyle->SetLineStyleString(2,"[12 12]");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
        
    h_outemom->SetLineColor(kRed);
    h_outemom->SetLineWidth(2);
    h_outcorremom->SetLineColor(kBlue);
    h_outcorremom->SetLineWidth(2);

    h_outcorremom -> GetXaxis()->SetTitle("E_{e}^{det}[GeV]");
    h_outcorremom -> GetYaxis()->SetTitle("#Events");
    h_outcorremom -> GetXaxis()->CenterTitle();
    h_outcorremom -> GetYaxis()->CenterTitle();
    
    h_outcorremom->Draw("hist");
    h_outemom->Draw("hist same");

    double LegXmin = 0.1, LegYmin = 0.65, YSpread = 0.25;
    TLegend* leg = new TLegend(LegXmin,LegYmin,LegXmin+0.9,LegYmin+YSpread);
    leg->AddEntry(h_outcorremom, "Vertex Electron");
    leg->AddEntry(h_outemom, "Detected Electron");
    leg->Draw();
    c->SaveAs(output_hist.c_str());
    delete c ;	       
  }

  out_file -> Close();
  
  delete out_file;
  delete event_holder;
  return 0 ; 
}

void StoreToGstFormat( Event & event, string output_file ) {
  //____________________________________________________________________________________
  // GENIE GHEP EVENT TREE FORMAT -> GENIE SUMMARY NTUPLE 
  // Note not all the information is propagated... 
  //____________________________________________________________________________________
  
  // Members for root file
  Int_t iev = event.GetEventID();
  Int_t neu = 0 ;
  Int_t fspl = 0 ;
  Int_t tgt = 0 ;
  Int_t A = 0 ;
  Int_t Z = 0 ;
  Int_t hitnuc = 0 ;
  Int_t hitqrk = 0 ;
  Bool_t qel = false ;
  Bool_t mec = false ;
  Bool_t res = false ;
  Bool_t dis = false ;
  Bool_t em = false ;
  Bool_t cc = false ;
  Bool_t nc = false ;
  Bool_t coh = false ;
  Bool_t dfr = false ;
  Bool_t imd = false ;
  Bool_t norm = false ; 
  Bool_t imdanh = false ;
  Bool_t singlek = false ;
  Bool_t nuel = false ;
  Double_t wght = 0 ;
  Double_t xs = 0 ;
  Double_t ys = 0 ;
  Double_t ts = 0 ;
  Double_t Q2s = 0 ;
  Double_t Ws = 0 ;
  Double_t x = 0 ;
  Double_t y = 0 ;
  Double_t t = 0 ;
  Double_t Q2 = 0 ;
  Double_t W = 0 ;
  Double_t Ev = 0 ;
  Double_t pxv = 0 ;
  Double_t pyv = 0 ;
  Double_t pzv = 0 ;
  Double_t En = 0 ;
  Double_t pxn = 0 ;
  Double_t pyn = 0 ;
  Double_t pzn = 0 ;
  Double_t El = 0 ;
  Double_t pxl = 0 ;
  Double_t pyl = 0 ;
  Double_t pzl = 0 ;
  Double_t pl = 0 ;
  Double_t cthl = 0 ;
  Double_t Evcorr = 0 ;
  Double_t pxvcorr = 0 ;
  Double_t pyvcorr = 0 ;
  Double_t pzvcorr = 0 ;
  Double_t Elcorr = 0 ;
  Double_t pxlcorr = 0 ;
  Double_t pylcorr  = 0 ;
  Double_t pzlcorr = 0 ;
  Int_t  nfp = 0 ;
  Int_t  nfn = 0 ;
  Int_t  nfpip = 0 ;
  Int_t  nfpim = 0 ;
  Int_t  nfpi0 = 0 ;
  Int_t  nfkp = 0 ;
  Int_t  nfkm = 0 ;
  Int_t  nfk0 = 0 ;
  Int_t  nfem = 0 ;
  Int_t  nfother = 0 ;
  Int_t nf = 0 ;
  Int_t pdgf[120] ;
  Double_t Ef[120] ;
  Double_t pxf[120] ;
  Double_t pyf[120] ;
  Double_t pzf[120] ;
  Double_t pf[120]; 
  Double_t cthf[120] ;
  Double_t vtxx = 0 ;
  Double_t vtxy = 0 ;
  Double_t vtxz = 0 ;
  Double_t vtxt = 0 ;
  Int_t resid = -1000; 
  Int_t nin = 0 ; 
  Int_t nipip = 0 ; 
  Int_t nipim = 0 ; 
  Int_t nipi0 = 0 ; 
  Int_t nikp = 0 ; 
  Int_t nikm = 0 ; 
  Int_t nik0 = 0 ; 
  Int_t niem = 0 ; 
  Int_t niother = 0 ; 
  Int_t ni = 0 ; 
  Int_t pdgi[120] ; 
  Double_t Ei[120] ;
  Double_t pxi[120] ;
  Double_t pyi[120] ;
  Double_t pzi[120] ;


  // Create tree branches
  if( event.GetEventID() == 0 ){
    s_tree->Branch("iev", &iev, "iev/I");
    s_tree->Branch("neu", &neu, "neu/I");
    s_tree->Branch("fspl", &fspl, "fspl/I");
    s_tree->Branch("tgt", &tgt, "tgt/I");
    s_tree->Branch("A", &A, "A/I");
    s_tree->Branch("Z", &Z, "Z/I");
    s_tree->Branch("hitnuc", &hitnuc, "hitnuc/I");
    s_tree->Branch("hitqrk", &hitqrk, "hitqrk/I");
    s_tree->Branch("qel", &qel, "qel/O");
    s_tree->Branch("mec", &mec, "mec/O");//go down from here
    s_tree->Branch("res", &res, "res/O");
    s_tree->Branch("dis", &dis, "dis/O");
    s_tree->Branch("em", &em, "em/O");
    s_tree->Branch("cc", &cc, "cc/O");
    s_tree->Branch("nc", &nc, "nc/O");
    s_tree->Branch("coh", &coh, "coh/O");
    s_tree->Branch("dfr", &dfr, "dfr/O");
    s_tree->Branch("imd", &imd, "imd/O");
    s_tree->Branch("norm", &norm, "norm/O");
    s_tree->Branch("imdanh", &imdanh, "imdanh/O");
    s_tree->Branch("singlek", &singlek, "singlek/O");
    s_tree->Branch("nuel", &nuel, "nuel/O");
    s_tree->Branch("wght", &wght, "wght/D");
    s_tree->Branch("xs", &xs, "xs/D");
    s_tree->Branch("ys", &ys, "ys/D");
    s_tree->Branch("ts", &ts, "ts/D");
    s_tree->Branch("Q2s", &Q2s, "Q2s/D");
    s_tree->Branch("Ws", &Ws, "Ws/D");
    s_tree->Branch("x", &x, "x/D");
    s_tree->Branch("y", &y, "y/D");
    s_tree->Branch("t", &t, "t/D");
    s_tree->Branch("Q2", &Q2, "Q2/D");
    s_tree->Branch("W", &W, "W/D");
    s_tree->Branch("Ev", &Ev, "Ev/D");
    s_tree->Branch("pxv", &pxv, "pxv/D");
    s_tree->Branch("pyv", &pyv, "pyv/D");
    s_tree->Branch("pzv", &pzv, "pzv/D");
    s_tree->Branch("En", &En, "En/D");
    s_tree->Branch("pxn", &pxn, "pxn/D");
    s_tree->Branch("pyn", &pyn, "pyn/D");
    s_tree->Branch("pzn", &pzn, "pzn/D");
    s_tree->Branch("El", &El, "El/D");
    s_tree->Branch("pxl", &pxl, "pxl/D");
    s_tree->Branch("pyl", &pyl, "pyl/D");
    s_tree->Branch("pzl", &pzl, "pzl/D");
    s_tree->Branch("pl", &pl, "pl/D");
    s_tree->Branch("cthl", &cthl, "cthl/D");
    s_tree->Branch("Evcorr", &Evcorr, "Evcorr/D");
    s_tree->Branch("pxvcorr", &pxvcorr, "pxvcorr/D");
    s_tree->Branch("pyvcorr", &pyvcorr, "pyvcorr/D");
    s_tree->Branch("pzvcorr", &pzvcorr, "pzvcorr/D");
    s_tree->Branch("Elcorr", &Elcorr, "Elcorr/D");
    s_tree->Branch("pxlcorr", &pxlcorr, "pxlcorr/D");
    s_tree->Branch("pylcorr", &pylcorr, "pylcorr/D");
    s_tree->Branch("pzlcorr", &pzlcorr, "pzlcorr/D");
    s_tree->Branch("nfp", &nfp, "nfp/I");
    s_tree->Branch("nfn", &nfn, "nfn/I");
    s_tree->Branch("nfpip", &nfpip, "nfpip/I");
    s_tree->Branch("nfpim", &nfpim, "nfpim/I");
    s_tree->Branch("nfpi0", &nfpi0, "nfpi0/I");
    s_tree->Branch("nfkp", &nfkp, "nfkp/I");
    s_tree->Branch("nfkm", &nfkm, "nfkm/I");
    s_tree->Branch("nfk0", &nfk0, "nfk0/I");
    s_tree->Branch("nfem", &nfem, "nfem/I");
    s_tree->Branch("nfother", &nfother, "nfother/I");
    s_tree->Branch("nf", &nf, "nf/I");
    s_tree->Branch("pdgf", pdgf, "pdgf[nf]/I");
    s_tree->Branch("Ef", Ef, "Ef[nf]/D");
    s_tree->Branch("pf", pf, "pf[nf]/D");
    s_tree->Branch("cthf", cthf, "cthf[nf]/D");
    s_tree->Branch("pxf", pxf, "pxf[nf]/D");
    s_tree->Branch("pyf", pyf, "pyf[nf]/D");
    s_tree->Branch("pzf", pzf, "pzf[nf]/D");
    s_tree->Branch("vtxx", &vtxx, "vtxx/D");
    s_tree->Branch("vtxy", &vtxy, "vtxy/D");
    s_tree->Branch("vtxz", &vtxz, "vtxz/D");
    s_tree->Branch("vtxt", &vtxt, "vtxt/D");
    s_tree->Branch("resid", &resid, "resid/I");
    s_tree->Branch("nin", &nin, "nin/I");
    s_tree->Branch("nipip", &nipip, "nipip/I");
    s_tree->Branch("nipim", &nipim, "nipim/I");
    s_tree->Branch("nipi0", &nipi0, "nipi0/I");
    s_tree->Branch("nikp", &nikp, "nikp/I");
    s_tree->Branch("nikm", &nikm, "nikm/I");
    s_tree->Branch("nik0", &nik0, "nik0/I");
    s_tree->Branch("niem", &niem, "niem/I");
    s_tree->Branch("niother", &niother, "niother/I");
    s_tree->Branch("ni", &ni, "ni/I");
    s_tree->Branch("pdgi", &pdgi, "pdgi[ni]/I");
    s_tree->Branch("Ei", &Ei, "Ei[ni]/D");
    s_tree->Branch("pxi", &pxi, "pxi[ni]/D");
    s_tree->Branch("pyi", &pyi, "pyi[ni]/D");
    s_tree->Branch("pzi", &pzi, "pzi[ni]/D");
  }

  tgt = event.GetTargetPdg();
  qel = event.IsQEL();
  mec = event.IsMEC();
  res = event.IsRES();
  dis = event.IsDIS();
  em  = event.IsEM();
  cc  = event.IsCC();
  nc  = event.IsNC();
  wght = event.GetEventWeight();
  xs   = event.GetTruexs();
  ys   = event.GetTrueys();
  Q2s  = event.GetTrueQ2s();
  Ws   = event.GetTrueWs();
  x    = event.GetTruex();
  y    = event.GetTruey();
  Q2   = event.GetTrueQ2();
  W    = event.GetTrueW();
  Ev   = event.GetInLepton4Mom().E();
  pxv  = event.GetInLepton4Mom().Px();
  pyv  = event.GetInLepton4Mom().Py();
  pyv  = event.GetInLepton4Mom().Pz();
  El   = event.GetOutLepton4Mom().E();
  pl   = event.GetOutLepton4Mom().P();
  cthl = event.GetOutLepton4Mom().CosTheta();  
  pxl  = event.GetOutLepton4Mom().Px();
  pyl  = event.GetOutLepton4Mom().Py();
  pzl  = event.GetOutLepton4Mom().Pz();
  Evcorr   = event.GetInCorrLepton4Mom().E();
  pxvcorr  = event.GetInCorrLepton4Mom().Px();
  pyvcorr  = event.GetInCorrLepton4Mom().Py();
  pyvcorr  = event.GetInCorrLepton4Mom().Pz();
  Elcorr   = event.GetOutCorrLepton4Mom().E();
  pxlcorr  = event.GetOutCorrLepton4Mom().Px();
  pylcorr  = event.GetOutCorrLepton4Mom().Py();
  pzlcorr  = event.GetOutCorrLepton4Mom().Pz();
  nfp  = event.GetTrueNProtons();
  nfn  = event.GetTrueNNeutrons();
  nfpip = event.GetTrueNPiP();
  nfpim = event.GetTrueNPiM();
  nfpi0 = event.GetTrueNPi0();
  nfkp  = event.GetTrueNKP();
  nfkm  = event.GetTrueNKM();
  nfk0  = event.GetTrueNK0();
  nfem  = event.GetTrueNEM();
  nfother = event.GetTrueNOther();

  nf = 0;
  for( auto k = 0 ; k < 120 ; ++k ) {
    pdgf[k] = 0;
    Ef[k] = 0;
    pxf[k] = 0;
    pyf[k] = 0; 
    pzf[k] = 0;
    pf[k] = 0;
    cthf[k] = 0;
  }

  ni = 0 ;
  for( auto k = 0 ; k < 120 ; ++k ) {
    pdgi[k] = 0;
    Ei[k] = 0;
    pxi[k] = 0;
    pyi[k] = 0; 
    pzi[k] = 0;
  }

  std::map<int,std::vector<TLorentzVector>> final_particles = event.GetFinalParticles4Mom();
  for ( auto it = final_particles.begin(); it != final_particles.end(); it++) {
    for ( unsigned int i = 0 ; i < (it->second).size() ; ++i ) { 
      pdgf[nf] = it->first;
      Ef[nf] = (it->second)[i].E();
      pxf[nf] = (it->second)[i].Px();
      pyf[nf] = (it->second)[i].Py();
      pzf[nf] = (it->second)[i].Pz();
      pf[nf] = (it->second)[i].P();
      cthf[nf] = (it->second)[i].CosTheta();
      ++nf ;
    }
  }
    
  vtxx = event.GetVertex().X();
  vtxy = event.GetVertex().Y();
  vtxz = event.GetVertex().Z();
  vtxt = event.GetVertex().T();
  resid = event.GetRESID();
  s_tree -> Fill();
 
}
