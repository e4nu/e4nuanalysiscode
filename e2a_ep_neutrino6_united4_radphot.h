#ifndef E2A_EP_H
#define E2A_EP_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "Fiducial.h"
#include "Subtraction.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class e2a_ep_neutrino6_united4_radphot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::string ftarget;    // The target name  // ------------------------------->>>>>>>>>>>>>Mariana
   std::string fbeam_en;   // The beam energy  // ------------------------------->>>>>>>>>>>>>Mariana

   Fiducial   *fiducialcut;
   Subtraction *rotation;
   int fTorusCurrent;
   std::string target_name;
   std::map<std::string,double> en_beam;
   std::map<std::string,double> en_beam_Ecal;
   std::map<std::string,double> en_beam_Eqe;

   // Declaration of leaf types
   UChar_t         npart;
   UInt_t          runnb;
   UInt_t          evntid;
   UChar_t         evstat;
   Char_t          evntype;
   Int_t           evntclas;
   Float_t         q_l;
   Float_t         t_l;
   Float_t         tr_time;
   Float_t         rf_time;
   Int_t           l2bit;
   Int_t           l3bit;
   UChar_t         helicity;
   Int_t           hlsc;
   Int_t           intt;
   Int_t           helicity_cor;
   Int_t           gpart;  //Number of Particles in Event
   Int_t           id[40];   //[gpart]
   Int_t           stat[40];   //[gpart]
   Int_t           dc[40];   //[gpart]
   Int_t           cc[40];   //[gpart]
   Int_t           sc[40];   //[gpart]
   Int_t           ec[40];   //[gpart]
   Int_t           lec[40];   //[gpart]
   Int_t           st[40];   //[gpart]
   Float_t         p[40];   //[gpart]
   Float_t         m[40];   //[gpart]
   Int_t           q[40];   //[gpart]
   Float_t         b[40];   //[gpart]
   Float_t         cx[40];   //[gpart]
   Float_t         cy[40];   //[gpart]
   Float_t         cz[40];   //[gpart]
   Float_t         vx[40];   //[gpart]
   Float_t         vy[40];   //[gpart]
   Float_t         vz[40];   //[gpart]
   Int_t           dc_part;
   Int_t           dc_sect[40];   //[dc_part]
   Int_t           dc_trk[40];   //[dc_part]
   Int_t           dc_stat[40];   //[dc_part]
   Int_t           tb_st[40];   //[dc_part]
   Float_t         dc_xsc[40];   //[dc_part]
   Float_t         dc_ysc[40];   //[dc_part]
   Float_t         dc_zsc[40];   //[dc_part]
   Float_t         dc_cxsc[40];   //[dc_part]
   Float_t         dc_cysc[40];   //[dc_part]
   Float_t         dc_czsc[40];   //[dc_part]
   Float_t         dc_vx[40];   //[dc_part]
   Float_t         dc_vy[40];   //[dc_part]
   Float_t         dc_vz[40];   //[dc_part]
   Float_t         dc_vr[40];   //[dc_part]
   Float_t         tl1_cx[40];   //[dc_part]
   Float_t         tl1_cy[40];   //[dc_part]
   Float_t         tl1_cz[40];   //[dc_part]
   Float_t         tl1_x[40];   //[dc_part]
   Float_t         tl1_y[40];   //[dc_part]
   Float_t         tl1_z[40];   //[dc_part]
   Float_t         tl1_r[40];   //[dc_part]
   Float_t         dc_c2[40];   //[dc_part]
   Int_t           ec_part;
   Int_t           ec_stat[40];   //[ec_part]
   Int_t           ec_sect[40];   //[ec_part]
   Int_t           ec_whol[40];   //[ec_part]
   Int_t           ec_inst[40];   //[ec_part]
   Int_t           ec_oust[40];   //[ec_part]
   Float_t         etot[40];   //[ec_part]
   Float_t         ec_ei[40];   //[ec_part]
   Float_t         ec_eo[40];   //[ec_part]
   Float_t         ec_t[40];   //[ec_part]
   Float_t         ec_r[40];   //[ec_part]
   Float_t         ech_x[40];   //[ec_part]
   Float_t         ech_y[40];   //[ec_part]
   Float_t         ech_z[40];   //[ec_part]
   Float_t         ec_m2[40];   //[ec_part]
   Float_t         ec_m3[40];   //[ec_part]
   Float_t         ec_m4[40];   //[ec_part]
   Float_t         ec_c2[40];   //[ec_part]
   Int_t           sc_part;
   Int_t           sc_sect[40];   //[sc_part]
   Int_t           sc_hit[40];   //[sc_part]
   Int_t           sc_pd[40];   //[sc_part]
   Int_t           sc_stat[40];   //[sc_part]
   Float_t         edep[40];   //[sc_part]
   Float_t         sc_t[40];   //[sc_part]
   Float_t         sc_r[40];   //[sc_part]
   Float_t         sc_c2[40];   //[sc_part]
   Int_t           cc_part;
   Int_t           cc_sect[40];   //[cc_part]
   Int_t           cc_hit[40];   //[cc_part]
   Int_t           cc_segm[40];   //[cc_part]
   Int_t           nphe[40];   //[cc_part]
   Float_t         cc_t[40];   //[cc_part]
   Float_t         cc_r[40];   //[cc_part]
   Float_t         cc_c2[40];   //[cc_part]
   Int_t           lac_part;
   Int_t           lec_sect[40];   //[lac_part]
   Int_t           lec_hit[40];   //[lac_part]
   Int_t           lec_stat[40];   //[lac_part]
   Float_t         lec_etot[40];   //[lac_part]
   Float_t         lec_ein[40];   //[lac_part]
   Float_t         lec_t[40];   //[lac_part]
   Float_t         lec_r[40];   //[lac_part]
   Float_t         lec_x[40];   //[lac_part]
   Float_t         lec_y[40];   //[lac_part]
   Float_t         lec_z[40];   //[lac_part]
   Float_t         lec_c2[40];   //[lac_part]

   // List of branches
   TBranch        *b_npart;   //!
   TBranch        *b_runnb;   //!
   TBranch        *b_evntid;   //!
   TBranch        *b_evstat;   //!
   TBranch        *b_evntype;   //!
   TBranch        *b_evntclas;   //!
   TBranch        *b_q_l;   //!
   TBranch        *b_t_l;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_rf_time;   //!
   TBranch        *b_l2bit;   //!
   TBranch        *b_l3bit;   //!
   TBranch        *b_helicity;   //!
   TBranch        *b_hlsc;   //!
   TBranch        *b_intt;   //!
   TBranch        *b_helicity_cor;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_id;   //!
   TBranch        *b_stat;   //!
   TBranch        *b_dc;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_ec;   //!
   TBranch        *b_lec;   //!
   TBranch        *b_st;   //!
   TBranch        *b_p;   //!
   TBranch        *b_m;   //!
   TBranch        *b_q;   //!
   TBranch        *b_b;   //!
   TBranch        *b_cx;   //!
   TBranch        *b_cy;   //!
   TBranch        *b_cz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_dc_part;   //!
   TBranch        *b_dc_sect;   //!
   TBranch        *b_dc_trk;   //!
   TBranch        *b_dc_stat;   //!
   TBranch        *b_tb_st;   //!
   TBranch        *b_dc_xsc;   //!
   TBranch        *b_dc_ysc;   //!
   TBranch        *b_dc_zsc;   //!
   TBranch        *b_dc_cxsc;   //!
   TBranch        *b_dc_cysc;   //!
   TBranch        *b_dc_czsc;   //!
   TBranch        *b_dc_vx;   //!
   TBranch        *b_dc_vy;   //!
   TBranch        *b_dc_vz;   //!
   TBranch        *b_dc_vr;   //!
   TBranch        *b_tl1_cx;   //!
   TBranch        *b_tl1_cy;   //!
   TBranch        *b_tl1_cz;   //!
   TBranch        *b_tl1_x;   //!
   TBranch        *b_tl1_y;   //!
   TBranch        *b_tl1_z;   //!
   TBranch        *b_tl1_r;   //!
   TBranch        *b_dc_c2;   //!
   TBranch        *b_ec_part;   //!
   TBranch        *b_ec_stat;   //!
   TBranch        *b_ec_sect;   //!
   TBranch        *b_ec_whol;   //!
   TBranch        *b_ec_inst;   //!
   TBranch        *b_ec_oust;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_ec_ei;   //!
   TBranch        *b_ec_eo;   //!
   TBranch        *b_ec_t;   //!
   TBranch        *b_ec_r;   //!
   TBranch        *b_ech_x;   //!
   TBranch        *b_ech_y;   //!
   TBranch        *b_ech_z;   //!
   TBranch        *b_ec_m2;   //!
   TBranch        *b_ec_m3;   //!
   TBranch        *b_ec_m4;   //!
   TBranch        *b_ec_c2;   //!
   TBranch        *b_sc_part;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_hit;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_sc_stat;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_sc_c2;   //!
   TBranch        *b_cc_part;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_hit;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_nphe;   //!
   TBranch        *b_cc_t;   //!
   TBranch        *b_cc_r;   //!
   TBranch        *b_cc_c2;   //!
   TBranch        *b_lac_part;   //!
   TBranch        *b_lec_sect;   //!
   TBranch        *b_lec_hit;   //!
   TBranch        *b_lec_stat;   //!
   TBranch        *b_lec_etot;   //!
   TBranch        *b_lec_ein;   //!
   TBranch        *b_lec_t;   //!
   TBranch        *b_lec_r;   //!
   TBranch        *b_lec_x;   //!
   TBranch        *b_lec_y;   //!
   TBranch        *b_lec_z;   //!
   TBranch        *b_lec_c2;   //!

   e2a_ep_neutrino6_united4_radphot(std::string, std::string,TTree *tree=0);
   virtual ~e2a_ep_neutrino6_united4_radphot();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   Bool_t EFiducialCut(std::string beam_en, TVector3 momentum) {
     return fiducialcut->EFiducialCut(beam_en, momentum);
   }
   Bool_t PFiducialCut(std::string beam_en, TVector3 momentum) {
     return fiducialcut->PFiducialCut(beam_en, momentum);
   }
   Bool_t PiplFiducialCut(std::string beam_en, TVector3 momentum,Float_t *philow,Float_t *phiup) {
     return fiducialcut->PiplFiducialCut(beam_en, momentum, philow, phiup);
   }
   Bool_t PimiFiducialCut(std::string beam_en, TVector3 momentum,Float_t *philow,Float_t *phiup) {
     return fiducialcut->PimiFiducialCut(beam_en, momentum, philow, phiup);
   }
   bool Phot_fid(TVector3 V3_phot) {
     return fiducialcut->Phot_fid(V3_phot);
   }
   bool Pi_phot_fid_united(std::string beam_en, TVector3 V3_pi_phot, int q_pi_phot) {
     return fiducialcut->Pi_phot_fid_united(beam_en,V3_pi_phot, q_pi_phot);
   }
   Bool_t GetEPhiLimits(std::string beam_en, Float_t momentum, Float_t theta, Int_t sector, Float_t *EPhiMin, Float_t *EPhiMax) {
     return fiducialcut->GetEPhiLimits(beam_en, momentum, theta, sector, EPhiMin, EPhiMax);
   }
   void prot3_rot_func(std::string beam_en, TVector3 V3q, TVector3 V3prot[3],TVector3  V3prot_uncorr[3],TLorentzVector V4el,double Ecal_3pto2p[][2],double  pmiss_perp_3pto2p[][2],double  P3pto2p[][2],double N_p1[3],double Ecal_3pto1p[3],double  pmiss_perp_3pto1p[3], double *N_p3det);
   void prot2_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot[2],TVector3  V3prot_uncorr[2],TLorentzVector V4el,double Ecal_2pto1p[2],double  pmiss_perp_2pto1p[2],double  P2pto1p[2], double *Nboth);
   void prot1_pi1_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi, int q_pi,bool radstat, double *N_pi_p,double *N_nopi_p);
   void prot1_pi2_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi[2], int q_pi[2],bool radstat[2], double *P_1p0pi,double P_1p1pi[2]);
   void prot1_pi3_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi[3], int q_pi[3],bool radstat[3],double *P_tot);
   void prot2_pi1_rot_func(std::string beam_en, TVector3 V3q,TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_1pi, int q_pi,bool radstat,TLorentzVector V4_el, double Ecal_2p1pi_to2p0pi[2],double p_miss_perp_2p1pi_to2p0pi[2],double P_2p1pito2p0pi[2],double P_2p1pito1p1pi[2],double P_2p1pito1p0pi[2],double *P_tot);
   void prot2_pi2_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_2pi[2], int q_pi[2],bool radstat[2],TLorentzVector V4_el, double Ecal_2p2pi[2],double p_miss_perp_2p2pi[2],double P_tot_2p[2]);
   void prot3_pi1_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_3prot_corr[3],TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, int q_pi,bool radstat,TLorentzVector V4_el, double Ecal_3p1pi[3],double p_miss_perp_3p1pi[3],double P_tot_3p[3]);
   void pi1_rot_func(std::string beam_en, TVector3 V3_pi, int q_pi,bool radstat,TVector3 V3_q,double *P_pi);
   void pi2_rot_func(std::string beam_en, TVector3 V3_pi[2], int q_pi[2],bool radstat[2], TVector3 V3_q,double *P_0pi,double P_1pi[2]);
   void pi3_rot_func(std::string beam_en, TVector3 V3_pi[3], int q_pi[3],bool radstat[3], TVector3 V3_q,double *P_0pi,double P_1pi[3],double P_320[3],double P_3210[][2]);
   void pi4_rot_func(std::string beam_en, TVector3 V3_pi[4], int q_pi[4],bool radstat[4], TVector3 V3_q,double *P_0pi,double *P_410,double *P_420,double *P_4210,double *P_430,double *P_4310,double *P_4320,double *P_43210);


};

#endif
#ifdef E2A_EP_C

e2a_ep_neutrino6_united4_radphot::e2a_ep_neutrino6_united4_radphot(std::string a_target,std::string a_beam_en, TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

     fiducialcut = new Fiducial();
     rotation = new Subtraction();
     ftarget = a_target;
     fbeam_en= a_beam_en;
     fTorusCurrent = 0;



#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("ch",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ch","e2a_ep_neutrino6_united4_radphot");
      chain->Add(Form("/work/clas/clase2/Mariana/data/e2a_%s_%s_v1/*.root/h10", ftarget.c_str(), fbeam_en.c_str()));
      //chain->Add("datafile.root/h10");

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

e2a_ep_neutrino6_united4_radphot::~e2a_ep_neutrino6_united4_radphot()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t e2a_ep_neutrino6_united4_radphot::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t e2a_ep_neutrino6_united4_radphot::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void e2a_ep_neutrino6_united4_radphot::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("runnb", &runnb, &b_runnb);
   fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
   fChain->SetBranchAddress("evstat", &evstat, &b_evstat);
   fChain->SetBranchAddress("evntype", &evntype, &b_evntype);
   fChain->SetBranchAddress("evntclas", &evntclas, &b_evntclas);
   fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
   fChain->SetBranchAddress("t_l", &t_l, &b_t_l);
   fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
   fChain->SetBranchAddress("rf_time", &rf_time, &b_rf_time);
   fChain->SetBranchAddress("l2bit", &l2bit, &b_l2bit);
   fChain->SetBranchAddress("l3bit", &l3bit, &b_l3bit);
   fChain->SetBranchAddress("helicity", &helicity, &b_helicity);
   fChain->SetBranchAddress("hlsc", &hlsc, &b_hlsc);
   fChain->SetBranchAddress("intt", &intt, &b_intt);
   fChain->SetBranchAddress("helicity_cor", &helicity_cor, &b_helicity_cor);
   fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("stat", stat, &b_stat);
   fChain->SetBranchAddress("dc", dc, &b_dc);
   fChain->SetBranchAddress("cc", cc, &b_cc);
   fChain->SetBranchAddress("sc", sc, &b_sc);
   fChain->SetBranchAddress("ec", ec, &b_ec);
   fChain->SetBranchAddress("lec", lec, &b_lec);
   fChain->SetBranchAddress("st", st, &b_st);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("b", b, &b_b);
   fChain->SetBranchAddress("cx", cx, &b_cx);
   fChain->SetBranchAddress("cy", cy, &b_cy);
   fChain->SetBranchAddress("cz", cz, &b_cz);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
   fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
   fChain->SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
   fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
   fChain->SetBranchAddress("tb_st", tb_st, &b_tb_st);
   fChain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
   fChain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
   fChain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
   fChain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
   fChain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
   fChain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
   fChain->SetBranchAddress("dc_vx", dc_vx, &b_dc_vx);
   fChain->SetBranchAddress("dc_vy", dc_vy, &b_dc_vy);
   fChain->SetBranchAddress("dc_vz", dc_vz, &b_dc_vz);
   fChain->SetBranchAddress("dc_vr", dc_vr, &b_dc_vr);
   fChain->SetBranchAddress("tl1_cx", tl1_cx, &b_tl1_cx);
   fChain->SetBranchAddress("tl1_cy", tl1_cy, &b_tl1_cy);
   fChain->SetBranchAddress("tl1_cz", tl1_cz, &b_tl1_cz);
   fChain->SetBranchAddress("tl1_x", tl1_x, &b_tl1_x);
   fChain->SetBranchAddress("tl1_y", tl1_y, &b_tl1_y);
   fChain->SetBranchAddress("tl1_z", tl1_z, &b_tl1_z);
   fChain->SetBranchAddress("tl1_r", tl1_r, &b_tl1_r);
   fChain->SetBranchAddress("dc_c2", dc_c2, &b_dc_c2);
   fChain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
   fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
   fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
   fChain->SetBranchAddress("ec_whol", ec_whol, &b_ec_whol);
   fChain->SetBranchAddress("ec_inst", ec_inst, &b_ec_inst);
   fChain->SetBranchAddress("ec_oust", ec_oust, &b_ec_oust);
   fChain->SetBranchAddress("etot", etot, &b_etot);
   fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
   fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
   fChain->SetBranchAddress("ec_t", ec_t, &b_ec_t);
   fChain->SetBranchAddress("ec_r", ec_r, &b_ec_r);
   fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
   fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
   fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
   fChain->SetBranchAddress("ec_m2", ec_m2, &b_ec_m2);
   fChain->SetBranchAddress("ec_m3", ec_m3, &b_ec_m3);
   fChain->SetBranchAddress("ec_m4", ec_m4, &b_ec_m4);
   fChain->SetBranchAddress("ec_c2", ec_c2, &b_ec_c2);
   fChain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
   fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
   fChain->SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
   fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
   fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
   fChain->SetBranchAddress("edep", edep, &b_edep);
   fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
   fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
   fChain->SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
   fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
   fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
   fChain->SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
   fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
   fChain->SetBranchAddress("nphe", nphe, &b_nphe);
   fChain->SetBranchAddress("cc_t", cc_t, &b_cc_t);
   fChain->SetBranchAddress("cc_r", cc_r, &b_cc_r);
   fChain->SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
   fChain->SetBranchAddress("lac_part", &lac_part, &b_lac_part);
   fChain->SetBranchAddress("lec_sect", lec_sect, &b_lec_sect);
   fChain->SetBranchAddress("lec_hit", lec_hit, &b_lec_hit);
   fChain->SetBranchAddress("lec_stat", lec_stat, &b_lec_stat);
   fChain->SetBranchAddress("lec_etot", lec_etot, &b_lec_etot);
   fChain->SetBranchAddress("lec_ein", lec_ein, &b_lec_ein);
   fChain->SetBranchAddress("lec_t", lec_t, &b_lec_t);
   fChain->SetBranchAddress("lec_r", lec_r, &b_lec_r);
   fChain->SetBranchAddress("lec_x", lec_x, &b_lec_x);
   fChain->SetBranchAddress("lec_y", lec_y, &b_lec_y);
   fChain->SetBranchAddress("lec_z", lec_z, &b_lec_z);
   fChain->SetBranchAddress("lec_c2", lec_c2, &b_lec_c2);
   Notify();
}

Bool_t e2a_ep_neutrino6_united4_radphot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void e2a_ep_neutrino6_united4_radphot::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t e2a_ep_neutrino6_united4_radphot::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef e2a_ep_neutrino6_united4_radphot_h
