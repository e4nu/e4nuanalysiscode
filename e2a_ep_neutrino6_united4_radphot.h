//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 12 15:09:51 2019 by ROOT version 5.34/36
// from TChain ch/e2a_ep_neutrino6_united4_radphot
//////////////////////////////////////////////////////////

#ifndef e2a_ep_neutrino6_united4_radphot_h
#define e2a_ep_neutrino6_united4_radphot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class e2a_ep_neutrino6_united4_radphot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

  std::string ftarget;    // The target name  // ------------------------------->>>>>>>>>>>>>Mariana
   std::string fbeam_en;   // The beam energy  // ------------------------------->>>>>>>>>>>>>Mariana

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
   Int_t           gpart;
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
};

#endif

#ifdef e2a_ep_neutrino6_united4_radphot_cxx
e2a_ep_neutrino6_united4_radphot::e2a_ep_neutrino6_united4_radphot(std::string a_target,std::string a_beam_en,TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {


     ftarget = a_target;
     fbeam_en=a_beam_en;


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
#endif // #ifdef e2a_ep_neutrino6_united4_radphot_cxx




const Float_t par_EcUVW[6][3] = {{60, 360, 400}, {55, 360, 400}, {50, 363, 400}, {52, 365, 396}, {60, 360, 398}, {50, 362, 398}};

const Float_t fgPar_2GeV_2250_Efid[6][6][9] = {{{62.2935, -92.5133, 87.0360, -38.4696, 6.3177, 0, 0, 0, 0}, {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0}, {-140.845, 1381.30, -4499.99, 7557.27, -7140.27, 3828.75, -1086.21, 126.468, 0}, {497.951, -1846.42, 2759.58, -1634.71, 345.006, 0, 0, 0, 0}, {9.40986, 180.752, -646.771, 1055.14, -909.094, 424.435, -99.8368, 9.02086, 0}, {288.485, -1016.03, 1463.72, -859.231, 185.976, 0, 0, 0, 0}}, {{61.1474, -88.768, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0}, {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0}, {21.3087, 138.975, -672.710, 1324.20, -1326.12, 714.866, -197.531, 21.9144, 0}, {375.091, -1411.50, 2082.58, -1192.17, 239.685, 0, 0, 0, 0}, {-121.816, 1182.59, -3800.98, 6319.82, -5937.33, 3179.37, -903.954, 105.764, 0}, {-4781.96, 43165.9, -159567, 318502, -376469, 271207, -116893, 27698.9, -2775.61}}, {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0}, {73.7620, -34.6321, -41.8796, 117.543, -81.2043, 17.1718, 0, 0, 0}, {157.046, -765.472, 1735.21, -2053.86, 1371.34, -515.214, 101.081, -8.07402, 0}, {-608.740, 4827.18, -13239.6, 17742.4, -12420.0, 4369.11, -607.877, 0, 0}, {-274.278, 2380.63, -7560.19, 12582.3, -11924.5, 6464.66, -1863.44, 221.134, 0}, {-1240.72, 8096.04, -19407.0, 23942.9, -16052.3, 5559.32, -776.123, 0, 0}}, {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0}, {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0}, {-71.2528, 879.668, -3027.37, 5226.61, -4999.19, 2689.35, -761.206, 88.1242, 0}, {-1269.89, 9486.25, -26103.8, 35581.2, -25373.0, 9062.87, -1277.60, 0, 0}, {-186.640, 1811.85, -6032.01, 10283.3, -9808.11, 5285.35, -1501.87, 174.799, 0}, {-530.826, 4643.56, -13864.2, 20580.2, -15898.0, 6106.69, -916.365, 0, 0}}, {{61.6665, -90.4268, 84.5606, -37.2240, 6.09207, 0, 0, 0, 0}, {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0}, {-1.53910, 216.936, -701.057, 1167.26, -1111.92, 615.364, -183.854, 22.8595, 0}, {-19.7415, 454.317, -1250.51, 1512.52, -762.408, 137.695, 0, 0, 0}, {-55.9612, 657.449, -2049.73, 3295.30, -2995.85, 1553.68, -427.764, 48.4324, 0}, {-522.682, 3356.77, -7535.50, 8756.49, -5518.61, 1795.60, -235.144, 0, 0}}, {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0}, {73.7620, -34.6321, -41.8796, 117.543, -81.2043, 17.1718, 0, 0, 0}, {-82.0368, 883.261, -2828.84, 4621.53, -4223.56, 2185.52, -598.218, 67.2908, 0}, {608.323, -2743.56, 4942.01, -4045.58, 1558.07, -226.240, 0, 0, 0}, {4.07203, 138.882, -321.983, 282.702, -12.9566, -129.159, 74.5884, -12.9994, 0}, {-866.737, 5984.13, -15129.6, 19134.6, -12757.7, 4276.79, -566.056, 0, 0}}};

const Float_t fgPar_2GeV_2250_EfidTheta_S3[4][8] = {{74.4893, -158.720, 251.241, -200.000, 52.2984, 25.4188, -18.8692, 3.27217}, {90.8413, -226.800, 358.487, -259.260, 30.9359, 68.7248, -38.9760, 6.47933}, {117.102, -429.455, 1208.29, -1922.72, 1791.40, -965.135, 277.459, -32.8536}, {55.0676, 91.0959, -444.252, 791.284, -717.492, 350.325, -87.3235, 8.68087}};

const Float_t fgPar_2GeV_2250_EfidTheta_S4[2][8] = {{77.7940, -192.492, 361.852, -394.127, 246.499, -84.6133, 13.9182, -0.713846}, {41.2902, 110.603, -586.690, 1130.70, -1137.27, 633.345, -185.038, 22.1482}};

const Float_t fgPar_2GeV_2250_EfidTheta_S5[8][8] = {{-12998.3, 57694.0, -109085, 114102, -71303.6, 26616.9, -5494.45, 483.756}, {17842.9, -74659.9, 133869, -133170, 79380.0, -28352.3, 5617.88, -476.314}, {65.5364, -99.8689, 88.5645, 24.8299, -121.327, 102.818, -37.8275, 5.28492}, {66.4049, -76.4096, -20.8674, 230.072, -318.905, 206.721, -66.3286, 8.48753}, {100.262, -358.882, 957.267, -1495.42, 1396.73, -765.881, 226.791, -27.9341}, {50.4447, 48.3032, -315.976, 580.141, -525.583, 252.075, -59.9294, 5.34805}, {78.5845, -155.728, 320.528, -420.296, 341.899, -164.626, 42.5274, -4.50224}, {95.9430, -221.787, 391.495, -350.033, 131.391, 13.2965, -24.0460, 4.92253}};

const Float_t fgPar_2GeV_2250_Pfid_For[6][4][7] = 
{{{60.2165, -189.720, 446.990, -523.122, 320.721, -97.8518, 11.5258}, 
  {-1457.16, 13814.2, -43182.7, 66646.0, -54355.1, 22423.5, -3683.76}, 
  {17.1086, 54.2974, -103.464, 111.325, -70.7673, 27.2551, -5.02858},
  {-2547.86, 22143.1, -66326.6, 101105.0, -82187.8, 33959.7, -5607.59}}, 
 {{65.7242, -246.922, 759.745, -1198.32, 1007.05, -428.060, 72.2644}, 
  {3384.16, -19353.1, 54083.5, -79843.4, 63870.2, -26079.2, 4250.29}, 
  {85.2489, -441.821, 1327.52, -1978.53, 1567.84, -633.530, 102.928}, 
  {411.998, -533.572, 599.925, 2099.52, -5061.48, 3701.58, -891.843}},
 {{110.022, -558.044, 1512.96, -2098.53, 1579.55, -613.478, 96.3279}, 
  {3937.29, -23745.1, 59651.0, -76988.6, 54276.0, -19900.2, 2974.95}, 
  {35.8488, -46.9595, 107.492, -93.9141, 10.5845, 26.1910, -9.89460}, 
  {-326.838, 4634.99, -11155.2, 11811.4, -5405.80, 554.030, 175.526}}, 
 {{38.9338, -62.8663, 118.218, -56.6953, -40.5083, 46.1782, -11.5822}, 
  {1864.83, -11735.6, 34175.4, -48928.5, 37315.8, -14496.1, 2254.05}, 
  {23.6892, 9.69854, 94.4521, -270.119, 288.132, -140.031, 25.9272}, 
  {-261.086, 4863.13, -11760.4, 13791.1, -8983.19, 3136.52, -457.183}}, 
 {{-11.0252, 348.901, -1172.63, 1980.73, -1759.08, 786.043, -139.299}, 
  {-2231.41, 23477.1, -78229.3, 129238.0, -111761.0, 48561.4, -8370.65}, 
  {104.415, -548.464, 1506.70, -2064.10, 1507.55, -561.677, 83.9247}, 
  {1402.87, -9008.78, 25660.0, -37543.3, 29860.8, -12238.4, 2019.03}}, 
 {{20.4577, 66.1373, -205.218, 372.864, -366.625, 177.596, -33.1168}, 
  {2059.77, -14468.3, 46492.9, -72168.2, 58275.9, -23615.8, 3800.60}, 
  {-18.9897, 392.519, -1234.31, 1950.24, -1623.01, 681.260, -113.806},
  {-3478.50, 32840.9, -104381.0, 167656.0, -143070.0, 61909.3, -10690.1}}};

const Float_t fgPar_2GeV_2250_Pfid_Bak[6][4][7] =
{{{110.007, 121.302, 97.8380, -1679.71, 4022.73, -3973.09, 1422.42}, 
  {69.7305, 359.843, -876.383, 649.612, 600.059, -1155.43, 472.866}, 
  {13.9334, -236.587, 810.783, -1614.65, 1851.97, -1125.48, 280.069}, 
  {10.1644, 51.7943, -527.843, 2071.12, -3480.34, 2663.52, -768.498}}, 
 {{161.555, -263.801, 770.924, -902.814, 503.641, -319.619, 171.147},
  {154.660, -619.711, 3444.65, -8994.29, 12253.9, -8439.82, 2321.14}, 
  {117.461, -1429.96, 6117.79, -13492.3, 16142.2, -9965.40, 2490.47}, 
  {7.77411, -17.3501, 279.462, -876.326, 1398.82, -1137.49, 365.383}}, 
 {{-31.1460, 1942.49, -9193.97, 21731.0, -26961.3, 16701.7, -4067.85}, 
  {154.660, -654.420, 3774.08, -9920.36, 13333.7, -8953.68, 2386.32}, 
  {63.2709, -867.859, 4000.97, -9557.57, 12215.1, -7926.91, 2052.90}, 
  {-28.1127, 484.636, -2665.71, 7484.94, -10740.7, 7561.79, -2076.70}}, 
 {{172.853, -656.312, 3768.76, -10243.0, 14600.3, -10616.3, 3095.27}, 
  {270.076, -1938.46, 9276.01, -21861.1, 27363.7, -17479.9, 4490.05}, 
  {32.2327, -432.593, 1666.57, -3491.43, 4031.58, -2406.30, 579.944},
  {-44.9153, 638.112, -2971.77, 7223.13, -9328.99, 6080.46, -1576.13}}, 
 {{45.7403, 875.133, -3646.85, 7848.52, -8905.36, 4914.78, -1010.91}, 
  {138.000, -449.485, 2806.13, -7725.44, 10777.3, -7482.95, 2056.80}, 
  {72.7551, -944.002, 4200.92, -9776.76, 12316.6, -7955.78, 2066.50},
  {-9.59531, 180.519, -795.797, 2124.85, -2978.29, 2040.14, -541.811}}, 
 {{77.5100, 494.571, -1625.99, 2397.48, -1177.99, -574.604, 530.446}, 
  {117.869, -56.8761, 330.252, -715.276, 807.257, -497.124, 133.989}, 
  {7.66164, -208.001, 996.883, -2772.33, 4100.81, -3008.90, 864.126}, 
  {-25.3497, 346.501, -1458.46, 3513.62, -4625.70, 3088.01, -818.696}}};

const Float_t fgPar_2GeV_2250_Pfid_ScpdS2[2][6] = 
{{-28.1486, 425.124, -935.693, 1065.39, -608.526, 137.658}, 
 {-15.2084, 345.466, -697.657, 751.738, -419.288, 95.2206}};

const Float_t fgPar_2GeV_2250_Pfid_ScpdS3[8][6] = 
{{17.1490, 294.605, -640.590, 707.758, -386.730, 83.2529},
 {35.9318, 204.580, -404.489, 413.240, -209.580, 41.7819}, 
 {47.6825, 274.777, -754.725, 1117.80, -846.816, 255.607}, 
 {44.7484, 344.543, -872.200, 1113.89, -694.736, 168.061}, 
 {-205.978, 828.617, -1199.65, 875.482, -317.846, 45.6938},
 {-240.595, 961.068, -1370.34, 977.625, -345.743, 48.3834},
 {-136.104, 479.276, -593.135, 374.730, -118.350, 14.7923},
 {-196.773, 700.974, -894.540, 577.460, -185.690, 23.6201}};

const Float_t fgPar_2GeV_2250_Pfid_ScpdS4[4][6] = 
{{81.8115, 139.810, -445.130, 804.212, -821.194, 364.924},
 {79.5053, 317.287, -1582.80, 3987.05, -4880.55, 2305.63}, 
 {-137.480, 633.288, -954.383, 721.057, -269.140, 39.4822}, 
 {-145.605, 697.662, -1088.74, 853.855, -330.883, 50.3421}};

const Float_t fgPar_2GeV_2250_Pfid_ScpdS5[8][6] =
{{-29.9426, 370.963, -714.697, 707.343, -348.995, 67.7647}, 
 {-27.4173, 372.536, -693.341, 652.792, -302.559, 54.7761}, 
 {-47.1617, 132.967, -104.776, 41.7673, -7.68238, 0.404311}, 
 {-54.5895, 149.685, -111.590, 41.2556, -6.93943, 0.301087}, 
 {-79.1386, 275.678, -341.972, 218.907, -69.5520, 8.66381}, 
 {-97.5794, 352.616, -468.487, 322.829, -111.159, 15.0975}, 
 {22.5823, -182.064, 365.317, -294.653, 108.779, -15.2712}, 
 {-7.59521, 2.91795, 31.6773, -28.3085, 10.5943, -1.57966}};




/*
//Parameters for 4 GeV proton's Fiducial Cut Rustam Niyazov
// <A HREF="http://www.physics.odu.edu/~rust/clas/fidp.html"</A> --Rustam Niyazov (ODU).

const Float_t fgPar_4Gev_2250_Pfidft1l[6][6]={
  {26.2564,0.441269,-29.7632,94.5137,7.71903,2.10915},
  {29.7455,-0.826489,4.09596,91.8187,8.38108,1.5016},
  {29.5399,-0.878321,43.1909,64.9772,11.1844,0.825411},
  {28.5857,0.4061,98.6296,95.5022,13.7297,0.415071},
  {31.9803,0.341766,257.124,103.504,14.2357,0.43387},
  {29.2846,-0.257616,51.1709,84.3207,10.2963,1.69991}};
const Float_t fgPar_4Gev_2250_Pfidft1r[6][6]={
  {34.7359,-1.45301,660.653,-79.1375,11.3239,1.05352},
  {30.6992,0.71858,442.087,4.20897,3.62722,3.35155},
  {19.1518,3.71404,-197.134,177.828,9.63173,1.35402},
  {23.9897,1.52101,23.9288,71.4476,8.89464,1.69512},
  {22.6619,2.4697,-54.5174,112.22,11.2561,0.687839},
  {20.9859,3.86504,-56.5229,230.635,13.6587,0.270987}}; 
const Float_t fgPar_4Gev_2250_Pfidft2l[6][6]={
  {24.683,0.470268,124.501,-9.04329,8.60129,1.66063},
  {26.2736,-0.591497,182.954,-51.059,7.65701,2.29757},
  {24.8681,1.15526,111.322,22.2304,9.46319,1.6834},
  {29.3639,1.307,282.797,89.5863,11.7162,0.376266},
  {36.8099,-0.785452,655.368,46.4935,12.0443,0.500522},
  {25.8401,0.899645,141.723,27.6687,9.62103,1.7379}}; 
const Float_t fgPar_4Gev_2250_Pfidft2r[6][6]={
  {32.9905,-0.580968,464.263,30.5379,11.7414,0.320415},
  {26.8867,0.748481,150.349,51.4182,8.70942,1.51013},
  {26.0729,0.357197,136.456,24.1839,6.70568,0.820883},
  {25.8339,1.018,149.648,38.7987,6.56928,0.527773},
  {27.997,0.0685368,268.87,-45.3343,5.26386,3.08026},
  {30.3568,1.60206,359.39,197.047,11.1523,0.451219}}; 
const Float_t fgPar_4Gev_2250_Pfidbt1l[6][6]= {
  {-24.4118,4.20154,-0.0480933,-0.0800641,0.000311929,0.000511191},
  {-34.5523,8.81812,0.221281,-0.203846,-0.00115322,0.00119883},
  {-29.4962,6.57417,0.0830637,-0.142094,-0.000271087,0.000801481},
  {-29.5177,6.23458,0.183415,-0.160458,-0.00121912,0.0010282},
  {-19.8091,4.37431,-0.046672,-0.124147,-7.21454e-05,0.000931229},
  {-38.1865,10.6462,0.363126,-0.267793,-0.00212252,0.00162732}}; 
const Float_t fgPar_4Gev_2250_Pfidbt1r[6][6]={
  {-15.6987,3.34818,-0.155291,-0.102923,0.000736214,0.000775517},
  {-15.9442,1.75807,-0.196246,-0.0524198,0.00118102,0.000398854},
  {-14.4453,1.65733,-0.269699,-0.0423913,0.00187485,0.000274252},
  {-18.5972,1.41622,-0.144491,-0.0369631,0.000874762,0.000326006},
  {-17.1008,0.577868,-0.173353,-0.021315,0.00108238,0.000189545},
  {2.21904,-3.38706,-0.636698,0.0953525,0.0038789,-0.000559086}}; 
const Float_t fgPar_4Gev_2250_Pfidbt2l[6][6]={
  {-13.7253,-1.53789,-0.296133,0.0648705,0.00269427,-0.000928492},
  {-12.356,-2.62192,-0.366191,0.115155,0.0033624,-0.00137599},
  {-2.52638,-9.6591,-0.743505,0.380195,0.0067055,-0.00369404},
  {-34.5804,15.3815,0.417723,-0.489802,-0.00337546,0.00370894},
  {1.87747,-7.70598,-0.919924,0.376373,0.00776553,-0.00354661},
  {-12.3968,-2.37408,-0.367352,0.114661,0.00352523,-0.00148841}}; 
const Float_t fgPar_4Gev_2250_Pfidbt2r[6][6]={
  {-29.5895,10.9088,0.248994,-0.326966,-0.00154954,0.00202508},
  {-7.20087,-6.19132,-0.568426,0.257971,0.00476513,-0.00236084},
  {-10.0076,-3.66545,-0.468027,0.163446,0.00421363,-0.00175242},
  {-9.03582,-5.14009,-0.515592,0.221044,0.00482855,-0.00237549},
  {-8.55955,-5.27785,-0.504058,0.201472,0.00404296,-0.00175892},
  {-21.122,5.19264,-0.0761427,-0.0826774,0.0018747,-0.000390706}}; 
const Float_t fgPar_4Gev_2250_Pfidbl[6][6]={
  {131.839,-6.64199,-22.8623,4.91185,126.5,20},
  {132.055,-5.2283,2.20945,-1.57951,128.429,11.4286},
  {137.945,-7.90553,-12.8716,3.94534,119.857,22.8571},
  {124.743,-3.54503,-22.8263,5.62231,130.429,11.4286},
  {136.455,-7.59559,-18.6847,4.52149,123.5,20},
  {126.556,-4.02284,-22.2328,5.23298,124.857,22.8571}};
const Float_t fgPar_4Gev_2250_Pfidbr[6][6]={
  {97.3917,2.99764,26.7715,-5.95695,126.5,20},
  {132.154,-6.60261,0.000146616,1.53542,128.429,11.4286},
  {113.746,-1.24667,32.0728,-9.35241,119.857,22.8571},
  {118.596,-2.44983,22.2973,-5.40976,130.429,11.4286},
  {125.129,-3.96273,21.6178,-5.86908,123.5,20},
  {111.201,-0.178015,25.1267,-6.55928,124.857,22.8571}}; 
*/




/*

// Steven McLauchlan (1/30/2001)
// Parameters for 1.1GeV fiducial cut
// 1500A torus current

const Double_t fgPar_1Gev_1500_Efid[6][5][6] = 
{{{-2.550493918887335,0.389508196279333,-0.001333283523319,0.000001947815908,-0.000000001316913,0.000000000000338},
  {739.665966640729152,-7.564586257707806,0.029070457076481,-0.000049407992830,0.000000038875378,-0.000000000011556},
  {-38.943294334931004,0.604798911243653,-0.002123924852338,0.000003596540690,-0.000000002945430,0.000000000000935},
  {1674.135084781084061,-15.638896622395071,0.056064673823181,-0.000094802481347,0.000000077287133,-0.000000000024356},
  {36.467163175325709,-0.149673554778768,0.000640450951938,-0.000001217358940,0.000000001070681,-0.000000000000351}},
 {{99.490760702509277,-0.454366869265921,0.001261304562511,-0.000001848741434,0.000000001351091,-0.000000000000387},
  {532.589474908597253,-1.974841260959739,0.007824354020230,-0.000016500091869,0.000000015314082,-0.000000000004989},
  {290.208151444308214,-1.855305261275745,0.005273424450297,-0.000007396122888,0.000000005060998,-0.000000000001343},
  {818.131759192972595,-6.630367625399796,0.020821102579573,-0.000027931394648,0.000000016344909,-0.000000000003256},
  {312.687789127750818,-2.178522834336562,0.006427983040174,-0.000009171395143,0.000000006328733,-0.000000000001690}},
 {{170.498641405867033,-0.966479071879565,0.002760453538862,-0.000003980182820,0.000000002815775,-0.000000000000777},
  {466.650046236903790,-3.304087162850300,0.008149251774215,-0.000006095026150,-0.000000001120372,0.000000000002095},
  {143.871532905251740,-0.916283465837924,0.002719612026644,-0.000003883439552,0.000000002675184,-0.000000000000709},
  {-2125.925066263838289,18.997533045118953,-0.064480904152770,0.000106234237021,-0.000000083803115,0.000000000025329},
  {41.722117894883461,-0.051366641486730,-0.000069860488095,0.000000472584061,-0.000000000621081,0.000000000000257}},
 {{171.754526904421567,-0.991058411485744,0.002787362029463,-0.000003945823511,0.000000002746342,-0.000000000000747},
  {-555.198014601586806,10.477079618697317,-0.037763564878876,0.000062643558690,-0.000000049665924,0.000000000015182},
  {288.616413125497161,-1.811215191482553,0.005136945426641,-0.000007195297301,0.000000004924735,-0.000000000001312},
  {-4643.290039163491201,38.734934810610625,-0.122303938566621,0.000188124697366,-0.000000140071715,0.000000000040371},
  {9.004203411705985,0.248267693775404,-0.001022198891318,0.000001864956947,-0.000000001576903,0.000000000000505}},
 {{42.952481198425112,0.023500888648231,-0.000278691000758,0.000000537010290,-0.000000000429457,0.000000000000126},
  {218.648760061770162,-2.778128633827438,0.012527499672852,-0.000022772919229,0.000000018984515,-0.000000000006007},
  {5.023697351767844,0.156293246543284,-0.000418113378155,0.000000518527417,-0.000000000294285,0.000000000000060},
  {-1006.509682268774100,9.875176831851448,-0.035000976170258,0.000058770106071,-0.000000046295933,0.000000000013793},
  {53.421679849517517,-0.120820531994036,0.000058685352898,0.000000344132876,-0.000000000514082,0.000000000000205}},
 {{79.766176403990343,-0.278810941289375,0.000645077454885,-0.000000784617681,0.000000000465126,-0.000000000000105},
  {3112.109016419252384,-25.096095304234694,0.081240967418761,-0.000127812099473,0.000000098172253,-0.000000000029353},
  {255.639083801755930,-1.655317133465159,0.004716500230189,-0.000006597502724,0.000000004522531,-0.000000000001211},
  {-1080.829142649767164,11.103472805372817,-0.040432491192383,0.000069496015412,-0.000000055728653,0.000000000016820},
  {212.220246941877747,-1.327983077827758,0.003603350662534,-0.000004720527866,0.000000003013837,-0.000000000000755}}};


// Steven McLauchlan (1/30/2001)
// Parameters for 1.1GeV fiducial cut
// 750A torus current
const Double_t fgPar_1Gev_750_Efid[6][5][6] = 
{{{53.448466137147108,-0.086287764363944,-0.000095208568679,0.000000573157491,-0.000000000693952,0.000000000000267},
  {2395.572627329481747,-19.579668396704726,0.064797456172212,-0.000100845569270,0.000000074006368,-0.000000000020579},
  {259.560198858548574,-1.844516389508475,0.005648994718534,-0.000008298117560,0.000000005840166,-0.000000000001577},
  {-1587.804278895457401,12.696688014266641,-0.035969419740232,0.000048344042919,-0.000000030096499,0.000000000006884},
  {31.244744412870517,-0.025807389788636,-0.000011577835718,0.000000214273351,-0.000000000339440,0.000000000000158}},
 {{34.465276195094624,0.045937262105175,-0.000434902242746,0.000000932745139,-0.000000000833190,0.000000000000271},
  {3635.381526237591515,-22.391133082958913,0.063385989108405,-0.000088140825578,0.000000058919265,-0.000000000015056},
  {198.778521131712267,-1.138544385556973,0.003128878030819,-0.000004213574978,0.000000002729891,-0.000000000000672},
  {1840.285580954853003,-14.731530568584006,0.047634304990508,-0.000070778318684,0.000000049201832,-0.000000000012891},
  {239.703202430365707,-1.665091837883518,0.005104440345427,-0.000007626633722,0.000000005541237,-0.000000000001565}},
 {{62.648795387438646,-0.202971507270251,0.000412901309633,-0.000000426657627,0.000000000204401,-0.000000000000033},
  {5505.758250992919784,-44.246117648832730,0.142983640614620,-0.000217740239693,0.000000157627772,-0.000000000043738},
  {137.610864617097349,-0.877776486033249,0.002764189351627,-0.000004129583160,0.000000002927356,-0.000000000000789},
  {344.946357574946262,5.621617559675438,-0.020344778542695,0.000033499083642,-0.000000026928489,0.000000000008469},
  {-35.029472261198009,0.766223788807457,-0.002745396684495,0.000004621528350,-0.000000003749284,0.000000000001183}},
 {{109.141173435470662,-0.618600665002767,0.001816612190204,-0.000002714071430,0.000000002008202,-0.000000000000584},
  {176.654641794276699,0.980557784282395,-0.004434087390540,0.000008006260292,-0.000000005770392,0.000000000001306},
  {-18.993624563198331,0.402038850300937,-0.001205887423007,0.000001590733648,-0.000000000912114,0.000000000000169},
  {1175.291585819301872,-3.528334832328303,0.008977772742253,-0.000011419463583,0.000000005926540,-0.000000000000681},
  {223.097532523005469,-1.419016100626559,0.004258429941881,-0.000006247651659,0.000000004429123,-0.000000000001208}},
 {{99.175888767618616,-0.484455189099180,0.001277836120915,-0.000001741410194,0.000000001190341,-0.000000000000324},
  {2108.985526175658379,-18.749089325155833,0.066411589994549,-0.000110694155495,0.000000088167629,-0.000000000026950},
  {102.031283321769024,-0.608825833913801,0.001948334892035,-0.000003043150820,0.000000002292849,-0.000000000000661},
  {-2834.362182720959026,23.552345615317957,-0.073523328825804,0.000111121579335,-0.000000080561949,0.000000000022405},
  {-67.714437911752213,0.790236817449855,-0.002559147773978,0.000003978412216,-0.000000002962100,0.000000000000848}},
 {{54.075117651289297,-0.167277937624739,0.000416438439455,-0.000000645102074,0.000000000538254,-0.000000000000179},
  {1674.291250256091644,-8.510168235424377,0.024353652698827,-0.000035048237822,0.000000024024525,-0.000000000006206},
  {194.725459887803083,-1.227022197653817,0.003754095508600,-0.000005645875233,0.000000004107310,-0.000000000001150},
  {580.651869204728655,-1.601754520529919,-0.000886817864826,0.000009335680173,-0.000000011083844,0.000000000003972},
  {55.543883831755714,-0.182118215858001,0.000519579714073,-0.000000834302749,0.000000000683447,-0.000000000000216}}};

*/
