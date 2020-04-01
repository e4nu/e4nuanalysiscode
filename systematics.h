#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "Fiducial.h"
#include "Subtraction.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class systematics {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::string ftarget;    // The target name  // ------------------------------->>>>>>>>>>>>>Mariana
   std::string fbeam_en;   // The beam energy  // ------------------------------->>>>>>>>>>>>>Mariana

   int N_tot;
   Fiducial   *fiducialcut;
   Subtraction *rotation;
   int fTorusCurrent;
   int fchoice;
   std::string target_name;
   std::map<std::string,double> en_beam;
   std::map<std::string,double> en_beam_Ecal;
   std::map<std::string,double> en_beam_Eqe;

   // Declaration of leaf types
   Int_t           iev;
   Int_t           neu;
   Int_t           fspl;
   Int_t           tgt;
   Int_t           Z;
   Int_t           A;
   Int_t           hitnuc;
   Int_t           hitqrk;
   Int_t           resid;
   Bool_t          sea;
   Bool_t          qel;
   Bool_t          mec;
   Bool_t          res;
   Bool_t          dis;
   Bool_t          coh;
   Bool_t          dfr;
   Bool_t          imd;
   Bool_t          imdanh;
   Bool_t          singlek;
   Bool_t          nuel;
   Bool_t          em;
   Bool_t          cc;
   Bool_t          nc;
   Bool_t          charm;
   Int_t           neut_code;
   Int_t           nuance_code;
   Double_t        wght;
   Double_t        xs;
   Double_t        ys;
   Double_t        ts;
   Double_t        Q2s;
   Double_t        Ws;
   Double_t        x;
   Double_t        y;
   Double_t        t;
   Double_t        Q2;
   Double_t        W;
   Double_t        EvRF;
   Double_t        Ev;
   Double_t        pxv;
   Double_t        pyv;
   Double_t        pzv;
   Double_t        En;
   Double_t        pxn;
   Double_t        pyn;
   Double_t        pzn;
   Double_t        El;
   Double_t        pxl;
   Double_t        pyl;
   Double_t        pzl;
   Double_t        pl;
   Double_t        cthl;
   Int_t           nfp;
   Int_t           nfn;
   Int_t           nfpip;
   Int_t           nfpim;
   Int_t           nfpi0;
   Int_t           nfkp;
   Int_t           nfkm;
   Int_t           nfk0;
   Int_t           nfem;
   Int_t           nfother;
   Int_t           nip;
   Int_t           nin;
   Int_t           nipip;
   Int_t           nipim;
   Int_t           nipi0;
   Int_t           nikp;
   Int_t           nikm;
   Int_t           nik0;
   Int_t           niem;
   Int_t           niother;
   Int_t           ni;
   Int_t           pdgi[2];   //[ni]
   Int_t           resc[1];   //[ni]
   Double_t        Ei[2];   //[ni]
   Double_t        pxi[2];   //[ni]
   Double_t        pyi[2];   //[ni]
   Double_t        pzi[2];   //[ni]
   Int_t           nf;
   Int_t           pdgf[120];   //[nf]
   Double_t        Ef[120];   //[nf]
   Double_t        pxf[120];   //[nf]
   Double_t        pyf[120];   //[nf]
   Double_t        pzf[120];   //[nf]
   Double_t        pf[120];   //[nf]
   Double_t        cthf[120];   //[nf]
   Double_t        vtxx;
   Double_t        vtxy;
   Double_t        vtxz;
   Double_t        vtxt;
   Double_t        sumKEf;
   Double_t        calresp0;

   // List of branches
   TBranch        *b_iev;   //!
   TBranch        *b_neu;   //!
   TBranch        *b_fspl;   //!
   TBranch        *b_tgt;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_A;   //!
   TBranch        *b_hitnuc;   //!
   TBranch        *b_hitqrk;   //!
   TBranch        *b_resid;   //!
   TBranch        *b_sea;   //!
   TBranch        *b_qel;   //!
   TBranch        *b_mec;   //!
   TBranch        *b_res;   //!
   TBranch        *b_dis;   //!
   TBranch        *b_coh;   //!
   TBranch        *b_dfr;   //!
   TBranch        *b_imd;   //!
   TBranch        *b_imdanh;   //!
   TBranch        *b_singlek;   //!
   TBranch        *b_nuel;   //!
   TBranch        *b_em;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_nc;   //!
   TBranch        *b_charm;   //!
   TBranch        *b_neut_code;   //!
   TBranch        *b_nuance_code;   //!
   TBranch        *b_wght;   //!
   TBranch        *b_xs;   //!
   TBranch        *b_ys;   //!
   TBranch        *b_ts;   //!
   TBranch        *b_Q2s;   //!
   TBranch        *b_Ws;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_t;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_EvRF;   //!
   TBranch        *b_Ev;   //!
   TBranch        *b_pxv;   //!
   TBranch        *b_pyv;   //!
   TBranch        *b_pzv;   //!
   TBranch        *b_En;   //!
   TBranch        *b_pxn;   //!
   TBranch        *b_pyn;   //!
   TBranch        *b_pzn;   //!
   TBranch        *b_El;   //!
   TBranch        *b_pxl;   //!
   TBranch        *b_pyl;   //!
   TBranch        *b_pzl;   //!
   TBranch        *b_pl;   //!
   TBranch        *b_cthl;   //!
   TBranch        *b_nfp;   //!
   TBranch        *b_nfn;   //!
   TBranch        *b_nfpip;   //!
   TBranch        *b_nfpim;   //!
   TBranch        *b_nfpi0;   //!
   TBranch        *b_nfkp;   //!
   TBranch        *b_nfkm;   //!
   TBranch        *b_nfk0;   //!
   TBranch        *b_nfem;   //!
   TBranch        *b_nfother;   //!
   TBranch        *b_nip;   //!
   TBranch        *b_nin;   //!
   TBranch        *b_nipip;   //!
   TBranch        *b_nipim;   //!
   TBranch        *b_nipi0;   //!
   TBranch        *b_nikp;   //!
   TBranch        *b_nikm;   //!
   TBranch        *b_nik0;   //!
   TBranch        *b_niem;   //!
   TBranch        *b_niother;   //!
   TBranch        *b_ni;   //!
   TBranch        *b_pdgi;   //!
   TBranch        *b_resc;   //!
   TBranch        *b_Ei;   //!
   TBranch        *b_pxi;   //!
   TBranch        *b_pyi;   //!
   TBranch        *b_pzi;   //!
   TBranch        *b_nf;   //!
   TBranch        *b_pdgf;   //!
   TBranch        *b_Ef;   //!
   TBranch        *b_pxf;   //!
   TBranch        *b_pyf;   //!
   TBranch        *b_pzf;   //!
   TBranch        *b_pf;   //!
   TBranch        *b_cthf;   //!
   TBranch        *b_vtxx;   //!
   TBranch        *b_vtxy;   //!
   TBranch        *b_vtxz;   //!
   TBranch        *b_vtxt;   //!
   TBranch        *b_sumKEf;   //!
   TBranch        *b_calresp0;   //!

   systematics(std::string, std::string, int, int, std::string, int ,TTree *tree=0 );
   virtual ~systematics();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t choice, std::string tweak, int sigma);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void SetFiducialCutParameters(std::string beam_en) {
     fiducialcut->SetFiducialCutParameters(beam_en);
   }
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
   double acceptance_c(double p, double cost, double phi, int id,TFile* file_acceptance);

};

#endif
#ifdef SYSTEMATICS_C

systematics::systematics(std::string a_target,std::string a_beam_en, int number_rotations, int choice, std::string tweak, int sigma,TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {


     fiducialcut = new Fiducial();
     rotation = new Subtraction();
     N_tot = number_rotations;
     ftarget = a_target;
     fbeam_en = a_beam_en;
     fTorusCurrent = 0;
     fchoice = choice;

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("gst",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("gst","systematics");
      if (fchoice == 1) { 
		chain->Add(Form("/pnfs/uboone/persistent/users/apapadop/GenieProduction/R-3_0_6/Rad-R-3_0_6_Clas_%sGeV/G18_10a_02_11a/%s/eresmaid_%s_%s_R-3_0_6-G18_10a_02_11a_RadCorr_*M.root", fbeam_en.c_str(), ftarget.c_str(), ftarget.c_str(), fbeam_en.c_str())); 
		//chain->Add(Form("/w/hallb-scifs17exp/clas/claseg2/apapadop/numaid_%s_%s_hA2018_LFG_FSI_NoRadCorr_3M.root", ftarget.c_str(), fbeam_en.c_str()));
	}
      if (fchoice == 0) { chain->Add(Form("/w/hallb-scifs17exp/clas/claseg2/apapadop/genie_filtered_data_e2a_ep_%s_%s_neutrino6_united4_radphot_test_100M.root",ftarget.c_str(), fbeam_en.c_str())); }

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

systematics::~systematics()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t systematics::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t systematics::LoadTree(Long64_t entry)
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

void systematics::Init(TTree *tree)
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


   fChain->SetBranchAddress("iev", &iev, &b_iev);
   fChain->SetBranchAddress("neu", &neu, &b_neu);
   fChain->SetBranchAddress("fspl", &fspl, &b_fspl);
   fChain->SetBranchAddress("tgt", &tgt, &b_tgt);
   fChain->SetBranchAddress("Z", &Z, &b_Z);
   fChain->SetBranchAddress("A", &A, &b_A);
   fChain->SetBranchAddress("hitnuc", &hitnuc, &b_hitnuc);
   fChain->SetBranchAddress("hitqrk", &hitqrk, &b_hitqrk);
   fChain->SetBranchAddress("resid", &resid, &b_resid);
   fChain->SetBranchAddress("sea", &sea, &b_sea);
   fChain->SetBranchAddress("qel", &qel, &b_qel);
   fChain->SetBranchAddress("mec", &mec, &b_mec);//go down from here
   fChain->SetBranchAddress("res", &res, &b_res);
   fChain->SetBranchAddress("dis", &dis, &b_dis);
   fChain->SetBranchAddress("coh", &coh, &b_coh);
   fChain->SetBranchAddress("dfr", &dfr, &b_dfr);
   fChain->SetBranchAddress("imd", &imd, &b_imd);
   fChain->SetBranchAddress("imdanh", &imdanh, &b_imdanh);
   fChain->SetBranchAddress("singlek", &singlek, &b_singlek);
   fChain->SetBranchAddress("nuel", &nuel, &b_nuel);
   fChain->SetBranchAddress("em", &em, &b_em);
   fChain->SetBranchAddress("cc", &cc, &b_cc);
   fChain->SetBranchAddress("nc", &nc, &b_nc);
   fChain->SetBranchAddress("charm", &charm, &b_charm);
   fChain->SetBranchAddress("neut_code", &neut_code, &b_neut_code);
   fChain->SetBranchAddress("nuance_code", &nuance_code, &b_nuance_code);
   fChain->SetBranchAddress("wght", &wght, &b_wght);
   fChain->SetBranchAddress("xs", &xs, &b_xs);
   fChain->SetBranchAddress("ys", &ys, &b_ys);
   fChain->SetBranchAddress("ts", &ts, &b_ts);
   fChain->SetBranchAddress("Q2s", &Q2s, &b_Q2s);
   fChain->SetBranchAddress("Ws", &Ws, &b_Ws);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("W", &W, &b_W);
   fChain->SetBranchAddress("EvRF", &EvRF, &b_EvRF);
   fChain->SetBranchAddress("Ev", &Ev, &b_Ev);
   fChain->SetBranchAddress("pxv", &pxv, &b_pxv);
   fChain->SetBranchAddress("pyv", &pyv, &b_pyv);
   fChain->SetBranchAddress("pzv", &pzv, &b_pzv);
   fChain->SetBranchAddress("En", &En, &b_En);
   fChain->SetBranchAddress("pxn", &pxn, &b_pxn);
   fChain->SetBranchAddress("pyn", &pyn, &b_pyn);
   fChain->SetBranchAddress("pzn", &pzn, &b_pzn);
   fChain->SetBranchAddress("El", &El, &b_El);
   fChain->SetBranchAddress("pxl", &pxl, &b_pxl);
   fChain->SetBranchAddress("pyl", &pyl, &b_pyl);
   fChain->SetBranchAddress("pzl", &pzl, &b_pzl);
   fChain->SetBranchAddress("pl", &pl, &b_pl);
   fChain->SetBranchAddress("cthl", &cthl, &b_cthl);
   fChain->SetBranchAddress("nfp", &nfp, &b_nfp);
   fChain->SetBranchAddress("nfn", &nfn, &b_nfn);
   fChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
   fChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
   fChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
   fChain->SetBranchAddress("nfkp", &nfkp, &b_nfkp);
   fChain->SetBranchAddress("nfkm", &nfkm, &b_nfkm);
   fChain->SetBranchAddress("nfk0", &nfk0, &b_nfk0);
   fChain->SetBranchAddress("nfem", &nfem, &b_nfem);
   fChain->SetBranchAddress("nfother", &nfother, &b_nfother);
   fChain->SetBranchAddress("nip", &nip, &b_nip);
   fChain->SetBranchAddress("nin", &nin, &b_nin);
   fChain->SetBranchAddress("nipip", &nipip, &b_nipip);
   fChain->SetBranchAddress("nipim", &nipim, &b_nipim);
   fChain->SetBranchAddress("nipi0", &nipi0, &b_nipi0);
   fChain->SetBranchAddress("nikp", &nikp, &b_nikp);
   fChain->SetBranchAddress("nikm", &nikm, &b_nikm);
   fChain->SetBranchAddress("nik0", &nik0, &b_nik0);
   fChain->SetBranchAddress("niem", &niem, &b_niem);
   fChain->SetBranchAddress("niother", &niother, &b_niother);
   fChain->SetBranchAddress("ni", &ni, &b_ni);
   fChain->SetBranchAddress("pdgi", pdgi, &b_pdgi);
   fChain->SetBranchAddress("resc", resc, &b_resc);
   fChain->SetBranchAddress("Ei", Ei, &b_Ei);
   fChain->SetBranchAddress("pxi", pxi, &b_pxi);
   fChain->SetBranchAddress("pyi", pyi, &b_pyi);
   fChain->SetBranchAddress("pzi", pzi, &b_pzi);
   fChain->SetBranchAddress("nf", &nf, &b_nf);
   fChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
   fChain->SetBranchAddress("Ef", Ef, &b_Ef);
   fChain->SetBranchAddress("pxf", pxf, &b_pxf);
   fChain->SetBranchAddress("pyf", pyf, &b_pyf);
   fChain->SetBranchAddress("pzf", pzf, &b_pzf);
   fChain->SetBranchAddress("pf", pf, &b_pf);
   fChain->SetBranchAddress("cthf", cthf, &b_cthf);
   fChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
   fChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
   fChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
   fChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);
   fChain->SetBranchAddress("sumKEf", &sumKEf, &b_sumKEf);
   fChain->SetBranchAddress("calresp0", &calresp0, &b_calresp0);



   Notify();
}

Bool_t systematics::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void systematics::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t systematics::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef systematics_h
