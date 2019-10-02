#ifndef SUBTRACTION_CXX
#define SUBTRACTION_CXX

#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TGraph.h>
#include "Subtraction.h"


void  Subtraction::prot3_rot_func(TVector3  V3prot[3],TVector3  V3prot_uncorr[3],TLorentzVector V4el,double Ecal_3pto2p[][2],double  pmiss_perp_3pto2p[][2],double  P3pto2p[][2],double N_p1[3],double Ecal_3pto1p[3],double  pmiss_perp_3pto1p[3], double *N_p3det){

    double m_prot=0.9382720813;
    const int N_3p=3, N_2p=2;
    double N_p2[N_3p]={0},N_p1det[3][N_3p]={0};
    TVector3 V3_3p_rot[N_3p],V3_2p_rot[N_3p],V3_prot_el_3pto2p[N_2p],V3_prot_el_3pto1p[N_3p];
    double rot_angle;
    double N_pthree=0;
    bool prot_stat[N_3p];
    int count =0;



    for(int g=0; g<N_tot; g++){

         rot_angle=gRandom->Uniform(0,2*TMath::Pi());

         for(int i=0;i<N_3p;i++) {
  	         V3_3p_rot[i]= V3prot_uncorr[i];
  	         V3_3p_rot[i].Rotate(rot_angle,V3q);
         }


         for(int ind_p=0;ind_p<N_3p;ind_p++) prot_stat[ind_p]=PFiducialCut(fbeam_en, V3_3p_rot[ind_p]);

         if( prot_stat[0]  && !prot_stat[1]  && !prot_stat[2])  N_p1[0]=N_p1[0]+1;
         if(!prot_stat[0] &&   prot_stat[1]  && !prot_stat[2])  N_p1[1]=N_p1[1]+1;
         if(!prot_stat[0] &&  !prot_stat[1]  &&  prot_stat[2])  N_p1[2]=N_p1[2]+1;
         if( prot_stat[0]  &&  prot_stat[1]  &&  prot_stat[2])  N_pthree=N_pthree+1;

         if(prot_stat[0] && prot_stat[1] && !prot_stat[2])   N_p2[0]=N_p2[0]+1;
         if(prot_stat[0] && !prot_stat[1] && prot_stat[2])   N_p2[1]=N_p2[1]+1;
         if(!prot_stat[0] && prot_stat[1] && prot_stat[2])   N_p2[2]=N_p2[2]+1;


     }//for loop of 3p rotations ends

   //-----------------------------------------  3p to 1p  -----------------------------------------------------------------------
     for(int j=0;j<N_3p;j++)    { //looping through 1p combinations out of 3protons

         V3_prot_el_3pto1p[j]=V4el.Vect()+ V3prot[j];
         Ecal_3pto1p[j]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[j].Mag()*V3prot[j].Mag())-m_prot+bind_en[target_name];
         pmiss_perp_3pto1p[j]=TMath::Sqrt(V3_prot_el_3pto1p[j].Px()*V3_prot_el_3pto1p[j].Px()+V3_prot_el_3pto1p[j].Py()*V3_prot_el_3pto1p[j].Py());
     }


   //-----------------------------------------  3p to 2p->1p  -----------------------------------------------------------------------
     for(int ind1=0;ind1<N_3p;ind1++){   //looping through 2p combinations  out of 3p
       for(int ind2=0;ind2<N_3p;ind2++){
  	    if(ind1!=ind2 && ind1<ind2){

  	     for(int g1=0; g1<N_tot; g1++){

  	       rot_angle=gRandom->Uniform(0,2*TMath::Pi());

  	       V3_2p_rot[ind1]=V3prot_uncorr[ind1];
  	       V3_2p_rot[ind2]=V3prot_uncorr[ind2];
  	       V3_2p_rot[ind1].Rotate(rot_angle,V3q);
  	       V3_2p_rot[ind2].Rotate(rot_angle,V3q);

  	       if(PFiducialCut(fbeam_en, V3_2p_rot[ind1])  && !PFiducialCut(fbeam_en, V3_2p_rot[ind2])) N_p1det[count][0]=N_p1det[count][0]+1;
  	       if(!PFiducialCut(fbeam_en, V3_2p_rot[ind1]) && PFiducialCut(fbeam_en, V3_2p_rot[ind2]))  N_p1det[count][1]=N_p1det[count][1]+1;
  	       if(PFiducialCut(fbeam_en, V3_2p_rot[ind1])  && PFiducialCut(fbeam_en, V3_2p_rot[ind2]))  N_p1det[count][2]=N_p1det[count][2]+1;
  	     }


  	     if( N_p1det[count][2]!=0   && N_pthree!=0){


  	       V3_prot_el_3pto2p[0]=V4el.Vect()+ V3prot[ind1];
  	       Ecal_3pto2p[count][0]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[ind1].Mag()*V3prot[ind1].Mag())-m_prot+bind_en[target_name];
  	       pmiss_perp_3pto2p[count][0]=TMath::Sqrt(V3_prot_el_3pto2p[0].Px()*V3_prot_el_3pto2p[0].Px()+V3_prot_el_3pto2p[0].Py()*V3_prot_el_3pto2p[0].Py());
  	       P3pto2p[count][0]=N_p1det[count][0]/N_p1det[count][2]*(N_p2[count]/N_pthree);

  	       V3_prot_el_3pto2p[1]=V4el.Vect()+ V3prot[ind2];
  	       Ecal_3pto2p[count][1]=V4el.E()+TMath::Sqrt(m_prot*m_prot+V3prot[ind2].Mag()*V3prot[ind2].Mag())-m_prot+bind_en[target_name];
  	       pmiss_perp_3pto2p[count][1]=TMath::Sqrt(V3_prot_el_3pto2p[1].Px()*V3_prot_el_3pto2p[1].Px()+V3_prot_el_3pto2p[1].Py()*V3_prot_el_3pto2p[1].Py());
  	       P3pto2p[count][1]=N_p1det[count][1]/N_p1det[count][2]*(N_p2[count]/N_pthree);
  	     }


       count=count +1;

  	   }
     }
    }
    *N_p3det=N_pthree;
}


void  Subtraction::prot2_rot_func(TVector3  V3prot[2],TVector3  V3prot_uncorr[2],TLorentzVector V4el,double Ecal_2pto1p[2],double  pmiss_perp_2pto1p[2],double  P2pto1p[2], double *Nboth){


    const int N2=2;
    double rot_angle, N_p2to1[N2]={0},m_prot=0.9382720813;
    TVector3 V3_2prot[N2], V3_prot_el_2pto1p[N2];
    double N_2=0;



    for(int g1=0; g1<N_tot; g1++){


      rot_angle=gRandom->Uniform(0,2*TMath::Pi());

      V3_2prot[0]=V3prot_uncorr[0];
      V3_2prot[1]=V3prot_uncorr[1];
      V3_2prot[0].Rotate(rot_angle,V3q);
      V3_2prot[1].Rotate(rot_angle,V3q);


      if(PFiducialCut(fbeam_en, V3_2prot[0])  && !PFiducialCut(fbeam_en, V3_2prot[1])) N_p2to1[0]=N_p2to1[0]+1;
      if(!PFiducialCut(fbeam_en, V3_2prot[0]) && PFiducialCut(fbeam_en, V3_2prot[1]))  N_p2to1[1]=N_p2to1[1]+1;
      if(PFiducialCut(fbeam_en, V3_2prot[0])  && PFiducialCut(fbeam_en, V3_2prot[1]))  N_2=N_2+1;
    }

     //-----------------------------------------  2p to 1p  -----------------------------------------------------------------------
      V3_prot_el_2pto1p[0]=V4el.Vect()+ V3prot[0];
      Ecal_2pto1p[0]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[0].Mag()*V3prot[0].Mag())-m_prot+bind_en[target_name];
      pmiss_perp_2pto1p[0]=TMath::Sqrt(V3_prot_el_2pto1p[0].Px()*V3_prot_el_2pto1p[0].Px()+V3_prot_el_2pto1p[0].Py()*V3_prot_el_2pto1p[0].Py());

      V3_prot_el_2pto1p[1]=V4el.Vect()+ V3prot[1];
      Ecal_2pto1p[1]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[1].Mag()*V3prot[1].Mag())-m_prot+bind_en[target_name];
      pmiss_perp_2pto1p[1]=TMath::Sqrt(V3_prot_el_2pto1p[1].Px()*V3_prot_el_2pto1p[1].Px()+V3_prot_el_2pto1p[1].Py()*V3_prot_el_2pto1p[1].Py());


    if( N_2!=0){
    P2pto1p[0]=N_p2to1[0]/N_2;
      P2pto1p[1]=N_p2to1[1]/N_2;
    }
    else{
    P2pto1p[0]=0;
      P2pto1p[1]=0;
    }

    *Nboth=N_2;
  }


void Subtraction::prot1_pi1_rot_func(TVector3  V3prot,TVector3 V3pi, int q_pi, double *N_pi_p,double *N_nopi_p){

    double rotation_ang;
    TVector3 V3_pi_rot, V3_p_rot;
    Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;
    bool pi_stat=true;

       double Npi_p = 0;
       double Nnopi_p = 0;

       for(int g=0; g<N_tot; g++){

         rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
         V3_p_rot= V3prot;
         V3_p_rot.Rotate(rotation_ang,V3q);


         V3_pi_rot=V3pi;
         V3_pi_rot.Rotate(rotation_ang,V3q);
         pi_stat=Pi_phot_fid_united(fbeam_en, V3_pi_rot,q_pi);


         if(PFiducialCut(fbeam_en, V3_p_rot)  && pi_stat) Npi_p=Npi_p+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !pi_stat) Nnopi_p=Nnopi_p+1;

       }
       *N_pi_p=Npi_p;
       *N_nopi_p=Nnopi_p;
  }


void Subtraction::prot1_pi2_rot_func(TVector3  V3prot,TVector3 V3pi[2], int q_pi[2], double *P_1p0pi,double P_1p1pi[2]){

    const int N_pi=2;
    double rotation_ang;
    TVector3 V3_rot_pi[2], V3_p_rot;
    bool status_pi[2]={true};

    double N_all = 0;
    double Nnopi = 0,N_1p1pi[2]={0};

       for(int g=0; g<N_tot; g++){

         rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
         V3_p_rot= V3prot;

         V3_p_rot.Rotate(rotation_ang,V3q);

         for(int i=0;i<N_pi;i++){

  	        V3_rot_pi[i]=V3pi[i];
  	        V3_rot_pi[i].Rotate(rotation_ang,V3q);
  	        status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);

         }

         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] ) N_1p1pi[0]=N_1p1pi[0]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] ) N_1p1pi[1]=N_1p1pi[1]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] ) N_all=N_all+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && !status_pi[1]) Nnopi=Nnopi+1;

       }


       if(N_all!=0){
         //----------------------1p2pi->1p0pi
         *P_1p0pi=Nnopi/N_all;

         //----------------------1p2pi->1p1pi->1p0pi
         double N_nopi_p = 0,N_pi_p=0;

         for(int h=0;h<N_pi;h++){

  	 prot1_pi1_rot_func(V3prot,V3pi[h],q_pi[h],&N_pi_p,&N_nopi_p);
  	 if(N_pi_p!=0) P_1p1pi[h]=(N_1p1pi[h]/N_all)*(N_nopi_p/N_pi_p);
  	 else  P_1p1pi[h]=0;
         }
       }   //N_all!=0 statement

       else{
         *P_1p0pi=0;
         P_1p1pi[0]=0;
         P_1p1pi[1]=0;
       }
  }


void Subtraction::prot1_pi3_rot_func(TVector3  V3prot,TVector3 V3pi[3], int q_pi[3], double *P_tot){

    *P_tot=0;
    const int N_pi=3;
    double rotation_ang;
    TVector3 V3_rot_pi[N_pi], V3_p_rot;
    bool status_pi[N_pi]={true};
    double N_all = 0;
    double Nnopi = 0,N_1p1pi[3]={0},N_1p2pi[3]={0};

    for(int g=0; g<N_tot; g++){

         rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
         V3_p_rot= V3prot;

         V3_p_rot.Rotate(rotation_ang,V3q);

         for(int i=0;i<N_pi;i++){

  	        V3_rot_pi[i]=V3pi[i];
  	        V3_rot_pi[i].Rotate(rotation_ang,V3q);
  	        status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
         }

         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] && !status_pi[2]) N_1p1pi[0]=N_1p1pi[0]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] && !status_pi[2]) N_1p1pi[1]=N_1p1pi[1]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && !status_pi[1] && status_pi[2]) N_1p1pi[2]=N_1p1pi[2]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] && !status_pi[2]) N_1p2pi[0]=N_1p2pi[0]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] && status_pi[2]) N_1p2pi[1]=N_1p2pi[1]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] && status_pi[2]) N_1p2pi[2]=N_1p2pi[2]+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] && status_pi[2])  N_all=N_all+1;
         if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && !status_pi[1] && !status_pi[2])Nnopi=Nnopi+1;
    }

       if(N_all!=0){
         //----------------------1p3pi->1p0pi
         double P_1p0pi=0;

         P_1p0pi=Nnopi/N_all;

         //----------------------1p3pi->1p1pi->1p0pi
         double N_nopi_p = 0,N_pi_p=0;
         const int N2pi=2;
         TVector3 V3_pion[N2pi];
         int q_pion[N2pi];

         double  P_1p1pi=0,P_1p2pi=0;

         for(int h=0;h<N_pi;h++){

          prot1_pi1_rot_func(V3prot,V3pi[h],q_pi[h],&N_pi_p,&N_nopi_p);
  	      if(N_pi_p!=0) P_1p1pi=P_1p1pi+(N_1p1pi[h]/N_all)*(N_nopi_p/N_pi_p);

    //----------------------1p3pi->1p2pi->1p0pi
  	     double P_1p1pion[N2pi]={0},P_1p0pion=0;

  	 if(h==0)   {
  	   V3_pion[0]=   V3pi[0];   V3_pion[1]=V3pi[1];
  	   q_pion[0]=    q_pi[0];    q_pion[1]=q_pi[1];

  	 }
  	 if(h==1)   {
  	   V3_pion[0]=   V3pi[0];    V3_pion[1]=V3pi[2];
  	   q_pion[0]=    q_pi[0];     q_pion[1]=q_pi[2];

  	 }
  	 if(h==2)   {
  	   V3_pion[0]=   V3pi[1];   V3_pion[1]=V3pi[2];
  	   q_pion[0]=    q_pi[1];    q_pion[1]=q_pi[2];

  	 }
   //----------------------1p3pi->1p2pi->1p1pi->1p0pi
  	 prot1_pi2_rot_func(V3prot,V3_pion, q_pion,&P_1p0pion,P_1p1pion);
  	 P_1p2pi=P_1p2pi+(N_1p2pi[h]/N_all)*(P_1p0pion-P_1p1pion[0]-P_1p1pion[1]);

       }//for loop ends

       *P_tot=P_1p2pi+P_1p1pi-P_1p0pi;
       }   //N_all!=0 statement

       else{
         *P_tot=0;
       }

  }


void Subtraction::prot2_pi1_rot_func(TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_1pi, int q_pi, TLorentzVector V4_el, double Ecal_2p1pi_to2p0pi[2],double p_miss_perp_2p1pi_to2p0pi[2],double P_2p1pito2p0pi[2],double P_2p1pito1p1pi[2],double P_2p1pito1p0pi[2],double *P_tot){

    const int N_2prot=2;
    TVector3 V3_2p_rotated[N_2prot],V3_1pirot;
    bool pi1_stat=true;
    double N_2p_0pi=0,N_all=0,N_1p_1pi[N_2prot]={0},N_1p_0pi[N_2prot]={0};
    double P_2pto1p[N_2prot]={0},N_2p_det=0;
    double   N_pidet=0,N_piundet=0,rot_angle;
    *P_tot=0;


       for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());


       V3_2p_rotated[0]=V3_2prot_uncorr[0];
       V3_2p_rotated[1]=V3_2prot_uncorr[1];
       V3_2p_rotated[0].Rotate(rot_angle,V3q);
       V3_2p_rotated[1].Rotate(rot_angle,V3q);

       V3_1pirot=V3_1pi;
       V3_1pirot.Rotate(rot_angle,V3q);
       pi1_stat=Pi_phot_fid_united(fbeam_en, V3_1pirot, q_pi);


       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi1_stat) N_2p_0pi=N_2p_0pi+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_1p_1pi[0]=N_1p_1pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_1p_1pi[1]=N_1p_1pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi1_stat) N_1p_0pi[0]=N_1p_0pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi1_stat) N_1p_0pi[1]=N_1p_0pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_all=N_all+1;
     }

   //---------------------------------- 2p 1pi ->2p 0pi   ----------------------------------------------
    if(N_all!=0){

      prot2_rot_func(V3_2prot_corr,V3_2prot_uncorr, V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,P_2pto1p ,&N_2p_det);

      for(int z=0;z<N_2prot;z++){

  	    P_2p1pito2p0pi[z]=(N_2p_0pi/N_all)*P_2pto1p[z];

        //---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

      	prot1_pi1_rot_func(V3_2prot_uncorr[z],V3_1pi, q_pi, &N_pidet,&N_piundet);
        if(N_pidet!=0) P_2p1pito1p1pi[z]=(N_1p_1pi[z]/N_all)*(N_piundet/N_pidet);
        else P_2p1pito1p1pi[z]=0;

        //---------------------------------- 2p 1pi ->1p 0pi   ----------------------------------------------
        P_2p1pito1p0pi[z]=(N_1p_0pi[z]/N_all);

        *P_tot=*P_tot+P_2p1pito2p0pi[z]+P_2p1pito1p1pi[z]-P_2p1pito1p0pi[z];

      }//looping through 2p
      }

  if(N_all==0){
      P_2p1pito2p0pi[0]=P_2p1pito2p0pi[1]=0;
      P_2p1pito1p1pi[0]= P_2p1pito1p1pi[1]=0;
      P_2p1pito1p0pi[0]=P_2p1pito1p0pi[1]=0;
      *P_tot=0;
    }

  }


void Subtraction::prot2_pi2_rot_func(TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_2pi[2], int q_pi[2], TLorentzVector V4_el, double Ecal_2p2pi[2],double p_miss_perp_2p2pi[2],double P_tot_2p[2]){

    const int N_2prot=2,N_2pi=2;
    TVector3 V3_2p_rotated[N_2prot],V3_2pirot[N_2pi];
    bool pi2_stat[N_2pi]={true};
    double   rot_angle;
    double N_2p_0pi=0,N_2p_1pi[N_2pi]={0},N_1p_2pi[N_2prot]={0},N_all=0,N_1p_1pi[N_2prot][N_2pi]={0},N_1p_0pi[N_2prot]={0};
    double   N_pidet=0,N_piundet=0;
    double P_2pto1p[N_2prot]={0},N_2p_det=0;
    double P_1p0pi=0,P_1p1pi[N_2pi]={0};
    double P_2p1pito2p0pi[2]={0}, P_2p1pito1p1pi[2]={0},P_2p1pito1p0pi[2]={0},Ptot=0;
    double P_2p2pito1p0pi[N_2prot]={0},P_2p2pito1p1pi[N_2prot]={0},P_2p2pito1p2pi[N_2prot]={0},P_2p2pito2p1pi[N_2prot]={0};
    P_tot_2p[0]=P_tot_2p[1]=0;

    for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());

       for(int k=0; k<N_2pi; k++){

         V3_2p_rotated[k]=V3_2prot_uncorr[k];
         V3_2p_rotated[k].Rotate(rot_angle,V3q);


  	     V3_2pirot[k]=V3_2pi[k];
  	     V3_2pirot[k].Rotate(rot_angle,V3q);
  	     pi2_stat[k]=Pi_phot_fid_united(fbeam_en, V3_2pirot[k], q_pi[k]);
       }

       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_2p_1pi[0]=N_2p_1pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_2p_1pi[1]=N_2p_1pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && !pi2_stat[1]) N_2p_0pi=N_2p_0pi+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[0]=N_1p_2pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[1]=N_1p_2pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[0][0]=N_1p_1pi[0][0]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[0][1]=N_1p_1pi[0][1]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[1][0]=N_1p_1pi[1][0]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[1][1]=N_1p_1pi[1][1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && !pi2_stat[1])  N_1p_0pi[0]=N_1p_0pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && !pi2_stat[1])  N_1p_0pi[1]=N_1p_0pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_all=N_all+1;
     }


    if(N_all!=0){

      prot2_rot_func(V3_2prot_corr,V3_2prot_uncorr, V4_el,Ecal_2p2pi,p_miss_perp_2p2pi,P_2pto1p ,&N_2p_det);

      for(int z=0;z<N_2prot;z++){
   //---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------

        P_2p2pito1p0pi[z]= N_1p_0pi[z]/N_all;

   //---------------------------------- 2p 2pi ->1p 1pi   ----------------------------------------------

        for(int k=0;k<N_2pi;k++){
          N_pidet=N_piundet=0;
          prot1_pi1_rot_func(V3_2prot_uncorr[z],V3_2pi[k], q_pi[k],&N_pidet,&N_piundet);
          if(N_pidet!=0) P_2p2pito1p1pi[z]=P_2p2pito1p1pi[z]+(N_1p_1pi[z][k]/N_all)*(N_piundet/N_pidet);
        }

   //---------------------------------- 2p 2pi ->2p 0pi   ----------------------------------------------

        if(N_2p_det!=0)P_2p1pito2p0pi[z]=(N_2p_0pi/N_all)*P_2pto1p[z];

        //---------------------------------- 2p 2pi ->1p 2pi   ----------------------------------------------

        P_1p0pi=P_1p1pi[0]=P_1p1pi[1]=0;
        prot1_pi2_rot_func(V3_2prot_uncorr[z],V3_2pi,q_pi,&P_1p0pi,P_1p1pi);
        P_2p2pito1p2pi[z]=(N_1p_2pi[z]/N_all)*(P_1p0pi-P_1p1pi[0]-P_1p1pi[1]);

   //---------------------------------- 2p 2pi ->2p 1pi   ----------------------------------------------

        P_2p1pito2p0pi[0]=P_2p1pito2p0pi[1]=0; P_2p1pito1p1pi[0]=P_2p1pito1p1pi[1]=0; P_2p1pito1p0pi[0]=P_2p1pito1p0pi[1]=0;Ptot=0;
        prot2_pi1_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi[z], q_pi[z],V4_el,Ecal_2p2pi,p_miss_perp_2p2pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

    // P_2p2pito2p1pi[z]=(N_2p_1pi[0]/N_all)*(-P_2p1pito2p0pi[z]- P_2p1pito1p1pi[z]+P_2p1pito1p0pi[z])+(N_2p_1pi[1]/N_all)*(-P_2p1pito2p0pi[z]- P_2p1pito1p1pi[z]+P_2p1pito1p0pi[z]);
   //P_tot_2p[z]=-P_2p2pito1p0pi[z]+P_2p2pito1p1pi[z]+P_2p1pito2p0pi[z]+P_2p2pito1p2pi[z]+P_2p2pito2p1pi[z];

   //P_2p2pito2p1pi[z]=(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[0]- P_2p1pito1p1pi[0]+P_2p1pito1p0pi[0])+(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[1]- P_2p1pito1p1pi[1]+P_2p1pito1p0pi[1]);

        P_2p2pito2p1pi[0]= P_2p2pito2p1pi[0]+(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[0]- P_2p1pito1p1pi[0]+P_2p1pito1p0pi[0]);
        P_2p2pito2p1pi[1]= P_2p2pito2p1pi[1]+(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[1]- P_2p1pito1p1pi[1]+P_2p1pito1p0pi[1]);

      }//looping through 2p

      P_tot_2p[0]=-P_2p2pito1p0pi[0]+P_2p2pito1p1pi[0]+P_2p1pito2p0pi[0]+P_2p2pito1p2pi[0]+P_2p2pito2p1pi[0];
      P_tot_2p[1]=-P_2p2pito1p0pi[1]+P_2p2pito1p1pi[1]+P_2p1pito2p0pi[1]+P_2p2pito1p2pi[1]+P_2p2pito2p1pi[1];

    }//N_all!=0

    if(N_all==0){
      P_tot_2p[0]= P_tot_2p[1]=0;
    }

}


void Subtraction::prot3_pi1_rot_func(TVector3 V3_3prot_corr[3],TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, int q_pi, TLorentzVector V4_el, double Ecal_3p1pi[3],double p_miss_perp_3p1pi[3],double P_tot_3p[3]){

    const int N_3prot=3;
    TVector3 V3_3p_rotated[N_3prot],V3_pirot;
    bool pi_stat=true;
    double   rot_angle;

    double N_1p0pi[N_3prot]={0},N_all=0,N_1p1pi[N_3prot]={0},N_2p0pi[N_3prot]={0},N_2p1pi[N_3prot]={0},N_3p0pi=0;
    double  P_3p1pito1p0pi[N_3prot]={0},P_3p1pito1p1pi[N_3prot]={0};
    double   N_pidet=0,N_piundet=0;
    TVector3 V3_2p_corr[N_3prot],V3_2p_uncorr[N_3prot];
    double Ecal_2p[2],p_miss_perp_2p[2],P_2pto1p[2]={0},N_2p_det=0,P_3p1pito2p0pi[N_3prot]={0};
    int count=0;
    double N_p1[N_3prot]={0},N_p_three=0;
    double E_cal_3pto2p[N_3prot][2]={0},p_miss_perp_3pto2p[N_3prot][2]={0},P_3pto2p[N_3prot][2]={0},P_3p1pito3p0pi[N_3prot]={0};
    double Ecal_2p1pi[2],p_miss_perp_2p1pi[2],P_3p1pito2p1pi[N_3prot]={0};
    double P_2p1pito2p0pi[2]={0},P_2p1pito1p1pi[2]={0},P_2p1pito1p0pi[2]={0},Ptot=0;

       for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());
       for(int k=0; k<N_3prot; k++){

         V3_3p_rotated[k]=V3_3prot_uncorr[k];
         V3_3p_rotated[k].Rotate(rot_angle,V3q);
       }

       V3_pirot=V3_pi;
       V3_pirot.Rotate(rot_angle,V3q);
       pi_stat=Pi_phot_fid_united(fbeam_en, V3_pirot, q_pi);


       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_1p0pi[0]=N_1p0pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_1p0pi[1]=N_1p0pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_1p0pi[2]=N_1p0pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[0]=N_1p1pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[1]=N_1p1pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[2]=N_1p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_2p0pi[0]=N_2p0pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_2p0pi[1]=N_2p0pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_2p0pi[2]=N_2p0pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[0]=N_2p1pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[1]=N_2p1pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[2]=N_2p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_3p0pi=N_3p0pi+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_all=N_all+1;
     }


    if(N_all!=0){

  //----------------------------------3p 1pi ->3p 0pi->1p0pi   ----------------------------------------------
         prot3_rot_func(V3_3prot_uncorr,V3_3prot_corr,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, Ecal_3p1pi,p_miss_perp_3p1pi,&N_p_three);

         if(N_p_three!=0){
  	       P_3p1pito3p0pi[0]= (N_3p0pi/N_all)*(N_p1[0]/N_p_three);
  	       P_3p1pito3p0pi[1]= (N_3p0pi/N_all)*(N_p1[1]/N_p_three);
  	       P_3p1pito3p0pi[2]= (N_3p0pi/N_all)*(N_p1[2]/N_p_three);
  	     }


         for(int z=0;z<N_3prot;z++){
   //---------------------------------- 3p 1pi ->1p 0pi   ----------------------------------------------

         P_3p1pito1p0pi[z]= -N_1p0pi[z]/N_all;

   //---------------------------------- 3p 1pi ->1p 1pi   ----------------------------------------------

         N_pidet=N_piundet=0;
         prot1_pi1_rot_func(V3_3prot_uncorr[z],V3_pi, q_pi, &N_pidet,&N_piundet);
         if(N_pidet!=0) P_3p1pito1p1pi[z]=(N_1p1pi[z]/N_all)*(N_piundet/N_pidet);

   //---------------------------------- 3p 1pi ->2p 0pi   ----------------------------------------------

  		   for(int i=0;i<N_3prot;i++){       //looping through 2p combinations  out of 3p
  		  	if(z!=i && z<i){               // 3 pairs of 2proton combinations with z, i indexes(z<i)

  			  V3_2p_corr[0]=V3_3prot_corr[z];V3_2p_corr[1]=V3_3prot_corr[i];
  			  V3_2p_uncorr[0]=V3_3prot_uncorr[z];V3_2p_uncorr[1]=V3_3prot_uncorr[i];

  			  P_2pto1p[0]=0;P_2pto1p[1]=0;N_2p_det=0;
  			  prot2_rot_func(V3_2p_corr,V3_2p_uncorr, V4_el,Ecal_2p,p_miss_perp_2p,P_2pto1p ,&N_2p_det);
  			  if(N_2p_det!=0){
  			    P_3p1pito2p0pi[z]=P_3p1pito2p0pi[z]+(N_2p0pi[count]/N_all)*P_2pto1p[0];
  			    P_3p1pito2p0pi[i]=P_3p1pito2p0pi[i]+(N_2p0pi[count]/N_all)*P_2pto1p[1];
  			  }


        //---------------------------------- 3p 1pi ->3p 0pi->2p0pi   ----------------------------------------------

  	      if(N_p_three!=0){
  	       P_3p1pito3p0pi[z]= P_3p1pito3p0pi[z]+(N_3p0pi/N_all)*(-P_3pto2p[count][0]);
  	       P_3p1pito3p0pi[i]= P_3p1pito3p0pi[i]+(N_3p0pi/N_all)*(-P_3pto2p[count][1]);
  	      }

   //---------------------------------- 3p 1pi ->2p 1pi   ----------------------------------------------
          P_2p1pito2p0pi[0]=P_2p1pito2p0pi[1]=0; P_2p1pito1p1pi[0]=P_2p1pito1p1pi[1]=0; P_2p1pito1p0pi[0]=P_2p1pito1p0pi[1]=0;Ptot=0;

          prot2_pi1_rot_func(V3_2p_corr,V3_2p_uncorr,V3_pi, q_pi,V4_el,Ecal_2p1pi,p_miss_perp_2p1pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

   // P_3p1pito2p1pi[z]=(N_2p1pi[count]/N_all)*(-P_2p1pito2p0pi[z]- P_2p1pito1p1pi[z]+P_2p1pito1p0pi[z]);
   // P_3p1pito2p1pi[i]=(N_2p1pi[count]/N_all)*(-P_2p1pito2p0pi[i]- P_2p1pito1p1pi[i]+P_2p1pito1p0pi[i]);

          P_3p1pito2p1pi[z]= P_3p1pito2p1pi[z]+(N_2p1pi[count]/N_all)*(-P_2p1pito2p0pi[0]- P_2p1pito1p1pi[0]+P_2p1pito1p0pi[0]);
          P_3p1pito2p1pi[i]= P_3p1pito2p1pi[i]+(N_2p1pi[count]/N_all)*(-P_2p1pito2p0pi[1]- P_2p1pito1p1pi[1]+P_2p1pito1p0pi[1]);

  	      count=count+1;
  		  	}
  		  }

  		      // P_tot_3p[z]=P_3p1pito2p1pi[z]+P_3p1pito3p0pi[z]+P_3p1pito2p0pi[z]+P_3p1pito1p1pi[z]+ P_3p1pito1p0pi[z];

      }//looping through 3p

      P_tot_3p[0]=P_3p1pito2p1pi[0]+P_3p1pito3p0pi[0]+P_3p1pito2p0pi[0]+P_3p1pito1p1pi[0]+ P_3p1pito1p0pi[0];
      P_tot_3p[1]=P_3p1pito2p1pi[1]+P_3p1pito3p0pi[1]+P_3p1pito2p0pi[1]+P_3p1pito1p1pi[1]+ P_3p1pito1p0pi[1];
      P_tot_3p[2]=P_3p1pito2p1pi[2]+P_3p1pito3p0pi[2]+P_3p1pito2p0pi[2]+P_3p1pito1p1pi[2]+ P_3p1pito1p0pi[2];

      }

      if(N_all==0){
        P_tot_3p[0]= P_tot_3p[1]=P_tot_3p[2]=0;
      }

  }


void Subtraction::pi1_rot_func(TVector3 V3_pi, int q_pi, double *P_pi){


    double N_pion=0;
    double rot_angle;
    TVector3 V3_rot_pi;
    Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;


    for(int g=0; g<N_tot; g++){

      V3_rot_pi=V3_pi;
      rot_angle=gRandom->Uniform(0,2*TMath::Pi());
      V3_rot_pi.Rotate(rot_angle,V3q);
      if(Pi_phot_fid_united(fbeam_en, V3_rot_pi,q_pi)) N_pion=N_pion+1;
    }

    if(N_pion!=0)     *P_pi=(N_tot-N_pion)/N_pion;
    else *P_pi=0;

    // if(!radstat) cout<<"nereqev     "<<N_pion<<"   radstat"<<"     "<<radstat<<endl;
  }


void Subtraction::pi2_rot_func(TVector3 V3_pi[2], int q_pi[2], double *P_0pi,double P_1pi[2]){

    const int N_pi=2;
    TVector3 V3_rot_pi[N_pi];
    double rot_angle;
    bool status_pi[N_pi]={true};
    Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;
    double N_bothpi=0,N_nopi=0,N_1pi[N_pi]={0},P_pi1[N_pi]={0};


    for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());
       for(int i=0;i<N_pi;i++){

  	      V3_rot_pi[i]=V3_pi[i];
  	      V3_rot_pi[i].Rotate(rot_angle,V3q);
  	      status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);

       }
       if( status_pi[0] && !status_pi[1]) N_1pi[0]=N_1pi[0]+1;
       if(!status_pi[0] &&  status_pi[1]) N_1pi[1]=N_1pi[1]+1;
       if(!status_pi[0] && !status_pi[1]) N_nopi=N_nopi+1;
       if( status_pi[0] &&  status_pi[1]) N_bothpi=N_bothpi+1;
    }


    pi1_rot_func(V3_pi[0],  q_pi[0], P_pi1);
    pi1_rot_func(V3_pi[1],  q_pi[1], P_pi1+1);

    if(N_bothpi!=0){
  	     *P_0pi=N_nopi/N_bothpi;
  	     P_1pi[0]=N_1pi[0]/N_bothpi*P_pi1[0];
  	     P_1pi[1]=N_1pi[1]/N_bothpi*P_pi1[1];
    }
    else{
  	     *P_0pi=0;
  	     P_1pi[0]=0;
  	     P_1pi[1]=0;
    }
}

void Subtraction::pi3_rot_func(TVector3 V3_pi[3], int q_pi[3], double *P_0pi, double P_1pi[3],double P_320[3],double P_3210[][2]){

   const int N_pi=3;
   TVector3 V3_rot_pi[N_pi];
   double rot_angle;
   bool status_pi[N_pi]={true};
   Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;
   double N_1pi[N_pi]={0},N_allpi=0,N_nopi=0,N_2pi[N_pi]={0};



   for(int g=0; g<N_tot; g++){

      rot_angle=gRandom->Uniform(0,2*TMath::Pi());

      for(int i=0;i<N_pi;i++){


  	    V3_rot_pi[i]=V3_pi[i];
  	    V3_rot_pi[i].Rotate(rot_angle,V3q);
  	    status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);

      }

         if( status_pi[0]  && !status_pi[1] &&  !status_pi[2]) N_1pi[0]=N_1pi[0]+1;
         if(!status_pi[0] &&   status_pi[1]  && !status_pi[2]) N_1pi[1]=N_1pi[1]+1;
         if(!status_pi[0] &&  !status_pi[1] &&   status_pi[2]) N_1pi[2]=N_1pi[2]+1;
         if( status_pi[0]  &&  status_pi[1] &&  !status_pi[2]) N_2pi[0]=N_2pi[0]+1;
         if( status_pi[0] &&  !status_pi[1]  &&  status_pi[2]) N_2pi[1]=N_2pi[1]+1;
         if(!status_pi[0] &&   status_pi[1] &&   status_pi[2]) N_2pi[2]=N_2pi[2]+1;
         if(!status_pi[0] &&  !status_pi[1] &&  !status_pi[2]) N_nopi=N_nopi+1;
         if( status_pi[0]  &&  status_pi[1]  &&  status_pi[2]) N_allpi=N_allpi+1;
    }

    const int N_pi2=2;
    double P_pi=0;
    double P_1pion[N_pi2]={0},P_0pion=0;
    TVector3 V3_pion[N_pi2];
    int q_pion[N_pi2];

    if(N_allpi!=0){
   //---------------------------3pi->0pi----------------------------------------------
      *P_0pi=N_nopi/N_allpi;
   //---------------------------3pi->1pi->0pi----------------------------------------------
      for(int h=0;h<N_pi;h++){
        pi1_rot_func( V3_pi[h],q_pi[h],&P_pi);
        P_1pi[h]=P_pi*(N_1pi[h]/N_allpi);
      //---------------------------3pi->2pi->0pi----------------------------------------------

        if(h==0)   {
  	       V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[1];
  	       q_pion[0]=q_pi[0];q_pion[1]=q_pi[1];

        }
        if(h==1)   {
  	       V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[2];
  	       q_pion[0]=q_pi[0];q_pion[1]=q_pi[2];
        }
        if(h==2)   {
  	       V3_pion[0]=V3_pi[1];V3_pion[1]=V3_pi[2];
  	       q_pion[0]=q_pi[1];q_pion[1]=q_pi[2];
        }
        pi2_rot_func( V3_pion,q_pion,&P_0pion, P_1pion);
        P_320[h]=P_0pion*(N_2pi[h]/N_allpi);

  //---------------------------3pi->2pi->1pi->0pi----------------------------------------------

      P_3210[h][0]=P_1pion[0]*(N_2pi[h]/N_allpi);
      P_3210[h][1]=P_1pion[1]*(N_2pi[h]/N_allpi);

      }//end of 3p loop


    }// end of N_allpi!=0 statement


    else{
      for(int h=0;h<3;h++){
        P_320[h]=0;
        P_3210[h][0]=0;
        P_3210[h][1]=0;
      }
    }


  }


void Subtraction::pi4_rot_func(TVector3 V3_pi[4], int q_pi[4], double *P_0pi,double *P_410,double *P_420,double *P_4210,double *P_430,double *P_4310,double *P_4320,double *P_43210){

   const int N_pi=4;
    TVector3 V3_rot_pi[N_pi];
   double rot_angle;
   bool status_pi[N_pi]={true};
   double N_1pi[N_pi]={0},N_allpi=0,N_nopi=0,N_2pi[6]={0},N_3pi[4]={0};



    for(int g=0; g<N_tot; g++){

         rot_angle=gRandom->Uniform(0,2*TMath::Pi());

         for(int i=0;i<N_pi;i++){

  	        V3_rot_pi[i]=V3_pi[i];
  	        V3_rot_pi[i].Rotate(rot_angle,V3q);
  	        status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);

         }

         if( status_pi[0]  && !status_pi[1] &&  !status_pi[2]  &&  !status_pi[3]) N_1pi[0]=N_1pi[0]+1; //1pi or phot
         if( !status_pi[0]  && status_pi[1] &&  !status_pi[2]  &&  !status_pi[3]) N_1pi[1]=N_1pi[1]+1;
         if( !status_pi[0]  && !status_pi[1] &&  status_pi[2]  &&  !status_pi[3]) N_1pi[2]=N_1pi[2]+1;
         if( !status_pi[0]  && !status_pi[1] &&  !status_pi[2]  &&  status_pi[3]) N_1pi[3]=N_1pi[3]+1;

         if( status_pi[0]  && status_pi[1] &&  !status_pi[2]  &&  !status_pi[3]) N_2pi[0]=N_2pi[0]+1;//2pi or phot
         if( status_pi[0]  &&! status_pi[1] &&  status_pi[2]  &&  !status_pi[3]) N_2pi[1]=N_2pi[1]+1;
         if( status_pi[0]  &&! status_pi[1] &&  !status_pi[2]  &&  status_pi[3]) N_2pi[2]=N_2pi[2]+1;
         if( !status_pi[0]  && status_pi[1] &&  status_pi[2]  &&  !status_pi[3]) N_2pi[3]=N_2pi[3]+1;
         if( !status_pi[0]  && status_pi[1] &&  !status_pi[2]  &&  status_pi[3]) N_2pi[4]=N_2pi[4]+1;
         if( !status_pi[0]  && !status_pi[1] &&  status_pi[2]  &&  status_pi[3]) N_2pi[5]=N_2pi[5]+1;

         if( status_pi[0]  && status_pi[1] &&  status_pi[2]  &&  !status_pi[3]) N_3pi[0]=N_3pi[0]+1;//3pi or phot
         if( status_pi[0]  && status_pi[1] &&  !status_pi[2]  &&  status_pi[3]) N_3pi[1]=N_3pi[1]+1;
         if( status_pi[0]  && !status_pi[1] &&  status_pi[2]  &&  status_pi[3]) N_3pi[2]=N_3pi[2]+1;
         if( !status_pi[0]  && status_pi[1] &&  status_pi[2]  &&  status_pi[3]) N_3pi[3]=N_3pi[3]+1;

         if( !status_pi[0]  && !status_pi[1] &&  !status_pi[2]  &&  !status_pi[3]) N_nopi=N_nopi+1; //0 pi or phot
         if( status_pi[0]  && status_pi[1] &&  status_pi[2]  && status_pi[3]) N_allpi=N_allpi+1; //4pi or phot
    }


    double P_pi=0;
    const int N_pi3=3;
    TVector3 V3_pion[N_pi3];
    double P_1pion[N_pi3]={0},P_0pion=0, P_320_pion[3]={0}, P_3210_pion[3][2]={0};
    int q_pion[N_pi3];

    if(N_allpi!=0){
   //---------------------------4pi->0pi----------------------------------------------
      *P_0pi=N_nopi/N_allpi;
   //---------------------------4pi->1pi->0pi----------------------------------------------
      for(int h=0;h<N_pi;h++){
        pi1_rot_func(V3_pi[h],q_pi[h],&P_pi);
        *P_410=*P_410+P_pi*(N_1pi[h]/N_allpi);

      //---------------------------4pi->3pi->0pi----------------------------------------------
        if(h==0)   {
  	       V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[1];V3_pion[2]=V3_pi[2];
  	       q_pion[0]=q_pi[0];q_pion[1]=q_pi[1];q_pion[2]=q_pi[2];
        }
        if(h==1)   {
  	       V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[1];V3_pion[2]=V3_pi[3];
  	       q_pion[0]=q_pi[0];q_pion[1]=q_pi[1];q_pion[2]=q_pi[3];
        }
        if(h==2)   {
  	       V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[2];V3_pion[2]=V3_pi[3];
  	       q_pion[0]=q_pi[0];q_pion[1]=q_pi[2];q_pion[2]=q_pi[3];
        }
        if(h==3)   {
  	       V3_pion[0]=V3_pi[1];V3_pion[1]=V3_pi[2];V3_pion[2]=V3_pi[3];
  	       q_pion[0]=q_pi[1];q_pion[1]=q_pi[2];q_pion[2]=q_pi[3];
        }

        pi3_rot_func(V3_pion,q_pion,&P_0pion, P_1pion,P_320_pion,P_3210_pion);
        *P_430=*P_430+P_0pion*(N_3pi[h]/N_allpi);

  //---------------------------4pi->3pi->1pi->0pi----------------------------------------------

        *P_4310= *P_4310+(P_1pion[0]+P_1pion[1]+P_1pion[2])*(N_3pi[h]/N_allpi);

  //---------------------------4pi->3pi->2pi->0pi----------------------------------------------

        *P_4320=*P_4320+(P_320_pion[0]+P_320_pion[1]+P_320_pion[2])*(N_3pi[h]/N_allpi);

  //---------------------------4pi->3pi->2pi->1pi->0pi----------------------------------------------

        *P_43210=*P_43210+((P_3210_pion[0][0]+P_3210_pion[0][1])+(P_3210_pion[1][0]+P_3210_pion[1][1])+
  			 (P_3210_pion[2][0]+P_3210_pion[2][1]))*(N_3pi[h]/N_allpi);

      }//end of 4pi loop


  //---------------------------4pi->2pi->0pi----------------------------------------------
    const int N2pi=2;
    TVector3 V3_2pi[N2pi];
    int q_2pi[N2pi];
    double P_0pi=0, P_1pi[N2pi]={0};

    for(int h=0;h<6;h++){

      if(h==0)   {
         V3_2pi[0]=V3_pi[0];      V3_2pi[1]=V3_pi[1];
         q_2pi[0]=q_pi[0];        q_2pi[1]=q_pi[1];

      }
      if(h==1)   {
         V3_2pi[0]=V3_pi[0];      V3_2pi[1]=V3_pi[2];
         q_2pi[0]=q_pi[0];        q_2pi[1]=q_pi[2];
      }
      if(h==2)   {
         V3_2pi[0]=V3_pi[0];      V3_2pi[1]=V3_pi[3];
         q_2pi[0]=q_pi[0];        q_2pi[1]=q_pi[3];
      }
      if(h==3)   {
         V3_2pi[0]=V3_pi[1];      V3_2pi[1]=V3_pi[2];
         q_2pi[0]=q_pi[1];        q_2pi[1]=q_pi[2];
      }
      if(h==4)   {
         V3_2pi[0]=V3_pi[1];      V3_2pi[1]=V3_pi[3];
         q_2pi[0]=q_pi[1];        q_2pi[1]=q_pi[3];
      }
      if(h==5)   {
         V3_2pi[0]=V3_pi[2];      V3_2pi[1]=V3_pi[3];
         q_2pi[0]=q_pi[2];        q_2pi[1]=q_pi[3];
      }

      pi2_rot_func( V3_2pi,q_2pi,&P_0pi, P_1pi);

      *P_420=*P_420+P_0pi*(N_2pi[h]/N_allpi);

  //---------------------------4pi->2pi->1pi->0pi----------------------------------------------

      *P_4210=*P_4210+(P_1pi[0]+P_1pi[0])*(N_2pi[h]/N_allpi);

    }

   }// end of N_allpi!=0 statement

   else{
     *P_0pi=0;
     *P_410=0;
     *P_430=0;
     *P_4310=0;
     *P_4320=0;
     *P_43210=0;
     *P_420=0;
     *P_4210=0;
  }



}

#endif
