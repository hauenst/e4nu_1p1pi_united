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

void Subtraction::prot1_pi2_rot_func(TVector3 V3_el, TVector3  V3prot, TVector3 V3pi[2], int q_pi[2], double P_1p1pi[2]){

    const int N_pi=2;
    double rotation_ang;
    TVector3 V3_rot_pi[2], V3_p_rot;
    bool status_pi[2]={true};

    double N_all = 0;
    double N_1p1pi[2]={0};

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

       }


       if(N_all!=0){
         //----------------------1p2pi->1p1pi
         for(int h=0;h<N_pi;h++){
  	        P_1p1pi[h]=(N_1p1pi[h]/N_all);
         }
       }   //N_all!=0 statement

       else{
         P_1p1pi[0]=0;
         P_1p1pi[1]=0;
       }
  }

//This function below still incomplete. Need to call other functions and calculate weights.
void Subtraction::prot1_pi3_rot_func(TVector3  V3prot, TVector3 V3pi[3], int q_pi[3], double *P_tot){
    const int N_pi=3;
    double P_1p3pito1p1pi[3] = {0};
    double rotation_ang;
    TVector3 V3_rot_pi[N_pi], V3_p_rot;
    bool status_pi[N_pi]={true};
    double N_all = 0;
    double N_1p1pi[3]={0},N_1p2pi[3]={0};

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
    }
    for(int i=0;i<3;i++)
    {
      if(N_all!=0)
      {
        P_1p3pito1p1pi[i] = (N1pi1p[i]/N_all);
      }
      else
      {
        P_1p3pito1p1pi[i] = 0;
      }
      P_tot += P_1p3pito1p1pi[i];
    }
    //--------------1p2pi->1p1pi------
    double P_1p2pito1p1pi[6] = {0};
    int count = 0;
    double V3pi2[2] = {0};
    int q_pi2[2] = {0};
    double P_1p1pi[2] = {0};
    for(int i=0; i<3; i++)
        {
          for(int j=0; j<3; j++)
          {
            if(i<j)
            {
              P_1p1pi = {0};
              V3pi2[0] = V3pi[i];
              V3pi2[1] = V3pi[j];
              q_pi2[0] = q_pi[i];
              q_pi2[1] = q_pi[j];
              prot1_pi2_rot_func(V3prot, V3pi, q_pi, P_1p1pi);
              P_1p3pito1p1pi[count] = P_1p1pi[0];
              P_1p3pito1p1pi[count+1] = P_1p1pi[1];
              P_tot -= (P_1p3pito1p1pi[count] + P_1p3pito1p1pi[count+1]);
              count = count+2;
            }
          }
        }
}


void Subtraction::prot2_pi1_rot_func(TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_1pi, int q_pi, TLorentzVector V4_el, double Ecal_2p1pi_to2p0pi[2],double p_miss_perp_2p1pi_to2p0pi[2],double P_2p1pito2p0pi[2],double P_2p1pito1p1pi[2],double P_2p1pito1p0pi[2],double *P_tot){

    const int N_2prot=2;
    TVector3 V3_2p_rotated[N_2prot],V3_1pirot;
    bool pi1_stat=true;
    double N_all=0,N_1p_1pi[N_2prot]={0};
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



       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_1p_1pi[0]=N_1p_1pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_1p_1pi[1]=N_1p_1pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_all=N_all+1;
     }


  }


void Subtraction::prot2_pi2_rot_func(TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_2pi[2], int q_pi[2], TLorentzVector V4_el, double Ecal_2p2pi[2],double p_miss_perp_2p2pi[2],double P_tot_2p[2]){

    const int N_2prot=2,N_2pi=2;
    TVector3 V3_2p_rotated[N_2prot],V3_2pirot[N_2pi];
    bool pi2_stat[N_2pi]={true};
    double   rot_angle;
    double N_2p_1pi[N_2pi]={0},N_1p_2pi[N_2prot]={0},N_all=0,N_1p_1pi[N_2prot][N_2pi]={0};
    double   N_pidet=0,N_piundet=0;
    double P_2pto1p[N_2prot]={0},N_2p_det=0;
    double P_1p1pi[N_2pi]={0};
    double P_2p1pito1p1pi[2]={0},Ptot=0;
    double P_2p2pito1p1pi[N_2prot]={0},P_2p2pito1p2pi[N_2prot]={0},P_2p2pito2p1pi[N_2prot]={0};
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
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[0]=N_1p_2pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[1]=N_1p_2pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[0][0]=N_1p_1pi[0][0]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[0][1]=N_1p_1pi[0][1]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[1][0]=N_1p_1pi[1][0]+1;
       if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[1][1]=N_1p_1pi[1][1]+1;
       if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_all=N_all+1;
     }


    if(N_all!=0){

      prot2_rot_func(V3_2prot_corr,V3_2prot_uncorr, V4_el,Ecal_2p2pi,p_miss_perp_2p2pi,P_2pto1p ,&N_2p_det);

      for(int z=0;z<N_2prot;z++){
   //---------------------------------- 2p 2pi ->1p 1pi   ----------------------------------------------

        for(int k=0;k<N_2pi;k++){
          N_pidet=N_piundet=0;
          prot1_pi1_rot_func(V3_2prot_uncorr[z],V3_2pi[k], q_pi[k],&N_pidet,&N_piundet);
          if(N_pidet!=0) P_2p2pito1p1pi[z]=P_2p2pito1p1pi[z]+(N_1p_1pi[z][k]/N_all)*(N_piundet/N_pidet);
        }

   //---------------------------------- 2p 2pi ->2p 1pi   ----------------------------------------------

       P_2p1pito1p1pi[0]=P_2p1pito1p1pi[1]=0;Ptot=0;
        //prot2_pi1_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi[z], q_pi[z],V4_el,Ecal_2p2pi,p_miss_perp_2p2pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot); ---DID NOT WANT TO REMOVE UNSURE OF HOW TO PROCEED

    // P_2p2pito2p1pi[z]=(N_2p_1pi[0]/N_all)*(-P_2p1pito2p0pi[z]- P_2p1pito1p1pi[z]+P_2p1pito1p0pi[z])+(N_2p_1pi[1]/N_all)*(-P_2p1pito2p0pi[z]- P_2p1pito1p1pi[z]+P_2p1pito1p0pi[z]);
   //P_tot_2p[z]=-P_2p2pito1p0pi[z]+P_2p2pito1p1pi[z]+P_2p1pito2p0pi[z]+P_2p2pito1p2pi[z]+P_2p2pito2p1pi[z];

   //P_2p2pito2p1pi[z]=(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[0]- P_2p1pito1p1pi[0]+P_2p1pito1p0pi[0])+(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[1]- P_2p1pito1p1pi[1]+P_2p1pito1p0pi[1]);


      }//looping through 2p



    }//N_all!=0

    if(N_all==0){
      P_tot_2p[0]= P_tot_2p[1]=0;
    }

}


void Subtraction::prot3_pi1_rot_func(TVector3 V3_3prot_corr[3],TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, int q_pi, TLorentzVector V4_el,double P_tot_3p[3]){

    const int N_3prot=3;
    TVector3 V3_3p_rotated[N_3prot],V3_pirot;
    bool pi_stat=true;
    double   rot_angle;

    double N_all=0,N_1p1pi[N_3prot]={0},N_2p1pi[N_3prot]={0};
    double  P_3p1pito1p1pi[N_3prot]={0};
    double   N_pidet=0,N_piundet=0;
    TVector3 V3_2p_corr[N_3prot],V3_2p_uncorr[N_3prot];
    double P_2pto1p[2]={0},N_2p_det=0;
    int count=0;
    double N_p1[N_3prot]={0},N_p_three=0;
    double P_3pto2p[N_3prot][2]={0};
    double P_3p1pito2p1pi[N_3prot]={0};
    double P_2p1pito1p1pi[2]={0},Ptot=0;

       for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());
       for(int k=0; k<N_3prot; k++){

         V3_3p_rotated[k]=V3_3prot_uncorr[k];
         V3_3p_rotated[k].Rotate(rot_angle,V3q);
       }

       V3_pirot=V3_pi;
       V3_pirot.Rotate(rot_angle,V3q);
       pi_stat=Pi_phot_fid_united(fbeam_en, V3_pirot, q_pi);


       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[0]=N_1p1pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[1]=N_1p1pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[2]=N_1p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[0]=N_2p1pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[1]=N_2p1pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[2]=N_2p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_all=N_all+1;
     }


    if(N_all!=0)
    {
      for(int z=0;z<N_3prot;z++){

   //---------------------------------- 3p 1pi ->1p 1pi   ----------------------------------------------

         N_pidet=N_piundet=0;
         prot1_pi1_rot_func(V3_3prot_uncorr[z],V3_pi, q_pi, &N_pidet,&N_piundet);
         if(N_pidet!=0) P_3p1pito1p1pi[z]=(N_1p1pi[z]/N_all)*(N_piundet/N_pidet);

   //---------------------------------- 3p 1pi ->2p 0pi   ----------------------------------------------

  		   for(int i=0;i<N_3prot;i++){       //looping through 2p combinations  out of 3p
  		  	if(z!=i && z<i){               // 3 pairs of 2proton combinations with z, i indexes(z<i)

   //---------------------------------- 3p 1pi ->2p 1pi   ----------------------------------------------
          P_2p1pito1p1pi[0]=P_2p1pito1p1pi[1]=0;
          Ptot=0;

          prot2_pi1_rot_func(V3_2p_corr,V3_2p_uncorr,V3_pi, q_pi, V4_el, P_2p1pito1p1pi,&Ptot);

          P_3p1pito2p1pi[z]= P_3p1pito2p1pi[z]+(N_2p1pi[count]/N_all)*(-P_2p1pito1p1pi[0]);
          P_3p1pito2p1pi[i]= P_3p1pito2p1pi[i]+(N_2p1pi[count]/N_all)*(-P_2p1pito1p1pi[1]);

  	      count=count+1;
  		  	}
  		  }
      }//looping through 3p

      P_tot_3p[0]=P_3p1pito2p1pi[0]+P_3p1pito1p1pi[0];
      P_tot_3p[1]=P_3p1pito2p1pi[1]+P_3p1pito1p1pi[1];
      P_tot_3p[2]=P_3p1pito2p1pi[2]+P_3p1pito1p1pi[2];

      }

      if(N_all==0){
        P_tot_3p[0]= P_tot_3p[1]=P_tot_3p[2]=0;
      }

  }

#endif
