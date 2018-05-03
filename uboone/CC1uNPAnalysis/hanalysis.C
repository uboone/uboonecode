#define hanalysis_cxx
// The class definition in hanalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("hanalysis.C")
// root> T->Process("hanalysis.C","some options")
// root> T->Process("hanalysis.C+")
//


#include "hanalysis.h"
#include <TH2.h>
#include "TH2D.h"
#include "TH1D.h"
#include "TVector3.h"
#include <TStyle.h>

//histograms declarations here

TH1D *fvertex_x;
TH1D *fvertex_y;
TH1D *fvertex_z;

TH1D *fmucand_startx;
TH1D *fmucand_starty;
TH1D *fmucand_startz;
TH1D *fmucand_endx;
TH1D *fmucand_endy;
TH1D *fmucand_endz;

TH1D *fpcand_startx;
TH1D *fpcand_starty;
TH1D *fpcand_startz;
TH1D *fpcand_endx;
TH1D *fpcand_endy;
TH1D *fpcand_endz;

TH1D *fEvis;
TH1D *fQ2cal;

//%%%%%%%%%%%%%%%%%%
TH1D *h_thetapp;
TH1D *h_thetapps[9];
TH1D *h_thetapps_sig[9];
TH1D *h_phipp;
TH1D *h_phipps[9];
TH1D *h_phipps_sig[9];
//==================================================
TH2D *muon_lenvsdqdx;
TH2D *proton_lenvsdqdx;

TH1D *trunmean_cand; 
TH1D *trunmean_pcand;

TH1D *trunmean_muon[9];
TH1D *trunmean_proton[9];
TH1D *trunmean_muon_sig[9];
TH1D *trunmean_proton_sig[9];

////======================================
TH1D *trklen_cand;
TH1D *trklen_pcand;
TH1D *trklen_muon[9];
TH1D *trklen_proton[9];

TH1D *trklen_muon_sig[9];
TH1D *trklen_proton_sig[9];


//====================================
TH1D *thetalep;
TH1D *thetahad;
//=========================================
TH1D *philep;
TH1D *phihad;
TH1D *phi_muon[9];
TH1D *phi_proton[9];
TH1D *phi_muon_sig[9];
TH1D *phi_proton_sig[9];
//==================================

TH1D *thetalep_muon;
TH1D *thetahad_proton;
TH1D *philep_muon;
TH1D *phihad_proton;

//=====================================
TH1D *costhetalep;
TH1D *costhetahad;

TH1D *costheta_muon[9];
TH1D *costheta_proton[9];

TH1D *costheta_muon_sig[9];
TH1D *costheta_proton_sig[9];


TH1D *costhetatruelep;
TH1D *costhetatruehad;
TH1D *costhetaresolep;
TH1D *costhetaresohad;


TH2D *costheta_vs_trklen;

TH1D *ntrkspreco;
TH1D *ntrks_proton[9];
TH1D *ntrks_proton_sig[9];


TH1D *costhetatrue_muon[9];
TH1D *costhetatrue_proton[9];
TH1D *costhetatrue_muon_sig[9];
TH1D *costhetatrue_proton_sig[9];


TH1D *costhetareso_muon[9];
TH1D *costhetareso_proton[9];
TH1D *costhetareso_muon_sig[9];
TH1D *costhetareso_proton_sig[9];


TH1D *cosphiresolep;
TH1D *cosphiresohad;

TH1D *cosphireso_muon[9];
TH1D *cosphireso_proton[9];
TH1D *cosphireso_muon_sig[9];
TH1D *cosphireso_proton_sig[9];
///==================================================

//=======================================================
TH1D *trkmom_muon[9];
TH1D *trkmom_muon_sig[9];
TH1D *trkmomlep;

//========================================================= 
TH1D *Plep;
TH1D *Phad;


TH1D *Plep_muon;
TH1D *Phad_proton;

TH1D *Presolep;
TH1D *Presohad;

TH1D *fPmuon[9];
TH1D *fPproton[9];

TH1D *fPmuon_sig[9];
TH1D *fPproton_sig[9];


TH1D *fPresomuon[9];
TH1D *fPresoproton[9];

TH1D *fPresomuon_sig[9];
TH1D *fPresoproton_sig[9];
///===================================================

TH1D *h_thetamup;
TH1D *h_thetamups[9];
TH1D *h_thetamups_sig[9];
TH1D *h_phimup;
TH1D *h_phimups[9];
TH1D *h_phimups_sig[9];
//==========================================
TH1D *Nhitslep;
TH1D *Nhitshad;
TH1D *Nhitsmuon[9];
TH1D *Nhitsproton[9];

TH1D *Nhitsmuon_sig[9];
TH1D *Nhitsproton_sig[9];
//============================================

TH2D *Pmutruevsreco;
TH2D *Pptruevsreco;

TH2D *Costhetamutruevsreco;
TH2D *Cosphimutruevsreco;
TH2D *Costhetaptruevsreco;
TH2D *Cosphiptruevsreco;


TH2D *ntrksp_truevsreco;
//===================================================
TH1D *fvertexx[9];
TH1D *fvertexy[9];
TH1D *fvertexz[9];

TH1D *fmucandstartx[9];
TH1D *fmucandstarty[9];
TH1D *fmucandstartz[9];
TH1D *fmucandendx[9];
TH1D *fmucandendy[9];
TH1D *fmucandendz[9];

TH1D *fpcandstartx[9];
TH1D *fpcandstarty[9];
TH1D *fpcandstartz[9];
TH1D *fpcandendx[9];
TH1D *fpcandendy[9];
TH1D *fpcandendz[9];

TH1D *fEvis_all[9];
TH1D *fQ2cal_all[9];
//===========================================================
TH1D *fvertexx_sig[9];
TH1D *fvertexy_sig[9];
TH1D *fvertexz_sig[9];

TH1D *fmucandstartx_sig[9];
TH1D *fmucandstarty_sig[9];
TH1D *fmucandstartz_sig[9];
TH1D *fmucandendx_sig[9];
TH1D *fmucandendy_sig[9];
TH1D *fmucandendz_sig[9];

TH1D *fpcandstartx_sig[9];
TH1D *fpcandstarty_sig[9];
TH1D *fpcandstartz_sig[9];
TH1D *fpcandendx_sig[9];
TH1D *fpcandendy_sig[9];
TH1D *fpcandendz_sig[9];

TH1D *fEvis_all_sig[9];
TH1D *fQ2cal_all_sig[9];


//===========================================================



TFile *fHistFile;

double FVx = 256.35;
double FVy = 233;
double FVz = 1036.8;
double borderx = 10.;
double bordery = 20.;
double borderz = 10.;

//This function returns if a 3D point is within the fiducial volume
bool inFV(double x, double y, double z) {
    if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
    else return false;
}


double thetax(double theta,double phi){
  TVector3 v;
  v.SetMagThetaPhi(1,theta,phi);
  TVector3 x_axis(1,0,0);
  double theta_x = v.Angle(x_axis);
  return theta_x;
}



void hanalysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   TString option = GetOption();


   //printf("Starting hanalysis with process option:");
   //output root file for histograms
   fHistFile= new TFile("histAnalysis.root","RECREATE");
   //histogram definitions
   
   muon_lenvsdqdx= new TH2D("muonlenvsdqdx", "muonlenvsdqdx", 100, 0, 800, 100, 0, 800);
   proton_lenvsdqdx = new TH2D("protonlenvsdqdx", "protonlenvsdqdx", 100, 0, 1800, 100, 0, 160);

   Pmutruevsreco=new TH2D("Pmutruevsreco", "Pmutruevsreco", 100, 0.0, 3.0,100, 0.0, 3.0);
   Pptruevsreco=new TH2D("Pptruevsreco", "Pptruevsreco", 100, 0.0, 1.5,100, 0.0, 1.5);

   Costhetamutruevsreco=new TH2D("Costhetamutruevsreco", "Costhetamutruevsreco", 100, -1.0, 1.0,100, -1.0, 1.0);
   Cosphimutruevsreco=new TH2D("Cosphimutruevsreco", "Cosphimutruevsreco", 100, -3.14, 3.14,100, -3.14, 3.14);
   
   Costhetaptruevsreco=new TH2D("Costhetaptruevsreco", "Costhetaptruevsreco", 100, -1.0, 1.0,100, -1.0, 1.0);
   Cosphiptruevsreco=new TH2D("Cosphiptruevsreco", "Cosphiptruevsreco", 100, -3.14, 3.14,100, -3.14, 3.14);
     
   trklen_cand= new TH1D("trklen_cand", "trklen_cand", 100, 0.0, 800.0);
   trklen_pcand= new TH1D("trklen_pcand", "trklen_pcand", 500, 0.0, 100.0);


   trunmean_cand= new TH1D("trunmean_cand", "trunmean_cand", 100, 0.0, 800.0);
   trunmean_pcand= new TH1D("trunmean_pcand", "trunmean_pcand", 100, 0.0, 1800.0);
   
   thetalep=new TH1D("thetalep", "thetalep", 100, 0.0, 3.14);
   thetahad=new TH1D("thetahad", "thetahad", 100, 0.0, 3.14);

   thetalep_muon=new TH1D("thetalep_muon", "thetalep_muon", 100, -1.0, 1.0);
   thetahad_proton=new TH1D("thetahad_proton", "thetahad_proton", 100, 0.0, 3.14);


   costhetalep=new TH1D("costhetalep", "costhetalep", 100, -1.0, 1.0);
   costhetahad=new TH1D("costhetahad", "costhetahad", 100, -1.0, 1.0);
   costhetatruelep=new TH1D("costhetatruelep", "costhetatruelep", 100, -1.0, 1.0);
   costhetatruehad=new TH1D("costhetatruehad", "costhetatruehad", 100, -1.0, 1.0);

   costhetaresolep=new TH1D("costhetaresolep", "costhetaresolep", 100, -0.5, 0.5);
   costhetaresohad=new TH1D("costhetaresohad", "costhetaresohad", 100, -0.5, 0.5);



   philep= new TH1D("philep", "philep", 100, -3.14, 3.14);
   phihad= new TH1D("phihad", "phihad", 100, -3.14, 3.14);
   cosphiresolep=new TH1D("cosphiresolep", "cosphiresolep", 100, -1.0, 1.0);
   cosphiresohad=new TH1D("cosphiresohad", "cosphiresohad", 100, -1.0, 1.0);


   philep_muon= new TH1D("philep_muon", "philep_muon", 100, 0.0, 3.14);
   phihad_proton= new TH1D("phihad_proton", "phihad_proton", 100, 0.0, 3.14);

   costheta_vs_trklen= new TH2D("costheta_vs_trklen", "costheta_vs_trklen", 100, -1.0, 1.0, 100, 0., 800);

   ntrkspreco=new TH1D("ntrkspreco", "ntrkspreco", 6.0, -0.5, 5.5);
   trkmomlep= new TH1D("trkmomlep", "trkmomlep", 100, 0.0, 3.0);

   Plep= new TH1D("Plep", "Plep", 100, 0.0, 3.0);
   Phad= new TH1D("Phad", "Phad", 100, 0.0, 1.5);
   Presolep= new TH1D("Presolep", "Presolep", 100, -1.0, 1.0);
   Presohad= new TH1D("Presohad", "Presohad", 100, -1.0, 1.0);
   Plep_muon= new TH1D("Plep_muon", "Plep_muon", 100, 0.0, 3.0);
   Phad_proton= new TH1D("Phad_proton", "Phad_proton", 100, 0.0, 1.5);
 


    h_thetamup= new TH1D("h_thetamup", "h_thetamup", 100, 0.0, 3.14);
    h_phimup= new TH1D("h_phimup", "h_phimup", 100, 0.0, 3.14);

    h_thetapp= new TH1D("h_thetapp", "h_thetapp", 100, 0.0, 3.14);
    h_phipp= new TH1D("h_phipp", "h_phipp", 100, 0.0, 3.14);


    Nhitslep= new TH1D("Nhitslep", "Nhitslep", 100, 0.0,4000.0);
    Nhitshad= new TH1D("Nhitshad", "Nhitshad", 100, 0.0, 600.0); 

    ntrksp_truevsreco=new TH2D("ntrksp_truevsreco", "ntrksp_truevsreco", 3, 0.5, 3.5, 3, 0.5, 3.5);

    fvertex_x=new TH1D("fvertex_x", "fvertex_x", 100, 0.0, 260.0);
    fvertex_y=new TH1D("fvertex_y", "fvertex_y", 100, -120, 120.0);;
    fvertex_z=new TH1D("fvertex_z", "fvertex_z", 100, 0.0, 1200.0);;

    fmucand_startx= new TH1D("fmucand_startx", "fmucand_startx", 100, 0.0, 260.0);
    fmucand_starty= new TH1D("fmucand_starty", "fmucand_starty", 100, -120, 120.0);
    fmucand_startz= new TH1D("fmucand_startz", "fmucand_startz", 100, 0.0, 1200.0);
    fmucand_endx= new TH1D("fmucand_endx", "fmucand_endx", 100, 0.0, 260.0);
    fmucand_endy= new TH1D("fmucand_endy", "fmucand_endy", 100, -120, 120.0);
    fmucand_endz= new TH1D("fmucand_endz", "fmucand_endz", 100, 0.0, 1200.0);

    fpcand_startx= new TH1D("fpcand_startx", "fpcand_startx", 100, 0.0, 260.0);
    fpcand_starty= new TH1D("fpcand_starty", "fpcand_starty", 100, -120, 120.0);
    fpcand_startz= new TH1D("fpcand_startz", "fpcand_startz", 100, 0.0, 1200.0);
    fpcand_endx= new TH1D("fpcand_endx", "fpcand_endx", 100, 0.0, 260.0);
    fpcand_endy= new TH1D("fpcand_endy", "fpcand_endy", 100, -120, 120.0);
    fpcand_endz= new TH1D("fpcand_endz", "fpcand_endz", 100, 0.0, 1200.0);

    fEvis=new TH1D("fEvis", "fEvis", 100, 0.0, 3.0);
    fQ2cal=new TH1D("fQ2cal", "fQ2cal", 100, 0.0, 3.0);



   // this is for MC^^^^^^^^^^^^^^^^^^^^
   for(int i=0; i<9; i++){
     trklen_muon[i]=new TH1D(Form("trklen_muon_%d", i),  Form("trklen_muon_%d", i), 100, 0.0, 800.0);
     trklen_proton[i]=new TH1D(Form("trklen_proton_%d", i), Form("trklen_proton_%d", i),  500, 0.0, 100.0);

     costheta_muon[i]=new TH1D(Form("costheta_muon_%d", i),  Form("costheta_muon_%d", i), 100, -1.0, 1.0);
     costheta_proton[i]=new TH1D(Form("costheta_proton_%d", i), Form("costheta_proton_%d", i),  100, -1.0, 1.0);
     costhetatrue_muon[i]=new TH1D(Form("costhetatrue_muon_%d", i),  Form("costhetatrue_muon_%d", i), 100, -1.0, 1.0);
     costhetatrue_proton[i]=new TH1D(Form("costhetatrue_proton_%d", i), Form("costhetatrue_proton_%d", i),  100, -1.0, 1.0);

     costhetareso_muon[i]=new TH1D(Form("costhetareso_muon_%d", i),  Form("costhetareso_muon_%d", i), 100, -0.5, 0.5);
     costhetareso_proton[i]=new TH1D(Form("costhetareso_proton_%d", i), Form("costhetareso_proton_%d", i),  100, -0.5, 0.5);

     cosphireso_muon[i]=new TH1D(Form("cosphireso_muon_%d", i),  Form("cosphireso_muon_%d", i), 100, -1.0, 1.0);
     cosphireso_proton[i]=new TH1D(Form("cosphireso_proton_%d", i), Form("cosphireso_proton_%d", i),  100, -1.0, 1.0);



     trunmean_muon[i]= new TH1D(Form("trunmean_muon_%d", i), Form("trunmean_muon_%d", i), 100, 0.0, 800.0);
     trunmean_proton[i]= new TH1D(Form("trunmean_proton_%d", i), Form("trunmean_proton_%d", i), 100, 0.0, 1800.0);
     phi_muon[i]=new TH1D(Form("phi_muon_%d", i), Form("phi_muon_%d", i), 100, -3.14, 3.14);
     phi_proton[i]=new TH1D(Form("phi_proton_%d", i), Form("phi_proton_%d", i), 100, -3.14, 3.14);
     ntrks_proton[i]=new TH1D(Form("ntrks_proton_%d", i), Form("ntrks_proton_%d", i), 6.0, -0.5, 5.5);
     trkmom_muon[i]= new TH1D(Form("trkmom_muon_%d", i), Form("trkmom_muon_%d", i), 100.0, 0.0, 3.0);       
     fPmuon[i]= new TH1D(Form("fPmuon_%d", i), Form("fPmuon_%d", i), 100.0, 0.0, 3.0);       
     fPproton[i]= new TH1D(Form("fPproton_%d", i), Form("fPproton_%d", i), 100.0, 0.0, 1.5);       
     fPresomuon[i]= new TH1D(Form("fPresomuon_%d", i), Form("fPresomuon_%d", i), 100.0, -1.0, 1.0);       
     fPresoproton[i]= new TH1D(Form("fPresoproton_%d", i), Form("fPresoproton_%d", i), 100.0, -1.0, 1.0);       
     h_thetamups[i]= new TH1D(Form("h_thetamup_%d", i), Form("h_thetamup_%d", i), 100.0, 0.0, 3.14);
     h_phimups[i]= new TH1D(Form("h_phimup_%d", i), Form("h_phimup_%d", i), 100.0, 0.0, 3.14);

     h_thetapps[i]= new TH1D(Form("h_thetapp_%d", i), Form("h_thetapp_%d", i), 100.0, 0.0, 3.14);
     h_phipps[i]= new TH1D(Form("h_phipp_%d", i), Form("h_phipp_%d", i), 100.0, 0.0, 3.14);


     Nhitsmuon[i]= new TH1D(Form("Nhitsmuon_%d", i), Form("Nhitsmuon_%d", i), 100.0, 0.0, 4000.0);
     Nhitsproton[i]= new TH1D(Form("Nhitsproton_%d", i), Form("Nhitsproton_%d", i), 100.0, 0.0, 600.0);

     fvertexx[i]=new TH1D(Form("fvertexx_%d",i), Form("fvertexx_%d", i), 100, 0.0, 260.0);
     fvertexy[i]=new TH1D(Form("fvertexy_%d",i), Form("fvertexy_%d", i), 100, -120, 120.0);;
     fvertexz[i]=new TH1D(Form("fvertexz_%d",i), Form("fvertexz_%d", i), 100, 0.0, 1200.0);;

     fmucandstartx[i]= new TH1D(Form("fmucandstartx_%d",i), Form("fmucandstartx_%d",i), 100, 0.0, 260.0);
     fmucandstarty[i]= new TH1D(Form("fmucandstarty_%d",i), Form("fmucandstarty_%d",i), 100, -120, 120.0);
     fmucandstartz[i]= new TH1D(Form("fmucandstartz_%d",i), Form("fmucandstartz_%d",i), 100, 0.0, 1200.0);
     fmucandendx[i]= new TH1D(Form("fmucandendx_%d",i), Form("fmucandendx_%d",i), 100, 0.0, 260.0);
     fmucandendy[i]= new TH1D(Form("fmucandendy_%d",i), Form("fmucandendy_%d",i), 100, -120, 120.0);
     fmucandendz[i]= new TH1D(Form("fmucandendz_%d",i), Form("fmucandendz_%d",i), 100, 0.0, 1200.0);

     fpcandstartx[i]= new TH1D(Form("fpcandstartx_%d",i), Form("fpcandstartx_%d",i), 100, 0.0, 260.0);
     fpcandstarty[i]= new TH1D(Form("fpcandstarty_%d",i), Form("fpcandstarty_%d",i), 100, -120, 120.0);
     fpcandstartz[i]= new TH1D(Form("fpcandstartz_%d",i), Form("fpcandstartz_%d",i), 100, 0.0, 1200.0);
     fpcandendx[i]= new TH1D(Form("fpcandendx_%d",i), Form("fpcandendx_%d",i), 100, 0.0, 260.0);
     fpcandendy[i]= new TH1D(Form("fpcandendy_%d",i), Form("fpcandendy_%d",i), 100, -120, 120.0);
     fpcandendz[i]= new TH1D(Form("fpcandendz_%d",i), Form("fpcandendz_%d",i), 100, 0.0, 1200.0);

     fEvis_all[i]=new TH1D(Form("fEvis_all_%d",i), Form("fEvis_all_%d",i), 100, 0.0, 3.0);
     fQ2cal_all[i]=new TH1D(Form("fQ2cal_all_%d",i), Form("fQ2cal_all_%d",i), 100, 0.0, 3.0);


     trklen_muon_sig[i]=new TH1D(Form("trklen_muon_sig_%d", i),  Form("trklen_muon_sig_%d", i), 100, 0.0, 800.0);
     trklen_proton_sig[i]=new TH1D(Form("trklen_proton_sig_%d", i), Form("trklen_proton_sig_%d", i),  500, 0.0, 100.0);

     costheta_muon_sig[i]=new TH1D(Form("costheta_muon_sig_%d", i),  Form("costheta_muon_sig_%d", i), 100, -1.0, 1.0);
     costheta_proton_sig[i]=new TH1D(Form("costheta_proton_sig_%d", i), Form("costheta_proton_sig_%d", i),  100, -1.0, 1.0);
     costhetatrue_muon_sig[i]=new TH1D(Form("costhetatrue_muon_sig_%d", i),  Form("costhetatrue_muon_sig_%d", i), 100, -1.0, 1.0);
     costhetatrue_proton_sig[i]=new TH1D(Form("costhetatrue_proton_sig_%d", i), Form("costhetatrue_proton_sig_%d", i),  100, -1.0, 1.0);

     costhetareso_muon_sig[i]=new TH1D(Form("costhetareso_muon_sig_%d", i),  Form("costhetareso_muon_sig_%d", i), 100, -0.5, 0.5);
     costhetareso_proton_sig[i]=new TH1D(Form("costhetareso_proton_sig_%d", i), Form("costhetareso_proton_sig_%d", i),  100, -0.5, 0.5);

     cosphireso_muon_sig[i]=new TH1D(Form("cosphireso_muon_sig_%d", i),  Form("cosphireso_muon_sig_%d", i), 100, -1.0, 1.0);
     cosphireso_proton_sig[i]=new TH1D(Form("cosphireso_proton_sig_%d", i), Form("cosphireso_proton_sig_%d", i),  100, -1.0, 1.0);



     trunmean_muon_sig[i]= new TH1D(Form("trunmean_muon_sig_%d", i), Form("trunmean_muon_sig_%d", i), 100, 0.0, 800.0);
     trunmean_proton_sig[i]= new TH1D(Form("trunmean_proton_sig_%d", i), Form("trunmean_proton_sig_%d", i), 100, 0.0, 1800.0);
     phi_muon_sig[i]=new TH1D(Form("phi_muon_sig_%d", i), Form("phi_muon_sig_%d", i), 100, -3.14, 3.14);
     phi_proton_sig[i]=new TH1D(Form("phi_proton_sig_%d", i), Form("phi_proton_sig_%d", i), 100, -3.14, 3.14);
     ntrks_proton_sig[i]=new TH1D(Form("ntrks_proton_sig_%d", i), Form("ntrks_proton_sig_%d", i), 6.0, -0.5, 5.5);
     trkmom_muon_sig[i]= new TH1D(Form("trkmom_muon_sig_%d", i), Form("trkmom_muon_sig_%d", i), 100.0, 0.0, 3.0);       
     fPmuon_sig[i]= new TH1D(Form("fPmuon_sig_%d", i), Form("fPmuon_sig_%d", i), 100.0, 0.0, 3.0);       
     fPproton_sig[i]= new TH1D(Form("fPproton_sig_%d", i), Form("fPproton_sig_%d", i), 100.0, 0.0, 1.5);       
     fPresomuon_sig[i]= new TH1D(Form("fPresomuon_sig_%d", i), Form("fPresomuon_sig_%d", i), 100.0, -1.0, 1.0);       
     fPresoproton_sig[i]= new TH1D(Form("fPresoproton_sig_%d", i), Form("fPresoproton_sig_%d", i), 100.0, -1.0, 1.0);       
     h_thetamups_sig[i]= new TH1D(Form("h_thetamup_sig_%d", i), Form("h_thetamup_sig_%d", i), 100.0, 0.0, 3.14);
     h_phimups_sig[i]= new TH1D(Form("h_phimup_sig_%d", i), Form("h_phimup_sig_%d", i), 100.0, 0.0, 3.14);

     h_thetapps_sig[i]= new TH1D(Form("h_thetapp_sig_%d", i), Form("h_thetapp_sig_%d", i), 100.0, 0.0, 3.14);
     h_phipps_sig[i]= new TH1D(Form("h_phipp_sig_%d", i), Form("h_phipp_sig_%d", i), 100.0, 0.0, 3.14);


     Nhitsmuon_sig[i]= new TH1D(Form("Nhitsmuon_sig_%d", i), Form("Nhitsmuon_sig_%d", i), 100.0, 0.0, 4000.0);
     Nhitsproton_sig[i]= new TH1D(Form("Nhitsproton_sig_%d", i), Form("Nhitsproton_sig_%d", i), 100.0, 0.0, 600.0);

     fvertexx_sig[i]=new TH1D(Form("fvertexx_sig_%d",i), Form("fvertexx_sig_%d", i), 100, 0.0, 260.0);
     fvertexy_sig[i]=new TH1D(Form("fvertexy_sig_%d",i), Form("fvertexy_sig_%d", i), 100, -120, 120.0);;
     fvertexz_sig[i]=new TH1D(Form("fvertexz_sig_%d",i), Form("fvertexz_sig_%d", i), 100, 0.0, 1200.0);;

     fmucandstartx_sig[i]= new TH1D(Form("fmucandstartx_sig_%d",i), Form("fmucandstartx_sig_%d",i), 100, 0.0, 260.0);
     fmucandstarty_sig[i]= new TH1D(Form("fmucandstarty_sig_%d",i), Form("fmucandstarty_sig_%d",i), 100, -120, 120.0);
     fmucandstartz_sig[i]= new TH1D(Form("fmucandstartz_sig_%d",i), Form("fmucandstartz_sig_%d",i), 100, 0.0, 1200.0);
     fmucandendx_sig[i]= new TH1D(Form("fmucandendx_sig_%d",i), Form("fmucandendx_sig_%d",i), 100, 0.0, 260.0);
     fmucandendy_sig[i]= new TH1D(Form("fmucandendy_sig_%d",i), Form("fmucandendy_sig_%d",i), 100, -120, 120.0);
     fmucandendz_sig[i]= new TH1D(Form("fmucandendz_sig_%d",i), Form("fmucandendz_sig_%d",i), 100, 0.0, 1200.0);

     fpcandstartx_sig[i]= new TH1D(Form("fpcandstartx_sig_%d",i), Form("fpcandstartx_sig_%d",i), 100, 0.0, 260.0);
     fpcandstarty_sig[i]= new TH1D(Form("fpcandstarty_sig_%d",i), Form("fpcandstarty_sig_%d",i), 100, -120, 120.0);
     fpcandstartz_sig[i]= new TH1D(Form("fpcandstartz_sig_%d",i), Form("fpcandstartz_sig_%d",i), 100, 0.0, 1200.0);
     fpcandendx_sig[i]= new TH1D(Form("fpcandendx_sig_%d",i), Form("fpcandendx_sig_%d",i), 100, 0.0, 260.0);
     fpcandendy_sig[i]= new TH1D(Form("fpcandendy_sig_%d",i), Form("fpcandendy_sig_%d",i), 100, -120, 120.0);
     fpcandendz_sig[i]= new TH1D(Form("fpcandendz_sig_%d",i), Form("fpcandendz_sig_%d",i), 100, 0.0, 1200.0);

     fEvis_all_sig[i]=new TH1D(Form("fEvis_all_sig_%d",i), Form("fEvis_all_sig_%d",i), 100, 0.0, 3.0);
     fQ2cal_all_sig[i]=new TH1D(Form("fQ2cal_all_sig_%d",i), Form("fQ2cal_all_sig_%d",i), 100, 0.0, 3.0);




    }
   // this is for MC^^^^^^^^^^^^^^^^^^^^^^^^^^
  

}

void hanalysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t hanalysis::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);
   //if(*_fTruemode  !=10)
   //if(*_fTruemode != 10 && *fPhiLep> 1.0 && *fPhiLep<2.14 && *fPhiHad<-1.0 && *fPhiHad>-2.14)
   //if(*_fTruemode != 10 && (*fPhiLep< 1.0 || *fPhiLep>2.14) && (*fPhiHad>-1.0 || *fPhiHad<-2.14))
   //if(inFV(*trackendxcandidate,*trackendycandidate,*trackendzcandidate) &&inFV(*trackstartxcandidate,*trackstartycandidate,*trackstartzcandidate)) 
   //if(!inFV(*trackendxcandidate,*trackendycandidate,*trackendzcandidate) || !inFV(*trackstartxcandidate,*trackstartycandidate,*trackstartzcandidate)) 
   
   //cosmic-removal cut
   if (*trackendycandidate > 95 && *fPhiLep>0){return true;}
   // DIC cut(s)
   if (TMath::Cos(thetax(*fThetaLep, *fPhiLep)) > 0.8){return true;}
   if (TMath::Cos(thetax(*fThetaLep, *fPhiLep)) < -0.8){return true;}

   {

   fvertex_x->Fill(*fvtxx);
   fvertex_y->Fill(*fvtxy);
   fvertex_z->Fill(*fvtxz);;

   fmucand_startx->Fill(*trackstartxcandidate);
   fmucand_starty->Fill(*trackstartycandidate);
   fmucand_startz->Fill(*trackstartzcandidate);
   fmucand_endx->Fill(*trackendxcandidate);
   fmucand_endy->Fill(*trackendycandidate);
   fmucand_endz->Fill(*trackendzcandidate);

   //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   // loop over all the proton candidates and fill the histograms
   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcand_startx->Fill(protoncandidate_startx[allp]);
      fpcand_starty->Fill(protoncandidate_starty[allp]);
      fpcand_startz->Fill(protoncandidate_startz[allp]);
      fpcand_endx->Fill(protoncandidate_endx[allp]);
      fpcand_endy->Fill(protoncandidate_endy[allp]);
      fpcand_endz->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapp->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipp->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
   //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   fEvis->Fill(*Evis);
   fQ2cal->Fill(*Q2cal);





    
   muon_lenvsdqdx->Fill(*TrunMean_cand, *ftrklenmuoncand);
   proton_lenvsdqdx->Fill(*TrunMean_pcand, *ftrklenprotoncand);
   

   //calculate angle between muon and proton
   TVector3 muonP(*fPlep*(TMath::Sin(*fThetaLep))*(TMath::Cos(*fPhiLep)), *fPlep*(TMath::Sin(*fThetaLep))*(TMath::Sin(*fPhiLep)), *fPlep*(TMath::Cos(*fThetaLep)));
   TVector3 protonP(*fPhad*(TMath::Sin(*fThetaHad))*(TMath::Cos(*fPhiHad)), *fPhad*(TMath::Sin(*fThetaHad))*(TMath::Sin(*fPhiHad)), *fPhad*(TMath::Cos(*fThetaHad)));



   h_thetamup->Fill(muonP.Angle(protonP));
   //h_thetamup->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
   h_phimup->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

   //if(muonP.Angle(protonP)<2.5)  // cosmic study
   {

   trunmean_cand->Fill(*TrunMean_cand*198.0/198.0, 1.0);
   trunmean_pcand->Fill(*TrunMean_pcand*198.0/198.0, 1.0);


   trklen_cand->Fill(*ftrklenmuoncand);
   trklen_pcand->Fill(*ftrklenprotoncand);


   thetalep->Fill(*fThetaLep);
   thetahad->Fill(*fThetaHad);


   //std::cout<<TMath::Cos(*fThetaLep)<<std::endl;
   philep->Fill(*fPhiLep);
   phihad->Fill(*fPhiHad); 

   costhetalep->Fill(TMath::Cos(*fThetaLep));
   costhetahad->Fill(TMath::Cos(*fThetaHad));

   costhetatruelep->Fill(*trackcand_parCosTheta);
   costhetatruehad->Fill(*trackpcand_parCosTheta);

   costhetaresolep->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
   costhetaresohad->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

   Costhetamutruevsreco->Fill(*trackcand_parCosTheta,TMath::Cos(*fThetaLep));
   Costhetaptruevsreco->Fill(*trackpcand_parCosTheta, TMath::Cos(*fThetaHad));


   cosphiresolep->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
   cosphiresohad->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);

   Cosphimutruevsreco->Fill(*trackcand_parPhi, *fPhiLep);
   Cosphiptruevsreco->Fill(*trackpcand_parPhi, *fPhiHad);

   costheta_vs_trklen->Fill(TMath::Cos(*fThetaLep), *ftrklenmuoncand);




   Plep->Fill(*fPlep);
   Phad->Fill(*fPhad);

   Presolep->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
   Presohad->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
   //fill the histograms with signal
   if(*OOFVflag==0 && (*truthtop==2 ||*truthtop==3 ||*truthtop==4)){
     if(*trackcand_origin!=1 ||*trackpcand_origin!=1){
      std::cout<<"origin of the muon cand is "<<*trackcand_origin<<" origin of proton cand is"<<*trackpcand_origin<<std::endl;
     }
     if(*trackcand_origin==1 && *trackpcand_origin==1){

      if((abs(*trackcand_parPDG)==13 && abs(*trackpcand_parPDG)==2212) ||(abs(*trackcand_parPDG)==2212 && abs(*trackpcand_parPDG)==13)){
       if(abs(*trackcand_parPDG)==13 && abs(*trackpcand_parPDG)==2212){
         Plep_muon->Fill(sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
         Phad_proton->Fill(sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
         thetalep_muon->Fill(TMath::Cos(*fThetaLep));
         thetahad_proton->Fill(*fThetaHad);
         philep_muon->Fill(*fPhiLep);
         phihad_proton->Fill(*fPhiHad);
       }
       if(abs(*trackcand_parPDG)==2212 && abs(*trackpcand_parPDG)==13){
         Phad_proton->Fill(sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
         Plep_muon->Fill(sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
         thetahad_proton->Fill(*fThetaLep);
         thetalep_muon->Fill(TMath::Cos(*fThetaHad));
         phihad_proton->Fill(*fPhiLep);
         philep_muon->Fill(*fPhiHad);
       }
      }
     } 
   }



   Pmutruevsreco->Fill(sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass),*fPlep);
   Pptruevsreco->Fill(sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass),*fPhad);

   //get the number of true protons and plot the Ntruep vs NRecoP
   ntrksp_truevsreco->Fill(*fNTruePTrks, *fNRecoPTrks);
  
   ntrkspreco->Fill(*fNRecoPTrks);

   trkmomlep->Fill(*trackmomcandidate_mcs);

   Nhitslep->Fill(*Nhits_muoncand);
   Nhitshad->Fill(*Nhits_protoncand);

   //of all the track candidates, check if the tracks is true muon/proton and get the momentum and angles




   if(*OOFVflag==0 && (*truthtop==2 || *truthtop==3 || *truthtop==4) && (*trackcand_origin==1 && *trackpcand_origin==1)){
      trklen_muon[0]->Fill(*ftrklenmuoncand);
      trklen_proton[0]->Fill(*ftrklenprotoncand);
      costheta_muon[0]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[0]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[0]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[0]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[0]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[0]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[0]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[0]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);




      trunmean_muon[0]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[0]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[0]->Fill(*fPhiLep);
      phi_proton[0]->Fill(*fPhiHad);
      ntrks_proton[0]->Fill(*fNRecoPTrks);
      trkmom_muon[0]->Fill(*trackmomcandidate_mcs);
      fPmuon[0]->Fill(*fPlep);
      fPproton[0]->Fill(*fPhad);
      fPresomuon[0]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[0]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[0]->Fill(muonP.Angle(protonP));
      //h_thetamups[0]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[0]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));
      Nhitsmuon[0]->Fill(*Nhits_muoncand);
      Nhitsproton[0]->Fill(*Nhits_protoncand);

      fvertexx[0]->Fill(*fvtxx);
      fvertexy[0]->Fill(*fvtxy);
      fvertexz[0]->Fill(*fvtxz);;

      fmucandstartx[0]->Fill(*trackstartxcandidate);
      fmucandstarty[0]->Fill(*trackstartycandidate);
      fmucandstartz[0]->Fill(*trackstartzcandidate);
      fmucandendx[0]->Fill(*trackendxcandidate);
      fmucandendy[0]->Fill(*trackendycandidate);
      fmucandendz[0]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[0]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[0]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[0]->Fill(protoncandidate_startz[allp]);
      fpcandendx[0]->Fill(protoncandidate_endx[allp]);
      fpcandendy[0]->Fill(protoncandidate_endy[allp]);
      fpcandendz[0]->Fill(protoncandidate_endz[allp]);
   }

   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[0]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[0]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[0]->Fill(*Evis);
      fQ2cal_all[0]->Fill(*Q2cal);
    //get all the signals with interaction separated==============================================
    if(*_fTruemode==0){
      trklen_muon_sig[0]->Fill(*ftrklenmuoncand);
      trklen_proton_sig[0]->Fill(*ftrklenprotoncand);
      costheta_muon_sig[0]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton_sig[0]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon_sig[0]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton_sig[0]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon_sig[0]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton_sig[0]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon_sig[0]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton_sig[0]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);




      trunmean_muon_sig[0]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton_sig[0]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon_sig[0]->Fill(*fPhiLep);
      phi_proton_sig[0]->Fill(*fPhiHad);
      ntrks_proton_sig[0]->Fill(*fNRecoPTrks);
      trkmom_muon_sig[0]->Fill(*trackmomcandidate_mcs);
      fPmuon_sig[0]->Fill(*fPlep);
      fPproton_sig[0]->Fill(*fPhad);
      fPresomuon_sig[0]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton_sig[0]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups_sig[0]->Fill(muonP.Angle(protonP));
      //h_thetamups_sig[0]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups_sig[0]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));
      Nhitsmuon_sig[0]->Fill(*Nhits_muoncand);
      Nhitsproton_sig[0]->Fill(*Nhits_protoncand);

      fvertexx_sig[0]->Fill(*fvtxx);
      fvertexy_sig[0]->Fill(*fvtxy);
      fvertexz_sig[0]->Fill(*fvtxz);;

      fmucandstartx_sig[0]->Fill(*trackstartxcandidate);
      fmucandstarty_sig[0]->Fill(*trackstartycandidate);
      fmucandstartz_sig[0]->Fill(*trackstartzcandidate);
      fmucandendx_sig[0]->Fill(*trackendxcandidate);
      fmucandendy_sig[0]->Fill(*trackendycandidate);
      fmucandendz_sig[0]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx_sig[0]->Fill(protoncandidate_startx[allp]);
      fpcandstarty_sig[0]->Fill(protoncandidate_starty[allp]);
      fpcandstartz_sig[0]->Fill(protoncandidate_startz[allp]);
      fpcandendx_sig[0]->Fill(protoncandidate_endx[allp]);
      fpcandendy_sig[0]->Fill(protoncandidate_endy[allp]);
      fpcandendz_sig[0]->Fill(protoncandidate_endz[allp]);
   }

   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps_sig[0]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps_sig[0]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }

      fEvis_all_sig[0]->Fill(*Evis);
      fQ2cal_all_sig[0]->Fill(*Q2cal);

    } else if(*_fTruemode==1){
      trklen_muon_sig[1]->Fill(*ftrklenmuoncand);
      trklen_proton_sig[1]->Fill(*ftrklenprotoncand);
      costheta_muon_sig[1]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton_sig[1]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon_sig[1]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton_sig[1]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon_sig[1]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton_sig[1]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon_sig[1]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton_sig[1]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);




      trunmean_muon_sig[1]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton_sig[1]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon_sig[1]->Fill(*fPhiLep);
      phi_proton_sig[1]->Fill(*fPhiHad);
      ntrks_proton_sig[1]->Fill(*fNRecoPTrks);
      trkmom_muon_sig[1]->Fill(*trackmomcandidate_mcs);
      fPmuon_sig[1]->Fill(*fPlep);
      fPproton_sig[1]->Fill(*fPhad);
      fPresomuon_sig[1]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton_sig[1]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups_sig[1]->Fill(muonP.Angle(protonP));
      //h_thetamups_sig[0]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups_sig[1]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));
      Nhitsmuon_sig[1]->Fill(*Nhits_muoncand);
      Nhitsproton_sig[1]->Fill(*Nhits_protoncand);

      fvertexx_sig[1]->Fill(*fvtxx);
      fvertexy_sig[1]->Fill(*fvtxy);
      fvertexz_sig[1]->Fill(*fvtxz);;

      fmucandstartx_sig[1]->Fill(*trackstartxcandidate);
      fmucandstarty_sig[1]->Fill(*trackstartycandidate);
      fmucandstartz_sig[1]->Fill(*trackstartzcandidate);
      fmucandendx_sig[1]->Fill(*trackendxcandidate);
      fmucandendy_sig[1]->Fill(*trackendycandidate);
      fmucandendz_sig[1]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx_sig[1]->Fill(protoncandidate_startx[allp]);
      fpcandstarty_sig[1]->Fill(protoncandidate_starty[allp]);
      fpcandstartz_sig[1]->Fill(protoncandidate_startz[allp]);
      fpcandendx_sig[1]->Fill(protoncandidate_endx[allp]);
      fpcandendy_sig[1]->Fill(protoncandidate_endy[allp]);
      fpcandendz_sig[1]->Fill(protoncandidate_endz[allp]);
   }

   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps_sig[1]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps_sig[1]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      fEvis_all_sig[1]->Fill(*Evis);
      fQ2cal_all_sig[1]->Fill(*Q2cal);

    } else if(*_fTruemode==2){
      trklen_muon_sig[2]->Fill(*ftrklenmuoncand);
      trklen_proton_sig[2]->Fill(*ftrklenprotoncand);
      costheta_muon_sig[2]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton_sig[2]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon_sig[2]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton_sig[2]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon_sig[2]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton_sig[2]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon_sig[2]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton_sig[2]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);




      trunmean_muon_sig[2]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton_sig[2]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon_sig[2]->Fill(*fPhiLep);
      phi_proton_sig[2]->Fill(*fPhiHad);
      ntrks_proton_sig[2]->Fill(*fNRecoPTrks);
      trkmom_muon_sig[2]->Fill(*trackmomcandidate_mcs);
      fPmuon_sig[2]->Fill(*fPlep);
      fPproton_sig[2]->Fill(*fPhad);
      fPresomuon_sig[2]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton_sig[2]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups_sig[2]->Fill(muonP.Angle(protonP));
      //h_thetamups_sig[0]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups_sig[2]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));
      Nhitsmuon_sig[2]->Fill(*Nhits_muoncand);
      Nhitsproton_sig[2]->Fill(*Nhits_protoncand);

      fvertexx_sig[2]->Fill(*fvtxx);
      fvertexy_sig[2]->Fill(*fvtxy);
      fvertexz_sig[2]->Fill(*fvtxz);;

      fmucandstartx_sig[2]->Fill(*trackstartxcandidate);
      fmucandstarty_sig[2]->Fill(*trackstartycandidate);
      fmucandstartz_sig[2]->Fill(*trackstartzcandidate);
      fmucandendx_sig[2]->Fill(*trackendxcandidate);
      fmucandendy_sig[2]->Fill(*trackendycandidate);
      fmucandendz_sig[2]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx_sig[2]->Fill(protoncandidate_startx[allp]);
      fpcandstarty_sig[2]->Fill(protoncandidate_starty[allp]);
      fpcandstartz_sig[2]->Fill(protoncandidate_startz[allp]);
      fpcandendx_sig[2]->Fill(protoncandidate_endx[allp]);
      fpcandendy_sig[2]->Fill(protoncandidate_endy[allp]);
      fpcandendz_sig[2]->Fill(protoncandidate_endz[allp]);
   }

   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps_sig[2]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps_sig[2]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      
      fEvis_all_sig[2]->Fill(*Evis);
      fQ2cal_all_sig[2]->Fill(*Q2cal);

    } else if(*_fTruemode==10){
      trklen_muon_sig[3]->Fill(*ftrklenmuoncand);
      trklen_proton_sig[3]->Fill(*ftrklenprotoncand);
      costheta_muon_sig[3]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton_sig[3]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon_sig[3]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton_sig[3]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon_sig[3]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton_sig[3]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon_sig[3]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton_sig[3]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);




      trunmean_muon_sig[3]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton_sig[3]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon_sig[3]->Fill(*fPhiLep);
      phi_proton_sig[3]->Fill(*fPhiHad);
      ntrks_proton_sig[3]->Fill(*fNRecoPTrks);
      trkmom_muon_sig[3]->Fill(*trackmomcandidate_mcs);
      fPmuon_sig[3]->Fill(*fPlep);
      fPproton_sig[3]->Fill(*fPhad);
      fPresomuon_sig[3]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton_sig[3]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups_sig[3]->Fill(muonP.Angle(protonP));
      //h_thetamups_sig[0]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups_sig[3]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));
      Nhitsmuon_sig[3]->Fill(*Nhits_muoncand);
      Nhitsproton_sig[3]->Fill(*Nhits_protoncand);

      fvertexx_sig[3]->Fill(*fvtxx);
      fvertexy_sig[3]->Fill(*fvtxy);
      fvertexz_sig[3]->Fill(*fvtxz);;

      fmucandstartx_sig[3]->Fill(*trackstartxcandidate);
      fmucandstarty_sig[3]->Fill(*trackstartycandidate);
      fmucandstartz_sig[3]->Fill(*trackstartzcandidate);
      fmucandendx_sig[3]->Fill(*trackendxcandidate);
      fmucandendy_sig[3]->Fill(*trackendycandidate);
      fmucandendz_sig[3]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx_sig[3]->Fill(protoncandidate_startx[allp]);
      fpcandstarty_sig[3]->Fill(protoncandidate_starty[allp]);
      fpcandstartz_sig[3]->Fill(protoncandidate_startz[allp]);
      fpcandendx_sig[3]->Fill(protoncandidate_endx[allp]);
      fpcandendy_sig[3]->Fill(protoncandidate_endy[allp]);
      fpcandendz_sig[3]->Fill(protoncandidate_endz[allp]);
   }

   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps_sig[3]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps_sig[3]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
      
      fEvis_all_sig[3]->Fill(*Evis);
      fQ2cal_all_sig[3]->Fill(*Q2cal);
 
    }


    //=============================================================================================
   }
   
   if(*OOFVflag==0 && *truthtop==1){
      trklen_muon[1]->Fill(*ftrklenmuoncand);
      trklen_proton[1]->Fill(*ftrklenprotoncand);
      costheta_muon[1]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[1]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[1]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[1]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[1]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[1]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[1]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[1]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);


      trunmean_muon[1]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[1]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[1]->Fill(*fPhiLep);
      phi_proton[1]->Fill(*fPhiHad);
      ntrks_proton[1]->Fill(*fNRecoPTrks);
      trkmom_muon[1]->Fill(*trackmomcandidate_mcs);
      fPmuon[1]->Fill(*fPlep);
      fPproton[1]->Fill(*fPhad);
      fPresomuon[1]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[1]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[1]->Fill(muonP.Angle(protonP)); 
      //h_thetamups[1]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[1]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

      Nhitsmuon[1]->Fill(*Nhits_muoncand);
      Nhitsproton[1]->Fill(*Nhits_protoncand);
      fvertexx[1]->Fill(*fvtxx);
      fvertexy[1]->Fill(*fvtxy);
      fvertexz[1]->Fill(*fvtxz);;

      fmucandstartx[1]->Fill(*trackstartxcandidate);
      fmucandstarty[1]->Fill(*trackstartycandidate);
      fmucandstartz[1]->Fill(*trackstartzcandidate);
      fmucandendx[1]->Fill(*trackendxcandidate);
      fmucandendy[1]->Fill(*trackendycandidate);
      fmucandendz[1]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[1]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[1]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[1]->Fill(protoncandidate_startz[allp]);
      fpcandendx[1]->Fill(protoncandidate_endx[allp]);
      fpcandendy[1]->Fill(protoncandidate_endy[allp]);
      fpcandendz[1]->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[1]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[1]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[1]->Fill(*Evis);
      fQ2cal_all[1]->Fill(*Q2cal);



     }
   if(*OOFVflag==0 && *truthtop==5){
      trklen_muon[2]->Fill(*ftrklenmuoncand);
      trklen_proton[2]->Fill(*ftrklenprotoncand);
      costheta_muon[2]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[2]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[2]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[2]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[2]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[2]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[2]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[2]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);


      trunmean_muon[2]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[2]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[2]->Fill(*fPhiLep);
      phi_proton[2]->Fill(*fPhiHad);
      ntrks_proton[2]->Fill(*fNRecoPTrks);
      trkmom_muon[2]->Fill(*trackmomcandidate_mcs);
      fPmuon[2]->Fill(*fPlep);
      fPproton[2]->Fill(*fPhad);
      fPresomuon[2]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[2]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[2]->Fill(muonP.Angle(protonP)); 
      //h_thetamups[2]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[2]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

      Nhitsmuon[2]->Fill(*Nhits_muoncand);
      Nhitsproton[2]->Fill(*Nhits_protoncand);

      fvertexx[2]->Fill(*fvtxx);
      fvertexy[2]->Fill(*fvtxy);
      fvertexz[2]->Fill(*fvtxz);;

      fmucandstartx[2]->Fill(*trackstartxcandidate);
      fmucandstarty[2]->Fill(*trackstartycandidate);
      fmucandstartz[2]->Fill(*trackstartzcandidate);
      fmucandendx[2]->Fill(*trackendxcandidate);
      fmucandendy[2]->Fill(*trackendycandidate);
      fmucandendz[2]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[2]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[2]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[2]->Fill(protoncandidate_startz[allp]);
      fpcandendx[2]->Fill(protoncandidate_endx[allp]);
      fpcandendy[2]->Fill(protoncandidate_endy[allp]);
      fpcandendz[2]->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[2]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[2]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[2]->Fill(*Evis);
      fQ2cal_all[2]->Fill(*Q2cal);



      }
   if(*OOFVflag==0 && *truthtop==6){
      trklen_muon[3]->Fill(*ftrklenmuoncand);
      trklen_proton[3]->Fill(*ftrklenprotoncand);
      costheta_muon[3]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[3]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[3]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[3]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[3]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[3]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[3]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[3]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);


      trunmean_muon[3]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[3]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[3]->Fill(*fPhiLep);
      phi_proton[3]->Fill(*fPhiHad);
      ntrks_proton[3]->Fill(*fNRecoPTrks);
      trkmom_muon[3]->Fill(*trackmomcandidate_mcs);
      fPmuon[3]->Fill(*fPlep);
      fPproton[3]->Fill(*fPhad);
      fPresomuon[3]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[3]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[3]->Fill(muonP.Angle(protonP)); 
      //h_thetamups[3]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[3]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

      Nhitsmuon[3]->Fill(*Nhits_muoncand);
      Nhitsproton[3]->Fill(*Nhits_protoncand);

      fvertexx[3]->Fill(*fvtxx);
      fvertexy[3]->Fill(*fvtxy);
      fvertexz[3]->Fill(*fvtxz);;

      fmucandstartx[3]->Fill(*trackstartxcandidate);
      fmucandstarty[3]->Fill(*trackstartycandidate);
      fmucandstartz[3]->Fill(*trackstartzcandidate);
      fmucandendx[3]->Fill(*trackendxcandidate);
      fmucandendy[3]->Fill(*trackendycandidate);
      fmucandendz[3]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[3]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[3]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[3]->Fill(protoncandidate_startz[allp]);
      fpcandendx[3]->Fill(protoncandidate_endx[allp]);
      fpcandendy[3]->Fill(protoncandidate_endy[allp]);
      fpcandendz[3]->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[3]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[3]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[3]->Fill(*Evis);
      fQ2cal_all[3]->Fill(*Q2cal);


      }
   if(*OOFVflag==0 && *truthtop==7){
      trklen_muon[4]->Fill(*ftrklenmuoncand);
      trklen_proton[4]->Fill(*ftrklenprotoncand);
      costheta_muon[4]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[4]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[4]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[4]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[4]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[4]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[4]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[4]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);


      trunmean_muon[4]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[4]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[4]->Fill(*fPhiLep);
      phi_proton[4]->Fill(*fPhiHad);
      ntrks_proton[4]->Fill(*fNRecoPTrks);
      trkmom_muon[4]->Fill(*trackmomcandidate_mcs);
      fPmuon[4]->Fill(*fPlep);
      fPproton[4]->Fill(*fPhad);
      fPresomuon[4]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[4]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[4]->Fill(muonP.Angle(protonP)); 
      //h_thetamups[4]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[4]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

      Nhitsmuon[4]->Fill(*Nhits_muoncand);
      Nhitsproton[4]->Fill(*Nhits_protoncand);

      fvertexx[4]->Fill(*fvtxx);
      fvertexy[4]->Fill(*fvtxy);
      fvertexz[4]->Fill(*fvtxz);;

      fmucandstartx[4]->Fill(*trackstartxcandidate);
      fmucandstarty[4]->Fill(*trackstartycandidate);
      fmucandstartz[4]->Fill(*trackstartzcandidate);
      fmucandendx[4]->Fill(*trackendxcandidate);
      fmucandendy[4]->Fill(*trackendycandidate);
      fmucandendz[4]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[4]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[4]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[4]->Fill(protoncandidate_startz[allp]);
      fpcandendx[4]->Fill(protoncandidate_endx[allp]);
      fpcandendy[4]->Fill(protoncandidate_endy[allp]);
      fpcandendz[4]->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[4]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[4]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[4]->Fill(*Evis);
      fQ2cal_all[4]->Fill(*Q2cal);


      }
   if(*OOFVflag==0 && *truthtop==8){
      trklen_muon[5]->Fill(*ftrklenmuoncand);
      trklen_proton[5]->Fill(*ftrklenprotoncand);
      costheta_muon[5]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[5]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[5]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[5]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[5]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[5]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[5]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[5]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);


      trunmean_muon[5]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[5]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[5]->Fill(*fPhiLep);
      phi_proton[5]->Fill(*fPhiHad);
      ntrks_proton[5]->Fill(*fNRecoPTrks);
      trkmom_muon[5]->Fill(*trackmomcandidate_mcs);
      fPmuon[5]->Fill(*fPlep);
      fPproton[5]->Fill(*fPhad);
      fPresomuon[5]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[5]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[5]->Fill(muonP.Angle(protonP)); 
      //h_thetamups[5]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[5]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

      Nhitsmuon[5]->Fill(*Nhits_muoncand);
      Nhitsproton[5]->Fill(*Nhits_protoncand);

      fvertexx[5]->Fill(*fvtxx);
      fvertexy[5]->Fill(*fvtxy);
      fvertexz[5]->Fill(*fvtxz);;

      fmucandstartx[5]->Fill(*trackstartxcandidate);
      fmucandstarty[5]->Fill(*trackstartycandidate);
      fmucandstartz[5]->Fill(*trackstartzcandidate);
      fmucandendx[5]->Fill(*trackendxcandidate);
      fmucandendy[5]->Fill(*trackendycandidate);
      fmucandendz[5]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[5]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[5]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[5]->Fill(protoncandidate_startz[allp]);
      fpcandendx[5]->Fill(protoncandidate_endx[allp]);
      fpcandendy[5]->Fill(protoncandidate_endy[allp]);
      fpcandendz[5]->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[5]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[5]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[5]->Fill(*Evis);
      fQ2cal_all[5]->Fill(*Q2cal);


      }
   if(*OOFVflag==1){
      trklen_muon[6]->Fill(*ftrklenmuoncand);
      trklen_proton[6]->Fill(*ftrklenprotoncand);
      costheta_muon[6]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[6]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[6]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[6]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[6]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[6]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[6]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[6]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);


      trunmean_muon[6]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[6]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[6]->Fill(*fPhiLep);
      phi_proton[6]->Fill(*fPhiHad);
      ntrks_proton[6]->Fill(*fNRecoPTrks);
      trkmom_muon[6]->Fill(*trackmomcandidate_mcs);
      fPmuon[6]->Fill(*fPlep);
      fPproton[6]->Fill(*fPhad);
      fPresomuon[6]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[6]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[6]->Fill(muonP.Angle(protonP)); 
      //h_thetamups[6]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[6]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

      Nhitsmuon[6]->Fill(*Nhits_muoncand);
      Nhitsproton[6]->Fill(*Nhits_protoncand);

      fvertexx[6]->Fill(*fvtxx);
      fvertexy[6]->Fill(*fvtxy);
      fvertexz[6]->Fill(*fvtxz);;

      fmucandstartx[6]->Fill(*trackstartxcandidate);
      fmucandstarty[6]->Fill(*trackstartycandidate);
      fmucandstartz[6]->Fill(*trackstartzcandidate);
      fmucandendx[6]->Fill(*trackendxcandidate);
      fmucandendy[6]->Fill(*trackendycandidate);
      fmucandendz[6]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[6]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[6]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[6]->Fill(protoncandidate_startz[allp]);
      fpcandendx[6]->Fill(protoncandidate_endx[allp]);
      fpcandendy[6]->Fill(protoncandidate_endy[allp]);
      fpcandendz[6]->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[6]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[6]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[6]->Fill(*Evis);
      fQ2cal_all[6]->Fill(*Q2cal);



      }
   if(*OOFVflag==0 &&(*truthtop==2 ||*truthtop==3 ||*truthtop==4)&& (*trackpcand_origin !=1 && *trackcand_origin !=1)){
      trklen_muon[7]->Fill(*ftrklenmuoncand);
      trklen_proton[7]->Fill(*ftrklenprotoncand);
      costheta_muon[7]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[7]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[7]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[7]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[7]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[7]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[7]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[7]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);


      trunmean_muon[7]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[7]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[7]->Fill(*fPhiLep);
      phi_proton[7]->Fill(*fPhiHad);
      ntrks_proton[7]->Fill(*fNRecoPTrks);
      trkmom_muon[7]->Fill(*trackmomcandidate_mcs);
      fPmuon[7]->Fill(*fPlep);
      fPproton[7]->Fill(*fPhad);
      fPresomuon[7]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[7]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[7]->Fill(muonP.Angle(protonP)); 
      //h_thetamups[7]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[7]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

      Nhitsmuon[7]->Fill(*Nhits_muoncand);
      Nhitsproton[7]->Fill(*Nhits_protoncand);

      fvertexx[7]->Fill(*fvtxx);
      fvertexy[7]->Fill(*fvtxy);
      fvertexz[7]->Fill(*fvtxz);;

      fmucandstartx[7]->Fill(*trackstartxcandidate);
      fmucandstarty[7]->Fill(*trackstartycandidate);
      fmucandstartz[7]->Fill(*trackstartzcandidate);
      fmucandendx[7]->Fill(*trackendxcandidate);
      fmucandendy[7]->Fill(*trackendycandidate);
      fmucandendz[7]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[7]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[7]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[7]->Fill(protoncandidate_startz[allp]);
      fpcandendx[7]->Fill(protoncandidate_endx[allp]);
      fpcandendy[7]->Fill(protoncandidate_endy[allp]);
      fpcandendz[7]->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[7]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[7]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[7]->Fill(*Evis);
      fQ2cal_all[7]->Fill(*Q2cal);


      }
   
   if(*OOFVflag==0 && ((*trackpcand_origin ==1 && *trackcand_origin !=1) ||(*trackpcand_origin !=1 && *trackcand_origin ==1))){
      trklen_muon[8]->Fill(*ftrklenmuoncand);
      trklen_proton[8]->Fill(*ftrklenprotoncand);
      costheta_muon[8]->Fill(TMath::Cos(*fThetaLep));
      costheta_proton[8]->Fill(TMath::Cos(*fThetaHad));
      costhetatrue_muon[8]->Fill(*trackcand_parCosTheta);
      costhetatrue_proton[8]->Fill(*trackpcand_parCosTheta);

      costhetareso_muon[8]->Fill((TMath::Cos(*fThetaLep)-*trackcand_parCosTheta)/ *trackcand_parCosTheta);
      costhetareso_proton[8]->Fill((TMath::Cos(*fThetaHad)-*trackpcand_parCosTheta)/ *trackpcand_parCosTheta);

      cosphireso_muon[8]->Fill((TMath::Cos(*fPhiLep)-*trackcand_parCosPhi)/ *trackcand_parCosPhi);
      cosphireso_proton[8]->Fill((TMath::Cos(*fPhiHad)-*trackpcand_parCosPhi)/ *trackpcand_parCosPhi);


      trunmean_muon[8]->Fill(*TrunMean_cand*198.0/198.0);
      trunmean_proton[8]->Fill(*TrunMean_pcand*198.0/198.0);
      phi_muon[8]->Fill(*fPhiLep);
      phi_proton[8]->Fill(*fPhiHad);
      ntrks_proton[8]->Fill(*fNRecoPTrks);
      trkmom_muon[8]->Fill(*trackmomcandidate_mcs);
      fPmuon[8]->Fill(*fPlep);
      fPproton[8]->Fill(*fPhad);
      fPresomuon[8]->Fill((*fPlep- sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass))/sqrt(*trackcand_parE * *trackcand_parE- *trackcand_parMass * *trackcand_parMass));
      fPresoproton[8]->Fill((*fPhad- sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass))/sqrt(*trackpcand_parE * *trackpcand_parE- *trackpcand_parMass * *trackpcand_parMass));
      h_thetamups[8]->Fill(muonP.Angle(protonP)); 
      //h_thetamups[8]->Fill(TMath::Abs(*fThetaLep-*fThetaHad));
      h_phimups[8]->Fill(TMath::Abs(*fPhiLep-*fPhiHad));

      Nhitsmuon[8]->Fill(*Nhits_muoncand);
      Nhitsproton[8]->Fill(*Nhits_protoncand);

      fvertexx[8]->Fill(*fvtxx);
      fvertexy[8]->Fill(*fvtxy);
      fvertexz[8]->Fill(*fvtxz);;

      fmucandstartx[8]->Fill(*trackstartxcandidate);
      fmucandstarty[8]->Fill(*trackstartycandidate);
      fmucandstartz[8]->Fill(*trackstartzcandidate);
      fmucandendx[8]->Fill(*trackendxcandidate);
      fmucandendy[8]->Fill(*trackendycandidate);
      fmucandendz[8]->Fill(*trackendzcandidate);

   for(int allp=0; allp< protoncandidate_id.GetSize(); ++allp){
      fpcandstartx[8]->Fill(protoncandidate_startx[allp]);
      fpcandstarty[8]->Fill(protoncandidate_starty[allp]);
      fpcandstartz[8]->Fill(protoncandidate_startz[allp]);
      fpcandendx[8]->Fill(protoncandidate_endx[allp]);
      fpcandendy[8]->Fill(protoncandidate_endy[allp]);
      fpcandendz[8]->Fill(protoncandidate_endz[allp]);
   }
   if(protoncandidate_id.GetSize()==2){ 
      h_thetapps[8]->Fill(TMath::Abs(protoncandidate_theta[0]-protoncandidate_theta[1]));
      h_phipps[8]->Fill(TMath::Abs(protoncandidate_phi[0]-protoncandidate_phi[1]));
   }
 
      /*fpcandstartx;
      fpcandstarty;
      fpcandstartz;
      fpcandendx;
      fpcandendy;
      fpcandendz;
      */
      fEvis_all[8]->Fill(*Evis);
      fQ2cal_all[8]->Fill(*Q2cal);


      }

      //if(*trackcand_parPDG !=13) cout<<"the PDG code of the muon candidate is "<<*trackcand_parPDG<<endl;
      //
      //if(*trackpcand_parPDG !=2212) cout<<"the PDG code of the proton pcandidate is "<<*trackpcand_parPDG<<endl;

      } //end of if theta mup> 2.5
    } //end of if the muon track contained in FV 


   return kTRUE;
}

void hanalysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void hanalysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   //-------------------------------------------------------------------------
   //gStyle->SetOptFit();
   //..Write out histograms to file
   fHistFile->cd();
   fHistFile->Write();  
   fHistFile->Close();
   cout << "Output file written" << endl;
   /*
   */
   //--------------------------------------------------------------------------
}
