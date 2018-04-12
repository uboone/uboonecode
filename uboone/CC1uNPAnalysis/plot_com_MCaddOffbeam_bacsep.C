

void plot_com_MCaddOffbeam_bacsep(){
  //loadStyle();
  int tune=3;
  TFile *input0 = new TFile("histAnalysis_onbeam_wCosmicCut.root");
  TFile *input1 = new TFile("histAnalysis_offbeam_wCosmicCut.root");
  TFile *input2;
  
  if (tune==3){input2= new TFile("histAnalysis_Tune3_wCosmicCut.root");}
  else{input2= new TFile("histAnalysis_MC_wCosmicCut.root");}

  gROOT->SetBatch();
/*  TFile *input0 = new TFile("histAnalysis_onbeam_muinFV.root");
  TFile *input1 = new TFile("histAnalysis_offbeam_muinFV.root");
  TFile *input2 = new TFile("histAnalysis_MC_muinFV.root");  
*/
/*
  TFile *input0 = new TFile("histAnalysis_onbeam_muOOFV.root");
  TFile *input1 = new TFile("histAnalysis_offbeam_muOOFV.root");
  TFile *input2 = new TFile("histAnalysis_MC_muOOFV.root"); 
*/

  //----------------------------------------------------------------------
    gStyle->SetOptStat(0000);
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    //gStyle->SetPadTickX(1);
    //gStyle->SetPadTickY(1);
    gStyle->SetPadColor(kWhite);
    gStyle->SetStatY(0.90);
    gStyle->SetStatX(0.90);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.2);
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelFont(62,"X");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleFont(62,"X");
    gStyle->SetTitleOffset(0.85,"X");
    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetLabelFont(62,"Y");
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleFont(62,"Y");
    gStyle->SetTitleOffset(1.0,"Y");
    gStyle->SetTitleX(0.22);
    gStyle->SetTitleY(0.98);
    gStyle->SetTitleSize(0.04,"t");
    //gStyle->SetTitleTextColor(kRed);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    //gStyle->SetTitleFontSize(0);
    //gStyle->SetCanvasColor(kWhite);
    //gStyle->SetPadGridX(kTRUE);
    //gStyle->SetPadGridY(kTRUE);
    //gStyle->SetGridStyle();
  //----------------------------------------------------------------------
  TH1D                  *h_ntrksp_allsel[3];

  h_ntrksp_allsel[0]=(TH1D*)input0->Get("ntrkspreco");
  h_ntrksp_allsel[0]->Sumw2();
 

  h_ntrksp_allsel[1]=(TH1D*)input1->Get("ntrkspreco");
  h_ntrksp_allsel[1]->Sumw2();

  h_ntrksp_allsel[2]=(TH1D*)input2->Get("ntrkspreco");
  h_ntrksp_allsel[2]->Sumw2();

  TH1D               *h_ntrksp_sig[2];
  TH1D               *h_ntrksp_bac[8];
  h_ntrksp_sig[0]=(TH1D*)input2->Get("ntrks_proton_0");
  h_ntrksp_sig[0]->Sumw2();

  h_ntrksp_bac[0]=(TH1D*)input2->Get("ntrks_proton_1");
  h_ntrksp_bac[0]->Sumw2();
 
  h_ntrksp_bac[1]=(TH1D*)input2->Get("ntrks_proton_2");
  h_ntrksp_bac[1]->Sumw2();

  h_ntrksp_bac[2]=(TH1D*)input2->Get("ntrks_proton_3");
  h_ntrksp_bac[2]->Sumw2();
 
  h_ntrksp_bac[3]=(TH1D*)input2->Get("ntrks_proton_4");
  h_ntrksp_bac[3]->Sumw2();

  h_ntrksp_bac[4]=(TH1D*)input2->Get("ntrks_proton_5");
  h_ntrksp_bac[4]->Sumw2();
 
  h_ntrksp_bac[5]=(TH1D*)input2->Get("ntrks_proton_6");
  h_ntrksp_bac[5]->Sumw2();

  h_ntrksp_bac[6]=(TH1D*)input2->Get("ntrks_proton_7");
  h_ntrksp_bac[6]->Sumw2();

  h_ntrksp_bac[7]=(TH1D*)input2->Get("ntrks_proton_8");
  h_ntrksp_bac[7]->Sumw2();
   //==========================================================
  TH1D                  *h_Evis_allsel[3];

  h_Evis_allsel[0]=(TH1D*)input0->Get("fEvis");
  h_Evis_allsel[0]->Rebin(4);
  h_Evis_allsel[0]->Sumw2();
 

  h_Evis_allsel[1]=(TH1D*)input1->Get("fEvis");
  h_Evis_allsel[1]->Rebin(4);
  h_Evis_allsel[1]->Sumw2();

  h_Evis_allsel[2]=(TH1D*)input2->Get("fEvis");
  h_Evis_allsel[2]->Rebin(4);
  h_Evis_allsel[2]->Sumw2();

  TH1D               *h_Evis_sig[2];
  TH1D               *h_Evis_bac[8];
  h_Evis_sig[0]=(TH1D*)input2->Get("fEvis_all_0");
  h_Evis_sig[0]->Rebin(4);
  h_Evis_sig[0]->Sumw2();

  h_Evis_bac[0]=(TH1D*)input2->Get("fEvis_all_1");
  h_Evis_bac[0]->Rebin(4);
  h_Evis_bac[0]->Sumw2();
 
  h_Evis_bac[1]=(TH1D*)input2->Get("fEvis_all_2");
  h_Evis_bac[1]->Rebin(4);
  h_Evis_bac[1]->Sumw2();

  h_Evis_bac[2]=(TH1D*)input2->Get("fEvis_all_3");
  h_Evis_bac[2]->Rebin(4);
  h_Evis_bac[2]->Sumw2();
 
  h_Evis_bac[3]=(TH1D*)input2->Get("fEvis_all_4");
  h_Evis_bac[3]->Rebin(4);
  h_Evis_bac[3]->Sumw2();

  h_Evis_bac[4]=(TH1D*)input2->Get("fEvis_all_5");
  h_Evis_bac[4]->Rebin(4);
  h_Evis_bac[4]->Sumw2();
 
  h_Evis_bac[5]=(TH1D*)input2->Get("fEvis_all_6");
  h_Evis_bac[5]->Rebin(4);
  h_Evis_bac[5]->Sumw2();

  h_Evis_bac[6]=(TH1D*)input2->Get("fEvis_all_7");
  h_Evis_bac[6]->Rebin(4);
  h_Evis_bac[6]->Sumw2();

  h_Evis_bac[7]=(TH1D*)input2->Get("fEvis_all_8");
  h_Evis_bac[7]->Rebin(4);
  h_Evis_bac[7]->Sumw2();
   //==============================================================
  TH1D                  *h_Q2cal_allsel[3];

  h_Q2cal_allsel[0]=(TH1D*)input0->Get("fQ2cal");
  h_Q2cal_allsel[0]->Rebin(4);
  h_Q2cal_allsel[0]->Sumw2();
 

  h_Q2cal_allsel[1]=(TH1D*)input1->Get("fQ2cal");
  h_Q2cal_allsel[1]->Rebin(4);
  h_Q2cal_allsel[1]->Sumw2();

  h_Q2cal_allsel[2]=(TH1D*)input2->Get("fQ2cal");
  h_Q2cal_allsel[2]->Rebin(4);
  h_Q2cal_allsel[2]->Sumw2();

  TH1D               *h_Q2cal_sig[2];
  TH1D               *h_Q2cal_bac[8];
  h_Q2cal_sig[0]=(TH1D*)input2->Get("fQ2cal_all_0");
  h_Q2cal_sig[0]->Rebin(4);
  h_Q2cal_sig[0]->Sumw2();

  h_Q2cal_bac[0]=(TH1D*)input2->Get("fQ2cal_all_1");
  h_Q2cal_bac[0]->Rebin(4);
  h_Q2cal_bac[0]->Sumw2();
 
  h_Q2cal_bac[1]=(TH1D*)input2->Get("fQ2cal_all_2");
  h_Q2cal_bac[1]->Rebin(4);
  h_Q2cal_bac[1]->Sumw2();

  h_Q2cal_bac[2]=(TH1D*)input2->Get("fQ2cal_all_3");
  h_Q2cal_bac[2]->Rebin(4);
  h_Q2cal_bac[2]->Sumw2();
 
  h_Q2cal_bac[3]=(TH1D*)input2->Get("fQ2cal_all_4");
  h_Q2cal_bac[3]->Rebin(4);
  h_Q2cal_bac[3]->Sumw2();

  h_Q2cal_bac[4]=(TH1D*)input2->Get("fQ2cal_all_5");
  h_Q2cal_bac[4]->Rebin(4);
  h_Q2cal_bac[4]->Sumw2();
 
  h_Q2cal_bac[5]=(TH1D*)input2->Get("fQ2cal_all_6");
  h_Q2cal_bac[5]->Rebin(4);
  h_Q2cal_bac[5]->Sumw2();

  h_Q2cal_bac[6]=(TH1D*)input2->Get("fQ2cal_all_7");
  h_Q2cal_bac[6]->Rebin(4);
  h_Q2cal_bac[6]->Sumw2();

  h_Q2cal_bac[7]=(TH1D*)input2->Get("fQ2cal_all_8");
  h_Q2cal_bac[7]->Rebin(4);
  h_Q2cal_bac[7]->Sumw2();






  //===============================================================
  TH1D                  *h_range_allsel[3];

  h_range_allsel[0]=(TH1D*)input0->Get("trklen_cand");
  h_range_allsel[0]->Rebin(4); 
  h_range_allsel[0]->Sumw2();
 

  h_range_allsel[1]=(TH1D*)input1->Get("trklen_cand");
  h_range_allsel[1]->Rebin(4);
  h_range_allsel[1]->Sumw2();

  h_range_allsel[2]=(TH1D*)input2->Get("trklen_cand");
  h_range_allsel[2]->Rebin(4);
  h_range_allsel[2]->Sumw2();

  TH1D               *h_range_sig[2];
  TH1D               *h_range_bac[8];
  h_range_sig[0]=(TH1D*)input2->Get("trklen_muon_0");
  h_range_sig[0]->Rebin(4);
  h_range_sig[0]->Sumw2();

  h_range_bac[0]=(TH1D*)input2->Get("trklen_muon_1");
  h_range_bac[0]->Rebin(4);
  h_range_bac[0]->Sumw2();
 
  h_range_bac[1]=(TH1D*)input2->Get("trklen_muon_2");
  h_range_bac[1]->Rebin(4);
  h_range_bac[1]->Sumw2();

  h_range_bac[2]=(TH1D*)input2->Get("trklen_muon_3");
  h_range_bac[2]->Rebin(4);
  h_range_bac[2]->Sumw2();
 
  h_range_bac[3]=(TH1D*)input2->Get("trklen_muon_4");
  h_range_bac[3]->Rebin(4);
  h_range_bac[3]->Sumw2();

  h_range_bac[4]=(TH1D*)input2->Get("trklen_muon_5");
  h_range_bac[4]->Rebin(4);
  h_range_bac[4]->Sumw2();
 
  h_range_bac[5]=(TH1D*)input2->Get("trklen_muon_6");
  h_range_bac[5]->Rebin(4);
  h_range_bac[5]->Sumw2();

  h_range_bac[6]=(TH1D*)input2->Get("trklen_muon_7");
  h_range_bac[6]->Rebin(4);
  h_range_bac[6]->Sumw2();

  h_range_bac[7]=(TH1D*)input2->Get("trklen_muon_8");
  h_range_bac[7]->Rebin(4);
  h_range_bac[7]->Sumw2();
   //===============================================================
  TH1D                  *h_trkmom_allsel[3];

  h_trkmom_allsel[0]=(TH1D*)input0->Get("trkmomlep");
  h_trkmom_allsel[0]->Rebin(4); 
  h_trkmom_allsel[0]->Sumw2();
 

  h_trkmom_allsel[1]=(TH1D*)input1->Get("trkmomlep");
  h_trkmom_allsel[1]->Rebin(4);
  h_trkmom_allsel[1]->Sumw2();

  h_trkmom_allsel[2]=(TH1D*)input2->Get("trkmomlep");
  h_trkmom_allsel[2]->Rebin(4);
  h_trkmom_allsel[2]->Sumw2();

  TH1D               *h_trkmom_sig[2];
  TH1D               *h_trkmom_bac[8];
  h_trkmom_sig[0]=(TH1D*)input2->Get("trkmom_muon_0");
  h_trkmom_sig[0]->Rebin(4);
  h_trkmom_sig[0]->Sumw2();

  h_trkmom_bac[0]=(TH1D*)input2->Get("trkmom_muon_1");
  h_trkmom_bac[0]->Rebin(4);
  h_trkmom_bac[0]->Sumw2();
 
  h_trkmom_bac[1]=(TH1D*)input2->Get("trkmom_muon_2");
  h_trkmom_bac[1]->Rebin(4);
  h_trkmom_bac[1]->Sumw2();

  h_trkmom_bac[2]=(TH1D*)input2->Get("trkmom_muon_3");
  h_trkmom_bac[2]->Rebin(4);
  h_trkmom_bac[2]->Sumw2();
 
  h_trkmom_bac[3]=(TH1D*)input2->Get("trkmom_muon_4");
  h_trkmom_bac[3]->Rebin(4);
  h_trkmom_bac[3]->Sumw2();

  h_trkmom_bac[4]=(TH1D*)input2->Get("trkmom_muon_5");
  h_trkmom_bac[4]->Rebin(4);
  h_trkmom_bac[4]->Sumw2();
 
  h_trkmom_bac[5]=(TH1D*)input2->Get("trkmom_muon_6");
  h_trkmom_bac[5]->Rebin(4);
  h_trkmom_bac[5]->Sumw2();

  h_trkmom_bac[6]=(TH1D*)input2->Get("trkmom_muon_7");
  h_trkmom_bac[6]->Rebin(4);
  h_trkmom_bac[6]->Sumw2();

  h_trkmom_bac[7]=(TH1D*)input2->Get("trkmom_muon_8");
  h_trkmom_bac[7]->Rebin(4);
  h_trkmom_bac[7]->Sumw2();



  //=============================================================== 
  TH1D                  *h_prange_allsel[3];
  TH1D                  *h_prange_allzoom[3];
  h_prange_allsel[0]=(TH1D*)input0->Get("trklen_pcand");
  h_prange_allzoom[0]=(TH1D*)h_prange_allsel[0]->Clone();
  h_prange_allzoom[0]->Rebin(2);
  h_prange_allzoom[0]->Sumw2();
  h_prange_allsel[0]->Rebin(5); 
  h_prange_allsel[0]->Sumw2();
 

  h_prange_allsel[1]=(TH1D*)input1->Get("trklen_pcand");
  h_prange_allzoom[1]=(TH1D*)h_prange_allsel[1]->Clone();
  h_prange_allzoom[1]->Rebin(2);
  h_prange_allzoom[1]->Sumw2();
  h_prange_allsel[1]->Rebin(5);
  h_prange_allsel[1]->Sumw2();
 

  h_tmp = (TH1D*)input2->Get("trklen_pcand");
//  h_prange_allsel[2]=(TH1D*)input2->Get("trklen_pcand");
  h_prange_allsel[2]=(TH1D*)h_tmp->Clone();
  h_prange_allzoom[2]=(TH1D*)h_prange_allsel[2]->Clone();
  h_prange_allzoom[2]->Rebin(2);
  h_prange_allzoom[2]->Sumw2();
  h_prange_allsel[2]->Rebin(5);
  h_prange_allsel[2]->Sumw2();

  TH1D               *h_prange_sig[2];
  TH1D               *h_prange_bac[8];
  TH1D               *h_prange_sigzoom[2];
  TH1D               *h_prange_baczoom[8];
  h_prange_sig[0]=(TH1D*)input2->Get("trklen_proton_0");
  h_prange_sigzoom[0]=(TH1D*)h_prange_sig[0]->Clone();
  h_prange_sigzoom[0]->Rebin(2);
  h_prange_sigzoom[0]->Sumw2();
  h_prange_sig[0]->Rebin(5);
  h_prange_sig[0]->Sumw2();

  h_prange_bac[0]=(TH1D*)input2->Get("trklen_proton_1");
  h_prange_baczoom[0]=(TH1D*)h_prange_bac[0]->Clone();
  h_prange_baczoom[0]->Sumw2();
  h_prange_bac[0]->Rebin(5);
  h_prange_bac[0]->Sumw2();
 
  h_prange_bac[1]=(TH1D*)input2->Get("trklen_proton_2");
  h_prange_baczoom[1]=(TH1D*)h_prange_bac[1]->Clone();
  h_prange_baczoom[1]->Rebin(2);
  h_prange_baczoom[1]->Sumw2();
  h_prange_bac[1]->Rebin(5);
  h_prange_bac[1]->Sumw2();

  h_prange_bac[2]=(TH1D*)input2->Get("trklen_proton_3");
  h_prange_baczoom[2]=(TH1D*)h_prange_bac[2]->Clone();
  h_prange_baczoom[2]->Rebin(2);
  h_prange_baczoom[2]->Sumw2();
  h_prange_bac[2]->Rebin(5);
  h_prange_bac[2]->Sumw2();
 
  h_prange_bac[3]=(TH1D*)input2->Get("trklen_proton_4");
  h_prange_baczoom[3]=(TH1D*)h_prange_bac[3]->Clone();
  h_prange_baczoom[3]->Rebin(2);
  h_prange_baczoom[3]->Sumw2();
  h_prange_bac[3]->Rebin(5);
  h_prange_bac[3]->Sumw2();

  h_prange_bac[4]=(TH1D*)input2->Get("trklen_proton_5");
  h_prange_baczoom[4]=(TH1D*)h_prange_bac[4]->Clone();
  h_prange_baczoom[4]->Rebin(2);
  h_prange_baczoom[4]->Sumw2();
  h_prange_bac[4]->Rebin(5);
  h_prange_bac[4]->Sumw2();
 
  h_prange_bac[5]=(TH1D*)input2->Get("trklen_proton_6");
  h_prange_baczoom[5]=(TH1D*)h_prange_bac[5]->Clone();
  h_prange_baczoom[5]->Rebin(2);
  h_prange_baczoom[5]->Sumw2();
  h_prange_bac[5]->Rebin(5);
  h_prange_bac[5]->Sumw2();

  h_prange_bac[6]=(TH1D*)input2->Get("trklen_proton_7");
  h_prange_baczoom[6]=(TH1D*)h_prange_bac[6]->Clone();
  h_prange_baczoom[6]->Rebin(2);
  h_prange_baczoom[6]->Sumw2();
  h_prange_bac[6]->Rebin(5);
  h_prange_bac[6]->Sumw2();

  h_prange_bac[7]=(TH1D*)input2->Get("trklen_proton_8");
  h_prange_baczoom[7]=(TH1D*)h_prange_bac[7]->Clone();
  h_prange_baczoom[7]->Rebin(2);
  h_prange_baczoom[7]->Sumw2();
  h_prange_bac[7]->Rebin(5);
  h_prange_bac[7]->Sumw2();
   //=============================================================
  TH1D                  *h_phi_allsel[3];

  h_phi_allsel[0]=(TH1D*)input0->Get("philep");
  h_phi_allsel[0]->Rebin(6); 
  h_phi_allsel[0]->Sumw2();
 

  h_phi_allsel[1]=(TH1D*)input1->Get("philep");
  h_phi_allsel[1]->Rebin(4);
  h_phi_allsel[1]->Sumw2();

  h_phi_allsel[2]=(TH1D*)input2->Get("philep");
  h_phi_allsel[2]->Rebin(4);
  h_phi_allsel[2]->Sumw2();

  TH1D               *h_phi_sig[2];
  TH1D               *h_phi_bac[8];
  h_phi_sig[0]=(TH1D*)input2->Get("phi_muon_0");
  h_phi_sig[0]->Rebin(4);
  h_phi_sig[0]->Sumw2();

  h_phi_bac[0]=(TH1D*)input2->Get("phi_muon_1");
  h_phi_bac[0]->Rebin(4);
  h_phi_bac[0]->Sumw2();
 
  h_phi_bac[1]=(TH1D*)input2->Get("phi_muon_2");
  h_phi_bac[1]->Rebin(4);
  h_phi_bac[1]->Sumw2();

  h_phi_bac[2]=(TH1D*)input2->Get("phi_muon_3");
  h_phi_bac[2]->Rebin(4);
  h_phi_bac[2]->Sumw2();
 
  h_phi_bac[3]=(TH1D*)input2->Get("phi_muon_4");
  h_phi_bac[3]->Rebin(4);
  h_phi_bac[3]->Sumw2();

  h_phi_bac[4]=(TH1D*)input2->Get("phi_muon_5");
  h_phi_bac[4]->Rebin(4);
  h_phi_bac[4]->Sumw2();
 
  h_phi_bac[5]=(TH1D*)input2->Get("phi_muon_6");
  h_phi_bac[5]->Rebin(4);
  h_phi_bac[5]->Sumw2();

  h_phi_bac[6]=(TH1D*)input2->Get("phi_muon_7");
  h_phi_bac[6]->Rebin(4);
  h_phi_bac[6]->Sumw2();

  h_phi_bac[7]=(TH1D*)input2->Get("phi_muon_8");
  h_phi_bac[7]->Rebin(4);
  h_phi_bac[7]->Sumw2();
   //==================================================================== 
  TH1D                  *h_pphi_allsel[3];

  h_pphi_allsel[0]=(TH1D*)input0->Get("phihad");
  h_pphi_allsel[0]->Rebin(4); 
  h_pphi_allsel[0]->Sumw2();
 

  h_pphi_allsel[1]=(TH1D*)input1->Get("phihad");
  h_pphi_allsel[1]->Rebin(4);
  h_pphi_allsel[1]->Sumw2();

  h_pphi_allsel[2]=(TH1D*)input2->Get("phihad");
  h_pphi_allsel[2]->Rebin(4);
  h_pphi_allsel[2]->Sumw2();

  TH1D               *h_pphi_sig[2];
  TH1D               *h_pphi_bac[8];
  h_pphi_sig[0]=(TH1D*)input2->Get("phi_proton_0");
  h_pphi_sig[0]->Rebin(4);
  h_pphi_sig[0]->Sumw2();

  h_pphi_bac[0]=(TH1D*)input2->Get("phi_proton_1");
  h_pphi_bac[0]->Rebin(4);
  h_pphi_bac[0]->Sumw2();
 
  h_pphi_bac[1]=(TH1D*)input2->Get("phi_proton_2");
  h_pphi_bac[1]->Rebin(4);
  h_pphi_bac[1]->Sumw2();

  h_pphi_bac[2]=(TH1D*)input2->Get("phi_proton_3");
  h_pphi_bac[2]->Rebin(4);
  h_pphi_bac[2]->Sumw2();
 
  h_pphi_bac[3]=(TH1D*)input2->Get("phi_proton_4");
  h_pphi_bac[3]->Rebin(4);
  h_pphi_bac[3]->Sumw2();

  h_pphi_bac[4]=(TH1D*)input2->Get("phi_proton_5");
  h_pphi_bac[4]->Rebin(4);
  h_pphi_bac[4]->Sumw2();
 
  h_pphi_bac[5]=(TH1D*)input2->Get("phi_proton_6");
  h_pphi_bac[5]->Rebin(4);
  h_pphi_bac[5]->Sumw2();

  h_pphi_bac[6]=(TH1D*)input2->Get("phi_proton_7");
  h_pphi_bac[6]->Rebin(4);
  h_pphi_bac[6]->Sumw2();
 
  h_pphi_bac[7]=(TH1D*)input2->Get("phi_proton_8");
  h_pphi_bac[7]->Rebin(4);
  h_pphi_bac[7]->Sumw2();
   //================================================================= 
  TH1D                  *h_costheta_allsel[3];

  h_costheta_allsel[0]=(TH1D*)input0->Get("costhetalep");
  h_costheta_allsel[0]->Rebin(5); 
  h_costheta_allsel[0]->Sumw2();
 

  h_costheta_allsel[1]=(TH1D*)input1->Get("costhetalep");
  h_costheta_allsel[1]->Rebin(5);
  h_costheta_allsel[1]->Sumw2();

  h_costheta_allsel[2]=(TH1D*)input2->Get("costhetalep");
  h_costheta_allsel[2]->Rebin(5);
  h_costheta_allsel[2]->Sumw2();

  TH1D               *h_costheta_sig[2];
  TH1D               *h_costheta_bac[8];
  h_costheta_sig[0]=(TH1D*)input2->Get("costheta_muon_0");
  h_costheta_sig[0]->Rebin(5);
  h_costheta_sig[0]->Sumw2();

  h_costheta_bac[0]=(TH1D*)input2->Get("costheta_muon_1");
  h_costheta_bac[0]->Rebin(5);
  h_costheta_bac[0]->Sumw2();
 
  h_costheta_bac[1]=(TH1D*)input2->Get("costheta_muon_2");
  h_costheta_bac[1]->Rebin(5);
  h_costheta_bac[1]->Sumw2();

  h_costheta_bac[2]=(TH1D*)input2->Get("costheta_muon_3");
  h_costheta_bac[2]->Rebin(5);
  h_costheta_bac[2]->Sumw2();
 
  h_costheta_bac[3]=(TH1D*)input2->Get("costheta_muon_4");
  h_costheta_bac[3]->Rebin(5);
  h_costheta_bac[3]->Sumw2();

  h_costheta_bac[4]=(TH1D*)input2->Get("costheta_muon_5");
  h_costheta_bac[4]->Rebin(5);
  h_costheta_bac[4]->Sumw2();
 
  h_costheta_bac[5]=(TH1D*)input2->Get("costheta_muon_6");
  h_costheta_bac[5]->Rebin(5);
  h_costheta_bac[5]->Sumw2();

  h_costheta_bac[6]=(TH1D*)input2->Get("costheta_muon_7");
  h_costheta_bac[6]->Rebin(5);
  h_costheta_bac[6]->Sumw2();

  h_costheta_bac[7]=(TH1D*)input2->Get("costheta_muon_8");
  h_costheta_bac[7]->Rebin(5);
  h_costheta_bac[7]->Sumw2();
   //========================================================
  TH1D                  *h_pcostheta_allsel[3];

  h_pcostheta_allsel[0]=(TH1D*)input0->Get("costhetahad");
  h_pcostheta_allsel[0]->Rebin(5); 
  h_pcostheta_allsel[0]->Sumw2();
 

  h_pcostheta_allsel[1]=(TH1D*)input1->Get("costhetahad");
  h_pcostheta_allsel[1]->Rebin(5);
  h_pcostheta_allsel[1]->Sumw2();

  h_pcostheta_allsel[2]=(TH1D*)input2->Get("costhetahad");
  h_pcostheta_allsel[2]->Rebin(5);
  h_pcostheta_allsel[2]->Sumw2();

  TH1D               *h_pcostheta_sig[2];
  TH1D               *h_pcostheta_bac[8];
  h_pcostheta_sig[0]=(TH1D*)input2->Get("costheta_proton_0");
  h_pcostheta_sig[0]->Rebin(5);
  h_pcostheta_sig[0]->Sumw2();

  h_pcostheta_bac[0]=(TH1D*)input2->Get("costheta_proton_1");
  h_pcostheta_bac[0]->Rebin(5);
  h_pcostheta_bac[0]->Sumw2();
 
  h_pcostheta_bac[1]=(TH1D*)input2->Get("costheta_proton_2");
  h_pcostheta_bac[1]->Rebin(5);
  h_pcostheta_bac[1]->Sumw2();

  h_pcostheta_bac[2]=(TH1D*)input2->Get("costheta_proton_3");
  h_pcostheta_bac[2]->Rebin(5);
  h_pcostheta_bac[2]->Sumw2();
 
  h_pcostheta_bac[3]=(TH1D*)input2->Get("costheta_proton_4");
  h_pcostheta_bac[3]->Rebin(5);
  h_pcostheta_bac[3]->Sumw2();

  h_pcostheta_bac[4]=(TH1D*)input2->Get("costheta_proton_5");
  h_pcostheta_bac[4]->Rebin(5);
  h_pcostheta_bac[4]->Sumw2();
 
  h_pcostheta_bac[5]=(TH1D*)input2->Get("costheta_proton_6");
  h_pcostheta_bac[5]->Rebin(5);
  h_pcostheta_bac[5]->Sumw2();

  h_pcostheta_bac[6]=(TH1D*)input2->Get("costheta_proton_7");
  h_pcostheta_bac[6]->Rebin(5);
  h_pcostheta_bac[6]->Sumw2();

  h_pcostheta_bac[7]=(TH1D*)input2->Get("costheta_proton_8");
  h_pcostheta_bac[7]->Rebin(5);
  h_pcostheta_bac[7]->Sumw2();

  //=================================================================
  TH1D                  *h_trunmean_muon_allsel[3];

  h_trunmean_muon_allsel[0]=(TH1D*)input0->Get("trunmean_cand");
  h_trunmean_muon_allsel[0]->Rebin(5); 
  h_trunmean_muon_allsel[0]->Sumw2();
 

  h_trunmean_muon_allsel[1]=(TH1D*)input1->Get("trunmean_cand");
  h_trunmean_muon_allsel[1]->Rebin(5);
  h_trunmean_muon_allsel[1]->Sumw2();

  h_trunmean_muon_allsel[2]=(TH1D*)input2->Get("trunmean_cand");
  h_trunmean_muon_allsel[2]->Rebin(5);
  h_trunmean_muon_allsel[2]->Sumw2();

  TH1D               *h_trunmean_muon_sig[2];
  TH1D               *h_trunmean_muon_bac[8];
  h_trunmean_muon_sig[0]=(TH1D*)input2->Get("trunmean_muon_0");
  h_trunmean_muon_sig[0]->Rebin(5);
  h_trunmean_muon_sig[0]->Sumw2();

  h_trunmean_muon_bac[0]=(TH1D*)input2->Get("trunmean_muon_1");
  h_trunmean_muon_bac[0]->Rebin(5);
  h_trunmean_muon_bac[0]->Sumw2();
 
  h_trunmean_muon_bac[1]=(TH1D*)input2->Get("trunmean_muon_2");
  h_trunmean_muon_bac[1]->Rebin(5);
  h_trunmean_muon_bac[1]->Sumw2();

  h_trunmean_muon_bac[2]=(TH1D*)input2->Get("trunmean_muon_3");
  h_trunmean_muon_bac[2]->Rebin(5);
  h_trunmean_muon_bac[2]->Sumw2();
 
  h_trunmean_muon_bac[3]=(TH1D*)input2->Get("trunmean_muon_4");
  h_trunmean_muon_bac[3]->Rebin(5);
  h_trunmean_muon_bac[3]->Sumw2();

  h_trunmean_muon_bac[4]=(TH1D*)input2->Get("trunmean_muon_5");
  h_trunmean_muon_bac[4]->Rebin(5);
  h_trunmean_muon_bac[4]->Sumw2();
 
  h_trunmean_muon_bac[5]=(TH1D*)input2->Get("trunmean_muon_6");
  h_trunmean_muon_bac[5]->Rebin(5);
  h_trunmean_muon_bac[5]->Sumw2();

  h_trunmean_muon_bac[6]=(TH1D*)input2->Get("trunmean_muon_7");
  h_trunmean_muon_bac[6]->Rebin(5);
  h_trunmean_muon_bac[6]->Sumw2();

  h_trunmean_muon_bac[7]=(TH1D*)input2->Get("trunmean_muon_8");
  h_trunmean_muon_bac[7]->Rebin(5);
  h_trunmean_muon_bac[7]->Sumw2();
  //===================================================================
  TH1D                  *h_trunmean_proton_allsel[3];

  h_trunmean_proton_allsel[0]=(TH1D*)input0->Get("trunmean_pcand");
  h_trunmean_proton_allsel[0]->Rebin(5); 
  h_trunmean_proton_allsel[0]->Sumw2();
 

  h_trunmean_proton_allsel[1]=(TH1D*)input1->Get("trunmean_pcand");
  h_trunmean_proton_allsel[1]->Rebin(5);
  h_trunmean_proton_allsel[1]->Sumw2();

  h_trunmean_proton_allsel[2]=(TH1D*)input2->Get("trunmean_pcand");
  h_trunmean_proton_allsel[2]->Rebin(5);
  h_trunmean_proton_allsel[2]->Sumw2();

  TH1D               *h_trunmean_proton_sig[2];
  TH1D               *h_trunmean_proton_bac[8];
  h_trunmean_proton_sig[0]=(TH1D*)input2->Get("trunmean_proton_0");
  h_trunmean_proton_sig[0]->Rebin(5);
  h_trunmean_proton_sig[0]->Sumw2();

  h_trunmean_proton_bac[0]=(TH1D*)input2->Get("trunmean_proton_1");
  h_trunmean_proton_bac[0]->Rebin(5);
  h_trunmean_proton_bac[0]->Sumw2();
 
  h_trunmean_proton_bac[1]=(TH1D*)input2->Get("trunmean_proton_2");
  h_trunmean_proton_bac[1]->Rebin(5);
  h_trunmean_proton_bac[1]->Sumw2();

  h_trunmean_proton_bac[2]=(TH1D*)input2->Get("trunmean_proton_3");
  h_trunmean_proton_bac[2]->Rebin(5);
  h_trunmean_proton_bac[2]->Sumw2();
 
  h_trunmean_proton_bac[3]=(TH1D*)input2->Get("trunmean_proton_4");
  h_trunmean_proton_bac[3]->Rebin(5);
  h_trunmean_proton_bac[3]->Sumw2();

  h_trunmean_proton_bac[4]=(TH1D*)input2->Get("trunmean_proton_5");
  h_trunmean_proton_bac[4]->Rebin(5);
  h_trunmean_proton_bac[4]->Sumw2();
 
  h_trunmean_proton_bac[5]=(TH1D*)input2->Get("trunmean_proton_6");
  h_trunmean_proton_bac[5]->Rebin(5);
  h_trunmean_proton_bac[5]->Sumw2();

  h_trunmean_proton_bac[6]=(TH1D*)input2->Get("trunmean_proton_7");
  h_trunmean_proton_bac[6]->Rebin(5);
  h_trunmean_proton_bac[6]->Sumw2();

  h_trunmean_proton_bac[7]=(TH1D*)input2->Get("trunmean_proton_8");
  h_trunmean_proton_bac[7]->Rebin(5);
  h_trunmean_proton_bac[7]->Sumw2();
 



   //===============================================================
  TH1D                  *h_plep_allsel[3];

  h_plep_allsel[0]=(TH1D*)input0->Get("Plep");
  h_plep_allsel[0]->Rebin(5); 
  h_plep_allsel[0]->Sumw2();
 

  h_plep_allsel[1]=(TH1D*)input1->Get("Plep");
  h_plep_allsel[1]->Rebin(5);
  h_plep_allsel[1]->Sumw2();

  h_plep_allsel[2]=(TH1D*)input2->Get("Plep");
  h_plep_allsel[2]->Rebin(5);
  h_plep_allsel[2]->Sumw2();

  TH1D               *h_plep_sig[2];
  TH1D               *h_plep_bac[8];
  h_plep_sig[0]=(TH1D*)input2->Get("fPmuon_0");
  h_plep_sig[0]->Rebin(5);
  h_plep_sig[0]->Sumw2();

  h_plep_bac[0]=(TH1D*)input2->Get("fPmuon_1");
  h_plep_bac[0]->Rebin(5);
  h_plep_bac[0]->Sumw2();
 
  h_plep_bac[1]=(TH1D*)input2->Get("fPmuon_2");
  h_plep_bac[1]->Rebin(5);
  h_plep_bac[1]->Sumw2();

  h_plep_bac[2]=(TH1D*)input2->Get("fPmuon_3");
  h_plep_bac[2]->Rebin(5);
  h_plep_bac[2]->Sumw2();
 
  h_plep_bac[3]=(TH1D*)input2->Get("fPmuon_4");
  h_plep_bac[3]->Rebin(5);
  h_plep_bac[3]->Sumw2();

  h_plep_bac[4]=(TH1D*)input2->Get("fPmuon_5");
  h_plep_bac[4]->Rebin(5);
  h_plep_bac[4]->Sumw2();
 
  h_plep_bac[5]=(TH1D*)input2->Get("fPmuon_6");
  h_plep_bac[5]->Rebin(5);
  h_plep_bac[5]->Sumw2();

  h_plep_bac[6]=(TH1D*)input2->Get("fPmuon_7");
  h_plep_bac[6]->Rebin(5);
  h_plep_bac[6]->Sumw2();

  h_plep_bac[7]=(TH1D*)input2->Get("fPmuon_8");
  h_plep_bac[7]->Rebin(5);
  h_plep_bac[7]->Sumw2();
   //=============================================================== 
  TH1D                  *h_phad_allsel[3];

  h_phad_allsel[0]=(TH1D*)input0->Get("Phad");
  h_phad_allsel[0]->Rebin(5); 
  h_phad_allsel[0]->Sumw2();
 

  h_phad_allsel[1]=(TH1D*)input1->Get("Phad");
  h_phad_allsel[1]->Rebin(5);
  h_phad_allsel[1]->Sumw2();

  h_phad_allsel[2]=(TH1D*)input2->Get("Phad");
  h_phad_allsel[2]->Rebin(5);
  h_phad_allsel[2]->Sumw2();

  TH1D               *h_phad_sig[2];
  TH1D               *h_phad_bac[8];
  h_phad_sig[0]=(TH1D*)input2->Get("fPproton_0");
  h_phad_sig[0]->Rebin(5);
  h_phad_sig[0]->Sumw2();

  h_phad_bac[0]=(TH1D*)input2->Get("fPproton_1");
  h_phad_bac[0]->Rebin(5);
  h_phad_bac[0]->Sumw2();
 
  h_phad_bac[1]=(TH1D*)input2->Get("fPproton_2");
  h_phad_bac[1]->Rebin(5);
  h_phad_bac[1]->Sumw2();

  h_phad_bac[2]=(TH1D*)input2->Get("fPproton_3");
  h_phad_bac[2]->Rebin(5);
  h_phad_bac[2]->Sumw2();
 
  h_phad_bac[3]=(TH1D*)input2->Get("fPproton_4");
  h_phad_bac[3]->Rebin(5);
  h_phad_bac[3]->Sumw2();

  h_phad_bac[4]=(TH1D*)input2->Get("fPproton_5");
  h_phad_bac[4]->Rebin(5);
  h_phad_bac[4]->Sumw2();
 
  h_phad_bac[5]=(TH1D*)input2->Get("fPproton_6");
  h_phad_bac[5]->Rebin(5);
  h_phad_bac[5]->Sumw2();

  h_phad_bac[6]=(TH1D*)input2->Get("fPproton_7");
  h_phad_bac[6]->Rebin(5);
  h_phad_bac[6]->Sumw2();

  h_phad_bac[7]=(TH1D*)input2->Get("fPproton_8");
  h_phad_bac[7]->Rebin(5);
  h_phad_bac[7]->Sumw2();
   //=============================================================== 
  TH1D                  *h_thetamup_allsel[3];

  h_thetamup_allsel[0]=(TH1D*)input0->Get("h_thetamup");
  h_thetamup_allsel[0]->Rebin(4); 
  h_thetamup_allsel[0]->Sumw2();
 

  h_thetamup_allsel[1]=(TH1D*)input1->Get("h_thetamup");
  h_thetamup_allsel[1]->Rebin(4);
  h_thetamup_allsel[1]->Sumw2();

  h_thetamup_allsel[2]=(TH1D*)input2->Get("h_thetamup");
  h_thetamup_allsel[2]->Rebin(4);
  h_thetamup_allsel[2]->Sumw2();

  TH1D               *h_thetamup_sig[2];
  TH1D               *h_thetamup_bac[8];
  h_thetamup_sig[0]=(TH1D*)input2->Get("h_thetamup_0");
  h_thetamup_sig[0]->Rebin(4);
  h_thetamup_sig[0]->Sumw2();

  h_thetamup_bac[0]=(TH1D*)input2->Get("h_thetamup_1");
  h_thetamup_bac[0]->Rebin(4);
  h_thetamup_bac[0]->Sumw2();
 
  h_thetamup_bac[1]=(TH1D*)input2->Get("h_thetamup_2");
  h_thetamup_bac[1]->Rebin(4);
  h_thetamup_bac[1]->Sumw2();

  h_thetamup_bac[2]=(TH1D*)input2->Get("h_thetamup_3");
  h_thetamup_bac[2]->Rebin(4);
  h_thetamup_bac[2]->Sumw2();
 
  h_thetamup_bac[3]=(TH1D*)input2->Get("h_thetamup_4");
  h_thetamup_bac[3]->Rebin(4);
  h_thetamup_bac[3]->Sumw2();

  h_thetamup_bac[4]=(TH1D*)input2->Get("h_thetamup_5");
  h_thetamup_bac[4]->Rebin(4);
  h_thetamup_bac[4]->Sumw2();
 
  h_thetamup_bac[5]=(TH1D*)input2->Get("h_thetamup_6");
  h_thetamup_bac[5]->Rebin(4);
  h_thetamup_bac[5]->Sumw2();

  h_thetamup_bac[6]=(TH1D*)input2->Get("h_thetamup_7");
  h_thetamup_bac[6]->Rebin(4);
  h_thetamup_bac[6]->Sumw2();

  h_thetamup_bac[7]=(TH1D*)input2->Get("h_thetamup_8");
  h_thetamup_bac[7]->Rebin(4);
  h_thetamup_bac[7]->Sumw2();
   //======================================================= 
  TH1D                  *h_thetapp_allsel[3];

  h_thetapp_allsel[0]=(TH1D*)input0->Get("h_thetapp");
  h_thetapp_allsel[0]->Rebin(4); 
  h_thetapp_allsel[0]->Sumw2();
 

  h_thetapp_allsel[1]=(TH1D*)input1->Get("h_thetapp");
  h_thetapp_allsel[1]->Rebin(4);
  h_thetapp_allsel[1]->Sumw2();

  h_thetapp_allsel[2]=(TH1D*)input2->Get("h_thetapp");
  h_thetapp_allsel[2]->Rebin(4);
  h_thetapp_allsel[2]->Sumw2();

  TH1D               *h_thetapp_sig[2];
  TH1D               *h_thetapp_bac[8];
  h_thetapp_sig[0]=(TH1D*)input2->Get("h_thetapp_0");
  h_thetapp_sig[0]->Rebin(4);
  h_thetapp_sig[0]->Sumw2();

  h_thetapp_bac[0]=(TH1D*)input2->Get("h_thetapp_1");
  h_thetapp_bac[0]->Rebin(4);
  h_thetapp_bac[0]->Sumw2();
 
  h_thetapp_bac[1]=(TH1D*)input2->Get("h_thetapp_2");
  h_thetapp_bac[1]->Rebin(4);
  h_thetapp_bac[1]->Sumw2();

  h_thetapp_bac[2]=(TH1D*)input2->Get("h_thetapp_3");
  h_thetapp_bac[2]->Rebin(4);
  h_thetapp_bac[2]->Sumw2();
 
  h_thetapp_bac[3]=(TH1D*)input2->Get("h_thetapp_4");
  h_thetapp_bac[3]->Rebin(4);
  h_thetapp_bac[3]->Sumw2();

  h_thetapp_bac[4]=(TH1D*)input2->Get("h_thetapp_5");
  h_thetapp_bac[4]->Rebin(4);
  h_thetapp_bac[4]->Sumw2();
 
  h_thetapp_bac[5]=(TH1D*)input2->Get("h_thetapp_6");
  h_thetapp_bac[5]->Rebin(4);
  h_thetapp_bac[5]->Sumw2();

  h_thetapp_bac[6]=(TH1D*)input2->Get("h_thetapp_7");
  h_thetapp_bac[6]->Rebin(4);
  h_thetapp_bac[6]->Sumw2();
   
  h_thetapp_bac[7]=(TH1D*)input2->Get("h_thetapp_8");
  h_thetapp_bac[7]->Rebin(4);
  h_thetapp_bac[7]->Sumw2();
 
  //====================================================
  TH1D                  *h_phimup_allsel[3];

  h_phimup_allsel[0]=(TH1D*)input0->Get("h_phimup");
  h_phimup_allsel[0]->Rebin(4); 
  h_phimup_allsel[0]->Sumw2();
 

  h_phimup_allsel[1]=(TH1D*)input1->Get("h_phimup");
  h_phimup_allsel[1]->Rebin(4);
  h_phimup_allsel[1]->Sumw2();

  h_phimup_allsel[2]=(TH1D*)input2->Get("h_phimup");
  h_phimup_allsel[2]->Rebin(4);
  h_phimup_allsel[2]->Sumw2();

  TH1D               *h_phimup_sig[2];
  TH1D               *h_phimup_bac[8];
  h_phimup_sig[0]=(TH1D*)input2->Get("h_phimup_0");
  h_phimup_sig[0]->Rebin(4);
  h_phimup_sig[0]->Sumw2();

  h_phimup_bac[0]=(TH1D*)input2->Get("h_phimup_1");
  h_phimup_bac[0]->Rebin(4);
  h_phimup_bac[0]->Sumw2();
 
  h_phimup_bac[1]=(TH1D*)input2->Get("h_phimup_2");
  h_phimup_bac[1]->Rebin(4);
  h_phimup_bac[1]->Sumw2();

  h_phimup_bac[2]=(TH1D*)input2->Get("h_phimup_3");
  h_phimup_bac[2]->Rebin(4);
  h_phimup_bac[2]->Sumw2();
 
  h_phimup_bac[3]=(TH1D*)input2->Get("h_phimup_4");
  h_phimup_bac[3]->Rebin(4);
  h_phimup_bac[3]->Sumw2();

  h_phimup_bac[4]=(TH1D*)input2->Get("h_phimup_5");
  h_phimup_bac[4]->Rebin(4);
  h_phimup_bac[4]->Sumw2();
 
  h_phimup_bac[5]=(TH1D*)input2->Get("h_phimup_6");
  h_phimup_bac[5]->Rebin(4);
  h_phimup_bac[5]->Sumw2();

  h_phimup_bac[6]=(TH1D*)input2->Get("h_phimup_7");
  h_phimup_bac[6]->Rebin(4);
  h_phimup_bac[6]->Sumw2();

  h_phimup_bac[7]=(TH1D*)input2->Get("h_phimup_8");
  h_phimup_bac[7]->Rebin(4);
  h_phimup_bac[7]->Sumw2();

  //=======================================================
  TH1D                  *h_phipp_allsel[3];

  h_phipp_allsel[0]=(TH1D*)input0->Get("h_phipp");
  h_phipp_allsel[0]->Rebin(4); 
  h_phipp_allsel[0]->Sumw2();
 

  h_phipp_allsel[1]=(TH1D*)input1->Get("h_phipp");
  h_phipp_allsel[1]->Rebin(4);
  h_phipp_allsel[1]->Sumw2();

  h_phipp_allsel[2]=(TH1D*)input2->Get("h_phipp");
  h_phipp_allsel[2]->Rebin(4);
  h_phipp_allsel[2]->Sumw2();

  TH1D               *h_phipp_sig[2];
  TH1D               *h_phipp_bac[8];
  h_phipp_sig[0]=(TH1D*)input2->Get("h_phipp_0");
  h_phipp_sig[0]->Rebin(4);
  h_phipp_sig[0]->Sumw2();

  h_phipp_bac[0]=(TH1D*)input2->Get("h_phipp_1");
  h_phipp_bac[0]->Rebin(4);
  h_phipp_bac[0]->Sumw2();
 
  h_phipp_bac[1]=(TH1D*)input2->Get("h_phipp_2");
  h_phipp_bac[1]->Rebin(4);
  h_phipp_bac[1]->Sumw2();

  h_phipp_bac[2]=(TH1D*)input2->Get("h_phipp_3");
  h_phipp_bac[2]->Rebin(4);
  h_phipp_bac[2]->Sumw2();
 
  h_phipp_bac[3]=(TH1D*)input2->Get("h_phipp_4");
  h_phipp_bac[3]->Rebin(4);
  h_phipp_bac[3]->Sumw2();

  h_phipp_bac[4]=(TH1D*)input2->Get("h_phipp_5");
  h_phipp_bac[4]->Rebin(4);
  h_phipp_bac[4]->Sumw2();
 
  h_phipp_bac[5]=(TH1D*)input2->Get("h_phipp_6");
  h_phipp_bac[5]->Rebin(4);
  h_phipp_bac[5]->Sumw2();

  h_phipp_bac[6]=(TH1D*)input2->Get("h_phipp_7");
  h_phipp_bac[6]->Rebin(4);
  h_phipp_bac[6]->Sumw2();
   
  h_phipp_bac[7]=(TH1D*)input2->Get("h_phipp_8");
  h_phipp_bac[7]->Rebin(4);
  h_phipp_bac[7]->Sumw2();
   //========================================================

  TH1D                  *h_Nhitsmuon_allsel[3];

  h_Nhitsmuon_allsel[0]=(TH1D*)input0->Get("Nhitslep");
  h_Nhitsmuon_allsel[0]->Rebin(4); 
  h_Nhitsmuon_allsel[0]->Sumw2();
 

  h_Nhitsmuon_allsel[1]=(TH1D*)input1->Get("Nhitslep");
  h_Nhitsmuon_allsel[1]->Rebin(4);
  h_Nhitsmuon_allsel[1]->Sumw2();

  h_Nhitsmuon_allsel[2]=(TH1D*)input2->Get("Nhitslep");
  h_Nhitsmuon_allsel[2]->Rebin(4);
  h_Nhitsmuon_allsel[2]->Sumw2();

  TH1D               *h_Nhitsmuon_sig[2];
  TH1D               *h_Nhitsmuon_bac[8];
  h_Nhitsmuon_sig[0]=(TH1D*)input2->Get("Nhitsmuon_0");
  h_Nhitsmuon_sig[0]->Rebin(4);
  h_Nhitsmuon_sig[0]->Sumw2();

  h_Nhitsmuon_bac[0]=(TH1D*)input2->Get("Nhitsmuon_1");
  h_Nhitsmuon_bac[0]->Rebin(4);
  h_Nhitsmuon_bac[0]->Sumw2();
 
  h_Nhitsmuon_bac[1]=(TH1D*)input2->Get("Nhitsmuon_2");
  h_Nhitsmuon_bac[1]->Rebin(4);
  h_Nhitsmuon_bac[1]->Sumw2();

  h_Nhitsmuon_bac[2]=(TH1D*)input2->Get("Nhitsmuon_3");
  h_Nhitsmuon_bac[2]->Rebin(4);
  h_Nhitsmuon_bac[2]->Sumw2();
 
  h_Nhitsmuon_bac[3]=(TH1D*)input2->Get("Nhitsmuon_4");
  h_Nhitsmuon_bac[3]->Rebin(4);
  h_Nhitsmuon_bac[3]->Sumw2();

  h_Nhitsmuon_bac[4]=(TH1D*)input2->Get("Nhitsmuon_5");
  h_Nhitsmuon_bac[4]->Rebin(4);
  h_Nhitsmuon_bac[4]->Sumw2();
 
  h_Nhitsmuon_bac[5]=(TH1D*)input2->Get("Nhitsmuon_6");
  h_Nhitsmuon_bac[5]->Rebin(4);
  h_Nhitsmuon_bac[5]->Sumw2();

  h_Nhitsmuon_bac[6]=(TH1D*)input2->Get("Nhitsmuon_7");
  h_Nhitsmuon_bac[6]->Rebin(4);
  h_Nhitsmuon_bac[6]->Sumw2();

  h_Nhitsmuon_bac[7]=(TH1D*)input2->Get("Nhitsmuon_8");
  h_Nhitsmuon_bac[7]->Rebin(4);
  h_Nhitsmuon_bac[7]->Sumw2();
   //========================================================
  TH1D                  *h_Nhitsproton_allsel[3];

  h_Nhitsproton_allsel[0]=(TH1D*)input0->Get("Nhitshad");
  h_Nhitsproton_allsel[0]->Rebin(4); 
  h_Nhitsproton_allsel[0]->Sumw2();
 

  h_Nhitsproton_allsel[1]=(TH1D*)input1->Get("Nhitshad");
  h_Nhitsproton_allsel[1]->Rebin(4);
  h_Nhitsproton_allsel[1]->Sumw2();

  h_Nhitsproton_allsel[2]=(TH1D*)input2->Get("Nhitshad");
  h_Nhitsproton_allsel[2]->Rebin(4);
  h_Nhitsproton_allsel[2]->Sumw2();

  TH1D               *h_Nhitsproton_sig[2];
  TH1D               *h_Nhitsproton_bac[8];
  h_Nhitsproton_sig[0]=(TH1D*)input2->Get("Nhitsproton_0");
  h_Nhitsproton_sig[0]->Rebin(4);
  h_Nhitsproton_sig[0]->Sumw2();

  h_Nhitsproton_bac[0]=(TH1D*)input2->Get("Nhitsproton_1");
  h_Nhitsproton_bac[0]->Rebin(4);
  h_Nhitsproton_bac[0]->Sumw2();
 
  h_Nhitsproton_bac[1]=(TH1D*)input2->Get("Nhitsproton_2");
  h_Nhitsproton_bac[1]->Rebin(4);
  h_Nhitsproton_bac[1]->Sumw2();

  h_Nhitsproton_bac[2]=(TH1D*)input2->Get("Nhitsproton_3");
  h_Nhitsproton_bac[2]->Rebin(4);
  h_Nhitsproton_bac[2]->Sumw2();
 
  h_Nhitsproton_bac[3]=(TH1D*)input2->Get("Nhitsproton_4");
  h_Nhitsproton_bac[3]->Rebin(4);
  h_Nhitsproton_bac[3]->Sumw2();

  h_Nhitsproton_bac[4]=(TH1D*)input2->Get("Nhitsproton_5");
  h_Nhitsproton_bac[4]->Rebin(4);
  h_Nhitsproton_bac[4]->Sumw2();
 
  h_Nhitsproton_bac[5]=(TH1D*)input2->Get("Nhitsproton_6");
  h_Nhitsproton_bac[5]->Rebin(4);
  h_Nhitsproton_bac[5]->Sumw2();

  h_Nhitsproton_bac[6]=(TH1D*)input2->Get("Nhitsproton_7");
  h_Nhitsproton_bac[6]->Rebin(4);
  h_Nhitsproton_bac[6]->Sumw2();

  h_Nhitsproton_bac[7]=(TH1D*)input2->Get("Nhitsproton_8");
  h_Nhitsproton_bac[7]->Rebin(4);
  h_Nhitsproton_bac[7]->Sumw2();
   //===============================================================
  TH1D                  *h_vtxx_allsel[3];

  h_vtxx_allsel[0]=(TH1D*)input0->Get("fvertex_x");
  h_vtxx_allsel[0]->Rebin(4); 
  h_vtxx_allsel[0]->Sumw2();
 

  h_vtxx_allsel[1]=(TH1D*)input1->Get("fvertex_x");
  h_vtxx_allsel[1]->Rebin(4);
  h_vtxx_allsel[1]->Sumw2();

  h_vtxx_allsel[2]=(TH1D*)input2->Get("fvertex_x");
  h_vtxx_allsel[2]->Rebin(4);
  h_vtxx_allsel[2]->Sumw2();

  TH1D               *h_vtxx_sig[2];
  TH1D               *h_vtxx_bac[8];
  h_vtxx_sig[0]=(TH1D*)input2->Get("fvertexx_0");
  h_vtxx_sig[0]->Rebin(4);
  h_vtxx_sig[0]->Sumw2();

  h_vtxx_bac[0]=(TH1D*)input2->Get("fvertexx_1");
  h_vtxx_bac[0]->Rebin(4);
  h_vtxx_bac[0]->Sumw2();
 
  h_vtxx_bac[1]=(TH1D*)input2->Get("fvertexx_2");
  h_vtxx_bac[1]->Rebin(4);
  h_vtxx_bac[1]->Sumw2();

  h_vtxx_bac[2]=(TH1D*)input2->Get("fvertexx_3");
  h_vtxx_bac[2]->Rebin(4);
  h_vtxx_bac[2]->Sumw2();
 
  h_vtxx_bac[3]=(TH1D*)input2->Get("fvertexx_4");
  h_vtxx_bac[3]->Rebin(4);
  h_vtxx_bac[3]->Sumw2();

  h_vtxx_bac[4]=(TH1D*)input2->Get("fvertexx_5");
  h_vtxx_bac[4]->Rebin(4);
  h_vtxx_bac[4]->Sumw2();
 
  h_vtxx_bac[5]=(TH1D*)input2->Get("fvertexx_6");
  h_vtxx_bac[5]->Rebin(4);
  h_vtxx_bac[5]->Sumw2();

  h_vtxx_bac[6]=(TH1D*)input2->Get("fvertexx_7");
  h_vtxx_bac[6]->Rebin(4);
  h_vtxx_bac[6]->Sumw2();

  h_vtxx_bac[7]=(TH1D*)input2->Get("fvertexx_8");
  h_vtxx_bac[7]->Rebin(4);
  h_vtxx_bac[7]->Sumw2();
   //===============================================================
  TH1D                  *h_vtxy_allsel[3];

  h_vtxy_allsel[0]=(TH1D*)input0->Get("fvertex_y");
  h_vtxy_allsel[0]->Rebin(4); 
  h_vtxy_allsel[0]->Sumw2();
 

  h_vtxy_allsel[1]=(TH1D*)input1->Get("fvertex_y");
  h_vtxy_allsel[1]->Rebin(4);
  h_vtxy_allsel[1]->Sumw2();

  h_vtxy_allsel[2]=(TH1D*)input2->Get("fvertex_y");
  h_vtxy_allsel[2]->Rebin(4);
  h_vtxy_allsel[2]->Sumw2();

  TH1D               *h_vtxy_sig[2];
  TH1D               *h_vtxy_bac[8];
  h_vtxy_sig[0]=(TH1D*)input2->Get("fvertexy_0");
  h_vtxy_sig[0]->Rebin(4);
  h_vtxy_sig[0]->Sumw2();

  h_vtxy_bac[0]=(TH1D*)input2->Get("fvertexy_1");
  h_vtxy_bac[0]->Rebin(4);
  h_vtxy_bac[0]->Sumw2();
 
  h_vtxy_bac[1]=(TH1D*)input2->Get("fvertexy_2");
  h_vtxy_bac[1]->Rebin(4);
  h_vtxy_bac[1]->Sumw2();

  h_vtxy_bac[2]=(TH1D*)input2->Get("fvertexy_3");
  h_vtxy_bac[2]->Rebin(4);
  h_vtxy_bac[2]->Sumw2();
 
  h_vtxy_bac[3]=(TH1D*)input2->Get("fvertexy_4");
  h_vtxy_bac[3]->Rebin(4);
  h_vtxy_bac[3]->Sumw2();

  h_vtxy_bac[4]=(TH1D*)input2->Get("fvertexy_5");
  h_vtxy_bac[4]->Rebin(4);
  h_vtxy_bac[4]->Sumw2();
 
  h_vtxy_bac[5]=(TH1D*)input2->Get("fvertexy_6");
  h_vtxy_bac[5]->Rebin(4);
  h_vtxy_bac[5]->Sumw2();

  h_vtxy_bac[6]=(TH1D*)input2->Get("fvertexy_7");
  h_vtxy_bac[6]->Rebin(4);
  h_vtxy_bac[6]->Sumw2();

  h_vtxy_bac[7]=(TH1D*)input2->Get("fvertexy_8");
  h_vtxy_bac[7]->Rebin(4);
  h_vtxy_bac[7]->Sumw2();
   //===============================================================
  TH1D                  *h_vtxz_allsel[3];

  h_vtxz_allsel[0]=(TH1D*)input0->Get("fvertex_z");
  h_vtxz_allsel[0]->Rebin(4); 
  h_vtxz_allsel[0]->Sumw2();
 

  h_vtxz_allsel[1]=(TH1D*)input1->Get("fvertex_z");
  h_vtxz_allsel[1]->Rebin(4);
  h_vtxz_allsel[1]->Sumw2();

  h_vtxz_allsel[2]=(TH1D*)input2->Get("fvertex_z");
  h_vtxz_allsel[2]->Rebin(4);
  h_vtxz_allsel[2]->Sumw2();

  TH1D               *h_vtxz_sig[2];
  TH1D               *h_vtxz_bac[8];
  h_vtxz_sig[0]=(TH1D*)input2->Get("fvertexz_0");
  h_vtxz_sig[0]->Rebin(4);
  h_vtxz_sig[0]->Sumw2();

  h_vtxz_bac[0]=(TH1D*)input2->Get("fvertexz_1");
  h_vtxz_bac[0]->Rebin(4);
  h_vtxz_bac[0]->Sumw2();
 
  h_vtxz_bac[1]=(TH1D*)input2->Get("fvertexz_2");
  h_vtxz_bac[1]->Rebin(4);
  h_vtxz_bac[1]->Sumw2();

  h_vtxz_bac[2]=(TH1D*)input2->Get("fvertexz_3");
  h_vtxz_bac[2]->Rebin(4);
  h_vtxz_bac[2]->Sumw2();
 
  h_vtxz_bac[3]=(TH1D*)input2->Get("fvertexz_4");
  h_vtxz_bac[3]->Rebin(4);
  h_vtxz_bac[3]->Sumw2();

  h_vtxz_bac[4]=(TH1D*)input2->Get("fvertexz_5");
  h_vtxz_bac[4]->Rebin(4);
  h_vtxz_bac[4]->Sumw2();
 
  h_vtxz_bac[5]=(TH1D*)input2->Get("fvertexz_6");
  h_vtxz_bac[5]->Rebin(4);
  h_vtxz_bac[5]->Sumw2();

  h_vtxz_bac[6]=(TH1D*)input2->Get("fvertexz_7");
  h_vtxz_bac[6]->Rebin(4);
  h_vtxz_bac[6]->Sumw2();
 
  h_vtxz_bac[7]=(TH1D*)input2->Get("fvertexz_8");
  h_vtxz_bac[7]->Rebin(4);
  h_vtxz_bac[7]->Sumw2();
 
  //====================================================================== 
  TH1D                  *h_mustartx_allsel[3];

  h_mustartx_allsel[0]=(TH1D*)input0->Get("fmucand_startx");
  h_mustartx_allsel[0]->Rebin(4); 
  h_mustartx_allsel[0]->Sumw2();
 

  h_mustartx_allsel[1]=(TH1D*)input1->Get("fmucand_startx");
  h_mustartx_allsel[1]->Rebin(4);
  h_mustartx_allsel[1]->Sumw2();

  h_mustartx_allsel[2]=(TH1D*)input2->Get("fmucand_startx");
  h_mustartx_allsel[2]->Rebin(4);
  h_mustartx_allsel[2]->Sumw2();

  TH1D               *h_mustartx_sig[2];
  TH1D               *h_mustartx_bac[8];
  h_mustartx_sig[0]=(TH1D*)input2->Get("fmucandstartx_0");
  h_mustartx_sig[0]->Rebin(4);
  h_mustartx_sig[0]->Sumw2();

  h_mustartx_bac[0]=(TH1D*)input2->Get("fmucandstartx_1");
  h_mustartx_bac[0]->Rebin(4);
  h_mustartx_bac[0]->Sumw2();
 
  h_mustartx_bac[1]=(TH1D*)input2->Get("fmucandstartx_2");
  h_mustartx_bac[1]->Rebin(4);
  h_mustartx_bac[1]->Sumw2();

  h_mustartx_bac[2]=(TH1D*)input2->Get("fmucandstartx_3");
  h_mustartx_bac[2]->Rebin(4);
  h_mustartx_bac[2]->Sumw2();
 
  h_mustartx_bac[3]=(TH1D*)input2->Get("fmucandstartx_4");
  h_mustartx_bac[3]->Rebin(4);
  h_mustartx_bac[3]->Sumw2();

  h_mustartx_bac[4]=(TH1D*)input2->Get("fmucandstartx_5");
  h_mustartx_bac[4]->Rebin(4);
  h_mustartx_bac[4]->Sumw2();
 
  h_mustartx_bac[5]=(TH1D*)input2->Get("fmucandstartx_6");
  h_mustartx_bac[5]->Rebin(4);
  h_mustartx_bac[5]->Sumw2();

  h_mustartx_bac[6]=(TH1D*)input2->Get("fmucandstartx_7");
  h_mustartx_bac[6]->Rebin(4);
  h_mustartx_bac[6]->Sumw2();

  h_mustartx_bac[7]=(TH1D*)input2->Get("fmucandstartx_8");
  h_mustartx_bac[7]->Rebin(4);
  h_mustartx_bac[7]->Sumw2();
   //=============================================================================== 
  TH1D                  *h_mustarty_allsel[3];

  h_mustarty_allsel[0]=(TH1D*)input0->Get("fmucand_starty");
  h_mustarty_allsel[0]->Rebin(4); 
  h_mustarty_allsel[0]->Sumw2();
 

  h_mustarty_allsel[1]=(TH1D*)input1->Get("fmucand_starty");
  h_mustarty_allsel[1]->Rebin(4);
  h_mustarty_allsel[1]->Sumw2();

  h_mustarty_allsel[2]=(TH1D*)input2->Get("fmucand_starty");
  h_mustarty_allsel[2]->Rebin(4);
  h_mustarty_allsel[2]->Sumw2();

  TH1D               *h_mustarty_sig[2];
  TH1D               *h_mustarty_bac[8];
  h_mustarty_sig[0]=(TH1D*)input2->Get("fmucandstarty_0");
  h_mustarty_sig[0]->Rebin(4);
  h_mustarty_sig[0]->Sumw2();

  h_mustarty_bac[0]=(TH1D*)input2->Get("fmucandstarty_1");
  h_mustarty_bac[0]->Rebin(4);
  h_mustarty_bac[0]->Sumw2();
 
  h_mustarty_bac[1]=(TH1D*)input2->Get("fmucandstarty_2");
  h_mustarty_bac[1]->Rebin(4);
  h_mustarty_bac[1]->Sumw2();

  h_mustarty_bac[2]=(TH1D*)input2->Get("fmucandstarty_3");
  h_mustarty_bac[2]->Rebin(4);
  h_mustarty_bac[2]->Sumw2();
 
  h_mustarty_bac[3]=(TH1D*)input2->Get("fmucandstarty_4");
  h_mustarty_bac[3]->Rebin(4);
  h_mustarty_bac[3]->Sumw2();

  h_mustarty_bac[4]=(TH1D*)input2->Get("fmucandstarty_5");
  h_mustarty_bac[4]->Rebin(4);
  h_mustarty_bac[4]->Sumw2();
 
  h_mustarty_bac[5]=(TH1D*)input2->Get("fmucandstarty_6");
  h_mustarty_bac[5]->Rebin(4);
  h_mustarty_bac[5]->Sumw2();

  h_mustarty_bac[6]=(TH1D*)input2->Get("fmucandstarty_7");
  h_mustarty_bac[6]->Rebin(4);
  h_mustarty_bac[6]->Sumw2();

  h_mustarty_bac[7]=(TH1D*)input2->Get("fmucandstarty_8");
  h_mustarty_bac[7]->Rebin(4);
  h_mustarty_bac[7]->Sumw2();
   //===================================================================
  TH1D                  *h_mustartz_allsel[3];

  h_mustartz_allsel[0]=(TH1D*)input0->Get("fmucand_startz");
  h_mustartz_allsel[0]->Rebin(4); 
  h_mustartz_allsel[0]->Sumw2();
 

  h_mustartz_allsel[1]=(TH1D*)input1->Get("fmucand_startz");
  h_mustartz_allsel[1]->Rebin(4);
  h_mustartz_allsel[1]->Sumw2();

  h_mustartz_allsel[2]=(TH1D*)input2->Get("fmucand_startz");
  h_mustartz_allsel[2]->Rebin(4);
  h_mustartz_allsel[2]->Sumw2();

  TH1D               *h_mustartz_sig[2];
  TH1D               *h_mustartz_bac[8];
  h_mustartz_sig[0]=(TH1D*)input2->Get("fmucandstartz_0");
  h_mustartz_sig[0]->Rebin(4);
  h_mustartz_sig[0]->Sumw2();

  h_mustartz_bac[0]=(TH1D*)input2->Get("fmucandstartz_1");
  h_mustartz_bac[0]->Rebin(4);
  h_mustartz_bac[0]->Sumw2();
 
  h_mustartz_bac[1]=(TH1D*)input2->Get("fmucandstartz_2");
  h_mustartz_bac[1]->Rebin(4);
  h_mustartz_bac[1]->Sumw2();

  h_mustartz_bac[2]=(TH1D*)input2->Get("fmucandstartz_3");
  h_mustartz_bac[2]->Rebin(4);
  h_mustartz_bac[2]->Sumw2();
 
  h_mustartz_bac[3]=(TH1D*)input2->Get("fmucandstartz_4");
  h_mustartz_bac[3]->Rebin(4);
  h_mustartz_bac[3]->Sumw2();

  h_mustartz_bac[4]=(TH1D*)input2->Get("fmucandstartz_5");
  h_mustartz_bac[4]->Rebin(4);
  h_mustartz_bac[4]->Sumw2();
 
  h_mustartz_bac[5]=(TH1D*)input2->Get("fmucandstartz_6");
  h_mustartz_bac[5]->Rebin(4);
  h_mustartz_bac[5]->Sumw2();

  h_mustartz_bac[6]=(TH1D*)input2->Get("fmucandstartz_7");
  h_mustartz_bac[6]->Rebin(4);
  h_mustartz_bac[6]->Sumw2();

  h_mustartz_bac[7]=(TH1D*)input2->Get("fmucandstartz_8");
  h_mustartz_bac[7]->Rebin(4);
  h_mustartz_bac[7]->Sumw2();


  //==========================================================================
  TH1D                  *h_muendx_allsel[3];

  h_muendx_allsel[0]=(TH1D*)input0->Get("fmucand_endx");
  h_muendx_allsel[0]->Rebin(4); 
  h_muendx_allsel[0]->Sumw2();
 

  h_muendx_allsel[1]=(TH1D*)input1->Get("fmucand_endx");
  h_muendx_allsel[1]->Rebin(4);
  h_muendx_allsel[1]->Sumw2();

  h_muendx_allsel[2]=(TH1D*)input2->Get("fmucand_endx");
  h_muendx_allsel[2]->Rebin(4);
  h_muendx_allsel[2]->Sumw2();

  TH1D               *h_muendx_sig[2];
  TH1D               *h_muendx_bac[8];
  h_muendx_sig[0]=(TH1D*)input2->Get("fmucandendx_0");
  h_muendx_sig[0]->Rebin(4);
  h_muendx_sig[0]->Sumw2();

  h_muendx_bac[0]=(TH1D*)input2->Get("fmucandendx_1");
  h_muendx_bac[0]->Rebin(4);
  h_muendx_bac[0]->Sumw2();
 
  h_muendx_bac[1]=(TH1D*)input2->Get("fmucandendx_2");
  h_muendx_bac[1]->Rebin(4);
  h_muendx_bac[1]->Sumw2();

  h_muendx_bac[2]=(TH1D*)input2->Get("fmucandendx_3");
  h_muendx_bac[2]->Rebin(4);
  h_muendx_bac[2]->Sumw2();
 
  h_muendx_bac[3]=(TH1D*)input2->Get("fmucandendx_4");
  h_muendx_bac[3]->Rebin(4);
  h_muendx_bac[3]->Sumw2();

  h_muendx_bac[4]=(TH1D*)input2->Get("fmucandendx_5");
  h_muendx_bac[4]->Rebin(4);
  h_muendx_bac[4]->Sumw2();
 
  h_muendx_bac[5]=(TH1D*)input2->Get("fmucandendx_6");
  h_muendx_bac[5]->Rebin(4);
  h_muendx_bac[5]->Sumw2();

  h_muendx_bac[6]=(TH1D*)input2->Get("fmucandendx_7");
  h_muendx_bac[6]->Rebin(4);
  h_muendx_bac[6]->Sumw2();

  h_muendx_bac[7]=(TH1D*)input2->Get("fmucandendx_8");
  h_muendx_bac[7]->Rebin(4);
  h_muendx_bac[7]->Sumw2();
   //=============================================================================== 
  TH1D                  *h_muendy_allsel[3];

  h_muendy_allsel[0]=(TH1D*)input0->Get("fmucand_endy");
  h_muendy_allsel[0]->Rebin(4); 
  h_muendy_allsel[0]->Sumw2();
 

  h_muendy_allsel[1]=(TH1D*)input1->Get("fmucand_endy");
  h_muendy_allsel[1]->Rebin(4);
  h_muendy_allsel[1]->Sumw2();

  h_muendy_allsel[2]=(TH1D*)input2->Get("fmucand_endy");
  h_muendy_allsel[2]->Rebin(4);
  h_muendy_allsel[2]->Sumw2();

  TH1D               *h_muendy_sig[2];
  TH1D               *h_muendy_bac[8];
  h_muendy_sig[0]=(TH1D*)input2->Get("fmucandendy_0");
  h_muendy_sig[0]->Rebin(4);
  h_muendy_sig[0]->Sumw2();

  h_muendy_bac[0]=(TH1D*)input2->Get("fmucandendy_1");
  h_muendy_bac[0]->Rebin(4);
  h_muendy_bac[0]->Sumw2();
 
  h_muendy_bac[1]=(TH1D*)input2->Get("fmucandendy_2");
  h_muendy_bac[1]->Rebin(4);
  h_muendy_bac[1]->Sumw2();

  h_muendy_bac[2]=(TH1D*)input2->Get("fmucandendy_3");
  h_muendy_bac[2]->Rebin(4);
  h_muendy_bac[2]->Sumw2();
 
  h_muendy_bac[3]=(TH1D*)input2->Get("fmucandendy_4");
  h_muendy_bac[3]->Rebin(4);
  h_muendy_bac[3]->Sumw2();

  h_muendy_bac[4]=(TH1D*)input2->Get("fmucandendy_5");
  h_muendy_bac[4]->Rebin(4);
  h_muendy_bac[4]->Sumw2();
 
  h_muendy_bac[5]=(TH1D*)input2->Get("fmucandendy_6");
  h_muendy_bac[5]->Rebin(4);
  h_muendy_bac[5]->Sumw2();

  h_muendy_bac[6]=(TH1D*)input2->Get("fmucandendy_7");
  h_muendy_bac[6]->Rebin(4);
  h_muendy_bac[6]->Sumw2();

  h_muendy_bac[7]=(TH1D*)input2->Get("fmucandendy_8");
  h_muendy_bac[7]->Rebin(4);
  h_muendy_bac[7]->Sumw2();
 

  //===================================================================
  TH1D                  *h_muendz_allsel[3];

  h_muendz_allsel[0]=(TH1D*)input0->Get("fmucand_endz");
  h_muendz_allsel[0]->Rebin(4); 
  h_muendz_allsel[0]->Sumw2();
 

  h_muendz_allsel[1]=(TH1D*)input1->Get("fmucand_endz");
  h_muendz_allsel[1]->Rebin(4);
  h_muendz_allsel[1]->Sumw2();

  h_muendz_allsel[2]=(TH1D*)input2->Get("fmucand_endz");
  h_muendz_allsel[2]->Rebin(4);
  h_muendz_allsel[2]->Sumw2();

  TH1D               *h_muendz_sig[2];
  TH1D               *h_muendz_bac[8];
  h_muendz_sig[0]=(TH1D*)input2->Get("fmucandendz_0");
  h_muendz_sig[0]->Rebin(4);
  h_muendz_sig[0]->Sumw2();

  h_muendz_bac[0]=(TH1D*)input2->Get("fmucandendz_1");
  h_muendz_bac[0]->Rebin(4);
  h_muendz_bac[0]->Sumw2();
 
  h_muendz_bac[1]=(TH1D*)input2->Get("fmucandendz_2");
  h_muendz_bac[1]->Rebin(4);
  h_muendz_bac[1]->Sumw2();

  h_muendz_bac[2]=(TH1D*)input2->Get("fmucandendz_3");
  h_muendz_bac[2]->Rebin(4);
  h_muendz_bac[2]->Sumw2();
 
  h_muendz_bac[3]=(TH1D*)input2->Get("fmucandendz_4");
  h_muendz_bac[3]->Rebin(4);
  h_muendz_bac[3]->Sumw2();

  h_muendz_bac[4]=(TH1D*)input2->Get("fmucandendz_5");
  h_muendz_bac[4]->Rebin(4);
  h_muendz_bac[4]->Sumw2();
 
  h_muendz_bac[5]=(TH1D*)input2->Get("fmucandendz_6");
  h_muendz_bac[5]->Rebin(4);
  h_muendz_bac[5]->Sumw2();

  h_muendz_bac[6]=(TH1D*)input2->Get("fmucandendz_7");
  h_muendz_bac[6]->Rebin(4);
  h_muendz_bac[6]->Sumw2();

  h_muendz_bac[7]=(TH1D*)input2->Get("fmucandendz_8");
  h_muendz_bac[7]->Rebin(4);
  h_muendz_bac[7]->Sumw2();
     
  //=========================================================================== 
  TH1D                  *h_pstartx_allsel[3];

  h_pstartx_allsel[0]=(TH1D*)input0->Get("fpcand_startx");
  h_pstartx_allsel[0]->Rebin(4); 
  h_pstartx_allsel[0]->Sumw2();
 

  h_pstartx_allsel[1]=(TH1D*)input1->Get("fpcand_startx");
  h_pstartx_allsel[1]->Rebin(4);
  h_pstartx_allsel[1]->Sumw2();

  h_pstartx_allsel[2]=(TH1D*)input2->Get("fpcand_startx");
  h_pstartx_allsel[2]->Rebin(4);
  h_pstartx_allsel[2]->Sumw2();

  TH1D               *h_pstartx_sig[2];
  TH1D               *h_pstartx_bac[8];
  h_pstartx_sig[0]=(TH1D*)input2->Get("fpcandstartx_0");
  h_pstartx_sig[0]->Rebin(4);
  h_pstartx_sig[0]->Sumw2();

  h_pstartx_bac[0]=(TH1D*)input2->Get("fpcandstartx_1");
  h_pstartx_bac[0]->Rebin(4);
  h_pstartx_bac[0]->Sumw2();
 
  h_pstartx_bac[1]=(TH1D*)input2->Get("fpcandstartx_2");
  h_pstartx_bac[1]->Rebin(4);
  h_pstartx_bac[1]->Sumw2();

  h_pstartx_bac[2]=(TH1D*)input2->Get("fpcandstartx_3");
  h_pstartx_bac[2]->Rebin(4);
  h_pstartx_bac[2]->Sumw2();
 
  h_pstartx_bac[3]=(TH1D*)input2->Get("fpcandstartx_4");
  h_pstartx_bac[3]->Rebin(4);
  h_pstartx_bac[3]->Sumw2();

  h_pstartx_bac[4]=(TH1D*)input2->Get("fpcandstartx_5");
  h_pstartx_bac[4]->Rebin(4);
  h_pstartx_bac[4]->Sumw2();
 
  h_pstartx_bac[5]=(TH1D*)input2->Get("fpcandstartx_6");
  h_pstartx_bac[5]->Rebin(4);
  h_pstartx_bac[5]->Sumw2();

  h_pstartx_bac[6]=(TH1D*)input2->Get("fpcandstartx_7");
  h_pstartx_bac[6]->Rebin(4);
  h_pstartx_bac[6]->Sumw2();

  h_pstartx_bac[7]=(TH1D*)input2->Get("fpcandstartx_8");
  h_pstartx_bac[7]->Rebin(4);
  h_pstartx_bac[7]->Sumw2();
   //=============================================================================== 

  TH1D                  *h_pstarty_allsel[3];

  h_pstarty_allsel[0]=(TH1D*)input0->Get("fpcand_starty");
  h_pstarty_allsel[0]->Rebin(4); 
  h_pstarty_allsel[0]->Sumw2();
 

  h_pstarty_allsel[1]=(TH1D*)input1->Get("fpcand_starty");
  h_pstarty_allsel[1]->Rebin(4);
  h_pstarty_allsel[1]->Sumw2();

  h_pstarty_allsel[2]=(TH1D*)input2->Get("fpcand_starty");
  h_pstarty_allsel[2]->Rebin(4);
  h_pstarty_allsel[2]->Sumw2();

  TH1D               *h_pstarty_sig[2];
  TH1D               *h_pstarty_bac[8];
  h_pstarty_sig[0]=(TH1D*)input2->Get("fpcandstarty_0");
  h_pstarty_sig[0]->Rebin(4);
  h_pstarty_sig[0]->Sumw2();

  h_pstarty_bac[0]=(TH1D*)input2->Get("fpcandstarty_1");
  h_pstarty_bac[0]->Rebin(4);
  h_pstarty_bac[0]->Sumw2();
 
  h_pstarty_bac[1]=(TH1D*)input2->Get("fpcandstarty_2");
  h_pstarty_bac[1]->Rebin(4);
  h_pstarty_bac[1]->Sumw2();

  h_pstarty_bac[2]=(TH1D*)input2->Get("fpcandstarty_3");
  h_pstarty_bac[2]->Rebin(4);
  h_pstarty_bac[2]->Sumw2();
 
  h_pstarty_bac[3]=(TH1D*)input2->Get("fpcandstarty_4");
  h_pstarty_bac[3]->Rebin(4);
  h_pstarty_bac[3]->Sumw2();

  h_pstarty_bac[4]=(TH1D*)input2->Get("fpcandstarty_5");
  h_pstarty_bac[4]->Rebin(4);
  h_pstarty_bac[4]->Sumw2();
 
  h_pstarty_bac[5]=(TH1D*)input2->Get("fpcandstarty_6");
  h_pstarty_bac[5]->Rebin(4);
  h_pstarty_bac[5]->Sumw2();

  h_pstarty_bac[6]=(TH1D*)input2->Get("fpcandstarty_7");
  h_pstarty_bac[6]->Rebin(4);
  h_pstarty_bac[6]->Sumw2();

  h_pstarty_bac[7]=(TH1D*)input2->Get("fpcandstarty_8");
  h_pstarty_bac[7]->Rebin(4);
  h_pstarty_bac[7]->Sumw2();
   //===================================================================
  TH1D                  *h_pstartz_allsel[3];

  h_pstartz_allsel[0]=(TH1D*)input0->Get("fpcand_startz");
  h_pstartz_allsel[0]->Rebin(4); 
  h_pstartz_allsel[0]->Sumw2();
 

  h_pstartz_allsel[1]=(TH1D*)input1->Get("fpcand_startz");
  h_pstartz_allsel[1]->Rebin(4);
  h_pstartz_allsel[1]->Sumw2();

  h_pstartz_allsel[2]=(TH1D*)input2->Get("fpcand_startz");
  h_pstartz_allsel[2]->Rebin(4);
  h_pstartz_allsel[2]->Sumw2();

  TH1D               *h_pstartz_sig[2];
  TH1D               *h_pstartz_bac[8];
  h_pstartz_sig[0]=(TH1D*)input2->Get("fpcandstartz_0");
  h_pstartz_sig[0]->Rebin(4);
  h_pstartz_sig[0]->Sumw2();

  h_pstartz_bac[0]=(TH1D*)input2->Get("fpcandstartz_1");
  h_pstartz_bac[0]->Rebin(4);
  h_pstartz_bac[0]->Sumw2();
 
  h_pstartz_bac[1]=(TH1D*)input2->Get("fpcandstartz_2");
  h_pstartz_bac[1]->Rebin(4);
  h_pstartz_bac[1]->Sumw2();

  h_pstartz_bac[2]=(TH1D*)input2->Get("fpcandstartz_3");
  h_pstartz_bac[2]->Rebin(4);
  h_pstartz_bac[2]->Sumw2();
 
  h_pstartz_bac[3]=(TH1D*)input2->Get("fpcandstartz_4");
  h_pstartz_bac[3]->Rebin(4);
  h_pstartz_bac[3]->Sumw2();

  h_pstartz_bac[4]=(TH1D*)input2->Get("fpcandstartz_5");
  h_pstartz_bac[4]->Rebin(4);
  h_pstartz_bac[4]->Sumw2();
 
  h_pstartz_bac[5]=(TH1D*)input2->Get("fpcandstartz_6");
  h_pstartz_bac[5]->Rebin(4);
  h_pstartz_bac[5]->Sumw2();

  h_pstartz_bac[6]=(TH1D*)input2->Get("fpcandstartz_7");
  h_pstartz_bac[6]->Rebin(4);
  h_pstartz_bac[6]->Sumw2();

  h_pstartz_bac[7]=(TH1D*)input2->Get("fpcandstartz_8");
  h_pstartz_bac[7]->Rebin(4);
  h_pstartz_bac[7]->Sumw2();

  //==========================================================================
  TH1D                  *h_pendx_allsel[3];

  h_pendx_allsel[0]=(TH1D*)input0->Get("fpcand_endx");
  h_pendx_allsel[0]->Rebin(4); 
  h_pendx_allsel[0]->Sumw2();
 

  h_pendx_allsel[1]=(TH1D*)input1->Get("fpcand_endx");
  h_pendx_allsel[1]->Rebin(4);
  h_pendx_allsel[1]->Sumw2();

  h_pendx_allsel[2]=(TH1D*)input2->Get("fpcand_endx");
  h_pendx_allsel[2]->Rebin(4);
  h_pendx_allsel[2]->Sumw2();

  TH1D               *h_pendx_sig[2];
  TH1D               *h_pendx_bac[8];
  h_pendx_sig[0]=(TH1D*)input2->Get("fpcandendx_0");
  h_pendx_sig[0]->Rebin(4);
  h_pendx_sig[0]->Sumw2();

  h_pendx_bac[0]=(TH1D*)input2->Get("fpcandendx_1");
  h_pendx_bac[0]->Rebin(4);
  h_pendx_bac[0]->Sumw2();
 
  h_pendx_bac[1]=(TH1D*)input2->Get("fpcandendx_2");
  h_pendx_bac[1]->Rebin(4);
  h_pendx_bac[1]->Sumw2();

  h_pendx_bac[2]=(TH1D*)input2->Get("fpcandendx_3");
  h_pendx_bac[2]->Rebin(4);
  h_pendx_bac[2]->Sumw2();
 
  h_pendx_bac[3]=(TH1D*)input2->Get("fpcandendx_4");
  h_pendx_bac[3]->Rebin(4);
  h_pendx_bac[3]->Sumw2();

  h_pendx_bac[4]=(TH1D*)input2->Get("fpcandendx_5");
  h_pendx_bac[4]->Rebin(4);
  h_pendx_bac[4]->Sumw2();
 
  h_pendx_bac[5]=(TH1D*)input2->Get("fpcandendx_6");
  h_pendx_bac[5]->Rebin(4);
  h_pendx_bac[5]->Sumw2();

  h_pendx_bac[6]=(TH1D*)input2->Get("fpcandendx_7");
  h_pendx_bac[6]->Rebin(4);
  h_pendx_bac[6]->Sumw2();

  h_pendx_bac[7]=(TH1D*)input2->Get("fpcandendx_8");
  h_pendx_bac[7]->Rebin(4);
  h_pendx_bac[7]->Sumw2();
   //=============================================================================== 
  TH1D                  *h_pendy_allsel[3];

  h_pendy_allsel[0]=(TH1D*)input0->Get("fpcand_endy");
  h_pendy_allsel[0]->Rebin(4); 
  h_pendy_allsel[0]->Sumw2();
 

  h_pendy_allsel[1]=(TH1D*)input1->Get("fpcand_endy");
  h_pendy_allsel[1]->Rebin(4);
  h_pendy_allsel[1]->Sumw2();

  h_pendy_allsel[2]=(TH1D*)input2->Get("fpcand_endy");
  h_pendy_allsel[2]->Rebin(4);
  h_pendy_allsel[2]->Sumw2();

  TH1D               *h_pendy_sig[2];
  TH1D               *h_pendy_bac[8];
  h_pendy_sig[0]=(TH1D*)input2->Get("fpcandendy_0");
  h_pendy_sig[0]->Rebin(4);
  h_pendy_sig[0]->Sumw2();

  h_pendy_bac[0]=(TH1D*)input2->Get("fpcandendy_1");
  h_pendy_bac[0]->Rebin(4);
  h_pendy_bac[0]->Sumw2();
 
  h_pendy_bac[1]=(TH1D*)input2->Get("fpcandendy_2");
  h_pendy_bac[1]->Rebin(4);
  h_pendy_bac[1]->Sumw2();

  h_pendy_bac[2]=(TH1D*)input2->Get("fpcandendy_3");
  h_pendy_bac[2]->Rebin(4);
  h_pendy_bac[2]->Sumw2();
 
  h_pendy_bac[3]=(TH1D*)input2->Get("fpcandendy_4");
  h_pendy_bac[3]->Rebin(4);
  h_pendy_bac[3]->Sumw2();

  h_pendy_bac[4]=(TH1D*)input2->Get("fpcandendy_5");
  h_pendy_bac[4]->Rebin(4);
  h_pendy_bac[4]->Sumw2();
 
  h_pendy_bac[5]=(TH1D*)input2->Get("fpcandendy_6");
  h_pendy_bac[5]->Rebin(4);
  h_pendy_bac[5]->Sumw2();

  h_pendy_bac[6]=(TH1D*)input2->Get("fpcandendy_7");
  h_pendy_bac[6]->Rebin(4);
  h_pendy_bac[6]->Sumw2();

  h_pendy_bac[7]=(TH1D*)input2->Get("fpcandendy_8");
  h_pendy_bac[7]->Rebin(4);
  h_pendy_bac[7]->Sumw2();
   //===================================================================
  TH1D                  *h_pendz_allsel[3];

  h_pendz_allsel[0]=(TH1D*)input0->Get("fpcand_endz");
  h_pendz_allsel[0]->Rebin(4); 
  h_pendz_allsel[0]->Sumw2();
 

  h_pendz_allsel[1]=(TH1D*)input1->Get("fpcand_endz");
  h_pendz_allsel[1]->Rebin(4);
  h_pendz_allsel[1]->Sumw2();

  h_pendz_allsel[2]=(TH1D*)input2->Get("fpcand_endz");
  h_pendz_allsel[2]->Rebin(4);
  h_pendz_allsel[2]->Sumw2();

  TH1D               *h_pendz_sig[2];
  TH1D               *h_pendz_bac[8];
  h_pendz_sig[0]=(TH1D*)input2->Get("fpcandendz_0");
  h_pendz_sig[0]->Rebin(4);
  h_pendz_sig[0]->Sumw2();

  h_pendz_bac[0]=(TH1D*)input2->Get("fpcandendz_1");
  h_pendz_bac[0]->Rebin(4);
  h_pendz_bac[0]->Sumw2();
 
  h_pendz_bac[1]=(TH1D*)input2->Get("fpcandendz_2");
  h_pendz_bac[1]->Rebin(4);
  h_pendz_bac[1]->Sumw2();

  h_pendz_bac[2]=(TH1D*)input2->Get("fpcandendz_3");
  h_pendz_bac[2]->Rebin(4);
  h_pendz_bac[2]->Sumw2();
 
  h_pendz_bac[3]=(TH1D*)input2->Get("fpcandendz_4");
  h_pendz_bac[3]->Rebin(4);
  h_pendz_bac[3]->Sumw2();

  h_pendz_bac[4]=(TH1D*)input2->Get("fpcandendz_5");
  h_pendz_bac[4]->Rebin(4);
  h_pendz_bac[4]->Sumw2();
 
  h_pendz_bac[5]=(TH1D*)input2->Get("fpcandendz_6");
  h_pendz_bac[5]->Rebin(4);
  h_pendz_bac[5]->Sumw2();

  h_pendz_bac[6]=(TH1D*)input2->Get("fpcandendz_7");
  h_pendz_bac[6]->Rebin(4);
  h_pendz_bac[6]->Sumw2();
 
  h_pendz_bac[7]=(TH1D*)input2->Get("fpcandendz_8");
  h_pendz_bac[7]->Rebin(4);
  h_pendz_bac[7]->Sumw2();
 
  //=============================================================================  
  std::cout<<"libocheck1111111111111<<<<<<<<<<<<<<"<<std::endl; 
  ////~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  Double_t onbeam_data=h_range_allsel[0]->GetEntries();
  Double_t offbeam_data=h_range_allsel[1]->GetEntries();
  Double_t mc_entry=h_range_allsel[2]->GetEntries();


  Float_t tot_bnbcos=199800;
  Float_t tot_off=194747;
  Float_t tot_on= 185543;
  Float_t E1DCNT_wcut_bnb=10528552;
  Float_t E1DCNT_wcut_extbnb=12577618;
  

  Float_t int_per_1e20POT = 99035.2 ;

  //Float_t dataPOT = 0.495 * (tot_on) / 189130 ; // 0.495
  //Float_t mcbnbcos_POT = float(tot_bnbcos)/int_per_1e20POT ; // * dataPOT ; // 2.3 
  Float_t mcbnbcos_POT;
  if (tune==3){mcbnbcos_POT=4.0677e20;} // Tune1
  else{mcbnbcos_POT=2.0237e20;} // Tune1
  //Float_t mcbnbcos_POT=4.0677e20; // Tune3
  Float_t ondataPOT=4.714e19;
  Float_t offdataPOT=2.841e19;
  Float_t dataPOT=4.714e19;// ??????????????/

  Float_t RatioDatatoMC = dataPOT/mcbnbcos_POT ;
  Float_t RatioMCtoData = mcbnbcos_POT/dataPOT ;

  cout<<"RatioDatatoMC = "<<RatioDatatoMC<<endl;
  cout<<"RatioMCtoData = "<<RatioMCtoData<<endl;

  
  Double_t scalefac=E1DCNT_wcut_bnb/E1DCNT_wcut_extbnb;
  //Double_t scalefac=1.2300*(382718./tot_off)*(tot_on/547616.); //calculated by total number revents
  //Double_t scalefac=1.2300*(382718./17278997.)*(10686716./547616.); //calculated by total number of beam spills

  //Double_t normfac=(onbeam_data-offbeam_data*scalefac)/mc_entry;
  Double_t normfac=dataPOT/mcbnbcos_POT;
  cout<<"normalization factor for monte carlo sample is: "<<normfac<<endl;

  Double_t scale_onoffbeam=0.0;
  scale_onoffbeam=scalefac;
  cout<<"scale factor for extbnb is : "<<scale_onoffbeam<<endl;


//------------------------------------------------------------------
  TH1D *h_onoff_range=(TH1D*)h_range_allsel[1]->Clone(Form("%s_on-off", h_range_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;

  TCanvas *c1 = new TCanvas("c1"); 

  h_range_allsel[0]->SetLineColor(kBlack);
  h_range_allsel[0]->SetLineWidth(2);
  h_range_allsel[0]->SetLineStyle(1);
  h_range_allsel[0]->GetXaxis()->SetTitle("Track Length of Muon Candidate[cm]");
  h_range_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_range_allsel[0]->SetMaximum(500);
  h_range_allsel[0]->Draw();
  
  //h_range_allsel[1]->SetLineColor(kRed);
  //h_range_allsel[1]->SetLineWidth(2);
  //h_range_allsel[1]->SetLineStyle(1);

  h_onoff_range->Add(h_range_allsel[0],h_range_allsel[1],1,-scale_onoffbeam);
  h_onoff_range->SetLineColor(kBlack);
  h_onoff_range->SetLineWidth(2);
  h_onoff_range->SetLineStyle(1);
  h_onoff_range->GetXaxis()->SetTitle("Track Length of Muon Candidate[cm]");
  h_onoff_range->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_range->DrawCopy("");

  h_range_allsel[2]->SetLineColor(kRed);
  h_range_allsel[2]->SetLineWidth(2);
  h_range_allsel[2]->SetLineStyle(1);
  h_range_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_range = new THStack("hs_range","");

  h_range_sig[0]-> SetFillColor(2); 
  h_range_sig[0]->Scale(normfac);
  h_range_bac[0]-> SetFillColor(46);
  h_range_bac[0]->Scale(normfac);
  h_range_bac[1]-> SetFillColor(4);
  h_range_bac[1]->Scale(normfac);
  h_range_bac[2]-> SetFillColor(8);
  h_range_bac[2]->Scale(normfac);
  h_range_bac[3] -> SetFillColor(5);
  h_range_bac[3]->Scale(normfac);
  h_range_bac[4]-> SetFillColor(3);
  h_range_bac[4]->Scale(normfac);
  h_range_bac[5]-> SetFillColor(6);
  h_range_bac[5]->Scale(normfac);
  h_range_bac[6]-> SetFillColor(7);
  h_range_bac[6]->Scale(normfac);
  h_range_bac[7]-> SetFillColor(9);
  h_range_bac[7]->Scale(normfac);
  h_range_allsel[1]->SetFillStyle(3005);
  h_range_allsel[1]->SetFillColor(28);
  h_range_allsel[1]->Scale(scale_onoffbeam);


   hs_range -> Add(h_range_sig[0]);
   hs_range -> Add(h_range_bac[0]);
   hs_range -> Add(h_range_bac[1]);
   hs_range -> Add(h_range_bac[2]);
   hs_range -> Add(h_range_bac[3]);
   hs_range -> Add(h_range_bac[4]);
   hs_range -> Add(h_range_bac[5]);
   hs_range -> Add(h_range_bac[6]);
   hs_range -> Add(h_range_bac[7]);
   hs_range -> Add(h_range_allsel[1]);
   hs_range -> Draw("HIST,SAME");

   h_range_allsel[2]->Draw("same");
   h_range_allsel[0]->Draw("same");
   //h_onoff_range->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  TLegend *legend = new TLegend(.60, .65, .90, .90);

  legend->AddEntry(h_range_allsel[0], "on_beam");
  legend->AddEntry(h_range_allsel[1], "off_beam");
  //legend->AddEntry(h_onoff_range, "on beam - off beam");
  legend->AddEntry(h_range_allsel[2],"MC");
  legend -> AddEntry(h_range_sig[0], "1uNp (CC)", "f");
  legend -> AddEntry(h_range_bac[0], "CC0#pi0P", "f");
  legend -> AddEntry(h_range_bac[1], "CC1#piNP", "f");
  legend -> AddEntry(h_range_bac[2], "CCN#piNp", "f");
  legend -> AddEntry(h_range_bac[3], "CC#nu_{e}", "f");
  legend -> AddEntry(h_range_bac[4], "NC", "f");
  legend -> AddEntry(h_range_bac[5], "OOFV", "f");
  legend -> AddEntry(h_range_bac[6], "Cosmic", "f");
  legend -> AddEntry(h_range_bac[7], "Mixed", "f");
  legend->Draw("same");


  if(tune==3){c1->Print("figures/Tune3/BackSep/h_range_allsel.png");}
  else{c1->Print("figures/Tune1/BackSep/h_range_allsel.png");}

 //==============================================================================

  TH1D *h_onoff_trkmom=(TH1D*)h_trkmom_allsel[1]->Clone(Form("%s_on-off", h_trkmom_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_trkmom_allsel[0]->SetLineColor(kBlack);
  h_trkmom_allsel[0]->SetLineWidth(2);
  h_trkmom_allsel[0]->SetLineStyle(1);
  h_trkmom_allsel[0]->GetXaxis()->SetTitle("Track Momentum of Muon Candidate with MCS calculation[cm]");
  h_trkmom_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_trkmom_allsel[0]->SetMaximum(600);
  h_trkmom_allsel[0]->Draw();
  
  //h_trkmom_allsel[1]->SetLineColor(kRed);
  //h_trkmom_allsel[1]->SetLineWidth(2);
  //h_trkmom_allsel[1]->SetLineStyle(1);

  h_onoff_trkmom->Add(h_trkmom_allsel[0],h_trkmom_allsel[1],1,-scale_onoffbeam);
  h_onoff_trkmom->SetLineColor(kBlack);
  h_onoff_trkmom->SetLineWidth(2);
  h_onoff_trkmom->SetLineStyle(1);
  h_onoff_trkmom->GetXaxis()->SetTitle("Track Momentum of Muon Candidate with MCS calculation[cm]");
  h_onoff_trkmom->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_trkmom->DrawCopy("");

  h_trkmom_allsel[2]->SetLineColor(kRed);
  h_trkmom_allsel[2]->SetLineWidth(2);
  h_trkmom_allsel[2]->SetLineStyle(1);
  h_trkmom_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_trkmom = new THStack("hs_trkmom","");

  h_trkmom_sig[0]-> SetFillColor(2); 
  h_trkmom_sig[0]->Scale(normfac);
  h_trkmom_bac[0]-> SetFillColor(46);
  h_trkmom_bac[0]->Scale(normfac);
  h_trkmom_bac[1]-> SetFillColor(4);
  h_trkmom_bac[1]->Scale(normfac);
  h_trkmom_bac[2]-> SetFillColor(8);
  h_trkmom_bac[2]->Scale(normfac);
  h_trkmom_bac[3] -> SetFillColor(5);
  h_trkmom_bac[3]->Scale(normfac);
  h_trkmom_bac[4]-> SetFillColor(3);
  h_trkmom_bac[4]->Scale(normfac);
  h_trkmom_bac[5]-> SetFillColor(6);
  h_trkmom_bac[5]->Scale(normfac);
  h_trkmom_bac[6]-> SetFillColor(7);
  h_trkmom_bac[6]->Scale(normfac);
  h_trkmom_bac[7]-> SetFillColor(9);
  h_trkmom_bac[7]->Scale(normfac);
  h_trkmom_allsel[1]->SetFillStyle(3005);
  h_trkmom_allsel[1]->SetFillColor(28);
  h_trkmom_allsel[1]->Scale(scale_onoffbeam);



   hs_trkmom -> Add(h_trkmom_sig[0]);
   hs_trkmom -> Add(h_trkmom_bac[0]);
   hs_trkmom -> Add(h_trkmom_bac[1]);
   hs_trkmom -> Add(h_trkmom_bac[2]);
   hs_trkmom -> Add(h_trkmom_bac[3]);
   hs_trkmom -> Add(h_trkmom_bac[4]);
   hs_trkmom -> Add(h_trkmom_bac[5]);
   hs_trkmom -> Add(h_trkmom_bac[6]);
   hs_trkmom -> Add(h_trkmom_bac[7]);
   hs_trkmom -> Add(h_trkmom_allsel[1]);
   hs_trkmom -> Draw("HIST,SAME");

   h_trkmom_allsel[2]->Draw("same");
   h_trkmom_allsel[0]->Draw("same");
   //h_onoff_trkmom->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  
  legend->Draw("same");
  if(tune==3){c1->Print("figures/Tune3/BackSep/h_trkmom_allsel.png");}
  else{c1->Print("figures/Tune1/BackSep/h_trkmom_allsel.png");}

  //============================================================================

  TH1D *h_onoff_prange=(TH1D*)h_prange_allsel[1]->Clone(Form("%s_on-off", h_prange_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_prange_allsel[0]->SetLineColor(kBlack);
  h_prange_allsel[0]->SetLineWidth(2);
  h_prange_allsel[0]->SetLineStyle(1);
  h_prange_allsel[0]->GetXaxis()->SetTitle("Track Length of The Leading Proton Candidate[cm]");
  h_prange_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_prange_allsel[0]->SetMaximum(2.0*h_prange_allsel[0]->GetMaximum());
  h_prange_allsel[0]->Draw();
  //h_prange_allsel[1]->SetLineColor(kRed);
  //h_prange_allsel[1]->SetLineWidth(2);
  //h_prange_allsel[1]->SetLineStyle(1);

  h_onoff_prange->Add(h_prange_allsel[0],h_prange_allsel[1],1,-scale_onoffbeam);
  h_onoff_prange->SetLineColor(kBlack);
  h_onoff_prange->SetLineWidth(2);
  h_onoff_prange->SetLineStyle(1);
  h_onoff_prange->GetXaxis()->SetTitle("Track Length of The Leading Proton Candidate[cm]");
  h_onoff_prange->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_prange->DrawCopy("");

  h_prange_allsel[2]->SetLineColor(kRed);
  h_prange_allsel[2]->SetLineWidth(2);
  h_prange_allsel[2]->SetLineStyle(1);
  h_prange_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_prange = new THStack("hs_prange","");

  h_prange_sig[0]-> SetFillColor(2); 
  h_prange_sig[0]->Scale(normfac);
  h_prange_bac[0]-> SetFillColor(46);
  h_prange_bac[0]->Scale(normfac);
  h_prange_bac[1]-> SetFillColor(4);
  h_prange_bac[1]->Scale(normfac);
  h_prange_bac[2]-> SetFillColor(8);
  h_prange_bac[2]->Scale(normfac);
  h_prange_bac[3] -> SetFillColor(5);
  h_prange_bac[3]->Scale(normfac);
  h_prange_bac[4]-> SetFillColor(3);
  h_prange_bac[4]->Scale(normfac);
  h_prange_bac[5]-> SetFillColor(6);
  h_prange_bac[5]->Scale(normfac);
  h_prange_bac[6]-> SetFillColor(7);
  h_prange_bac[6]->Scale(normfac);
  h_prange_bac[7]-> SetFillColor(9);
  h_prange_bac[7]->Scale(normfac);
  h_prange_allsel[1]->SetFillStyle(3005);
  h_prange_allsel[1]->SetFillColor(28);
  h_prange_allsel[1]->Scale(scale_onoffbeam);



  hs_prange -> Add(h_prange_sig[0]);
  hs_prange -> Add(h_prange_bac[0]);
  hs_prange -> Add(h_prange_bac[1]);
  hs_prange -> Add(h_prange_bac[2]);
  hs_prange -> Add(h_prange_bac[3]);
  hs_prange -> Add(h_prange_bac[4]);
  hs_prange -> Add(h_prange_bac[5]);
  hs_prange -> Add(h_prange_bac[6]);
  hs_prange -> Add(h_prange_bac[7]);
  hs_prange -> Add(h_prange_allsel[1]);
  hs_prange -> Draw("HIST,SAME");

  //h_prange_allsel[2]->Draw("same"); // MC error bars, but done wrong
  h_prange_allsel[0]->Draw("same");
  //h_onoff_prange->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legend->Draw("same");
  if(tune==3){c1->Print("figures/Tune3/BackSep/h_prange_allsel.png");}
  else{c1->Print("figures/Tune1/BackSep/h_prange_allsel.png");}

  //=====================ZOOM in on proton length ============
  THStack *hs_prange_zoom = new THStack("hs_prange_zoom","");
  std::cout << "allzoom[0] " << h_prange_allzoom[0]->GetNbinsX() << std::endl;
  h_prange_allzoom[0]->SetLineColor(kBlack);
  h_prange_allzoom[0]->SetLineWidth(2);
  h_prange_allzoom[0]->SetLineStyle(1);
  h_prange_allzoom[0]->GetXaxis()->SetTitle("Track Length of The Leading Proton Candidate[cm]");
  h_prange_allzoom[0]->GetYaxis()->SetTitle("No. of Tracks");
  std::cout << "allzoom[0] " << h_prange_allzoom[0]->GetNbinsX() << std::endl;
  h_prange_allzoom[0]->SetMaximum(2.0*h_prange_allzoom[0]->GetMaximum());
  h_prange_allzoom[0]->GetXaxis()->SetRangeUser(0,10);
  h_prange_allzoom[0]->Draw();
  std::cout << "allzoom[0] " << h_prange_allzoom[0]->GetNbinsX() << std::endl;
  
  h_prange_sigzoom[0]-> SetFillColor(2); 
  h_prange_sigzoom[0]->Scale(normfac);
  h_prange_baczoom[0]-> SetFillColor(46);
  h_prange_baczoom[0]->Scale(normfac);
  h_prange_baczoom[1]-> SetFillColor(4);
  h_prange_baczoom[1]->Scale(normfac);
  h_prange_baczoom[2]-> SetFillColor(8);
  h_prange_baczoom[2]->Scale(normfac);
  h_prange_baczoom[3] -> SetFillColor(5);
  h_prange_baczoom[3]->Scale(normfac);
  h_prange_baczoom[4]-> SetFillColor(3);
  h_prange_baczoom[4]->Scale(normfac);
  h_prange_baczoom[5]-> SetFillColor(6);
  h_prange_baczoom[5]->Scale(normfac);
  h_prange_baczoom[6]-> SetFillColor(7);
  h_prange_baczoom[6]->Scale(normfac);
  h_prange_baczoom[7]-> SetFillColor(9);
  h_prange_baczoom[7]->Scale(normfac);
  h_prange_allzoom[1]->SetFillStyle(3005);
  h_prange_allzoom[1]->SetFillColor(28);
  h_prange_allzoom[1]->Scale(scale_onoffbeam);

  hs_prange_zoom -> Add(h_prange_sigzoom[0]);
  hs_prange_zoom -> Add(h_prange_baczoom[0]);
  hs_prange_zoom -> Add(h_prange_baczoom[1]);
  hs_prange_zoom -> Add(h_prange_baczoom[2]);
  hs_prange_zoom -> Add(h_prange_baczoom[3]);
  hs_prange_zoom -> Add(h_prange_baczoom[4]);
  hs_prange_zoom -> Add(h_prange_baczoom[5]);
  hs_prange_zoom -> Add(h_prange_baczoom[6]);
  hs_prange_zoom -> Add(h_prange_baczoom[7]);
  hs_prange_zoom -> Add(h_prange_allzoom[1]);
  hs_prange_zoom -> Draw("HIST SAME");
//  hs_prange_zoom -> GetXaxis()->SetRangeUser(0,10);
//  hs_prange_zoom -> Draw("HIST");

//  h_prange_zoom[2]->Draw("same");
  std::cout << "allzoom[0] " << h_prange_allzoom[0]->GetNbinsX() << ", h_prange_sigzoom[0] " << h_prange_sigzoom[0]->GetNbinsX() << std::endl;
  h_prange_allzoom[0]->Draw("same");
  
  if(tune==3){c1->Print("figures/Tune3/BackSep/h_prange_zoom.png");}
  else{c1->Print("figures/Tune1/BackSep/h_prange_zoom.png");}
  
  //===================================================================
  TH1D *h_onoff_costheta=(TH1D*)h_costheta_allsel[1]->Clone(Form("%s_on-off", h_costheta_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_costheta_allsel[0]->SetLineColor(kBlack);
  h_costheta_allsel[0]->SetLineWidth(2);
  h_costheta_allsel[0]->SetLineStyle(1);
  h_costheta_allsel[0]->GetXaxis()->SetTitle("cos#theta_{#mu}");
  h_costheta_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_costheta_allsel[0]->SetMaximum(800);
  h_costheta_allsel[0]->Draw();
  //h_costheta_allsel[1]->SetLineColor(kRed);
  //h_costheta_allsel[1]->SetLineWidth(2);
  //h_costheta_allsel[1]->SetLineStyle(1);

  h_onoff_costheta->Add(h_costheta_allsel[0],h_costheta_allsel[1],1,-scale_onoffbeam);
  h_onoff_costheta->SetLineColor(kBlack);
  h_onoff_costheta->SetLineWidth(2);
  h_onoff_costheta->SetLineStyle(1);
  h_onoff_costheta->GetXaxis()->SetTitle("cos#theta_{#mu}");
  h_onoff_costheta->GetYaxis()->SetTitle("No. of Tracks");
  h_onoff_costheta->SetMaximum(500);
  //h_onoff_costheta->DrawCopy("");

  h_costheta_allsel[2]->SetLineColor(kRed);
  h_costheta_allsel[2]->SetLineWidth(2);
  h_costheta_allsel[2]->SetLineStyle(1);
  h_costheta_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_costheta = new THStack("hs_costheta","");

  h_costheta_sig[0]-> SetFillColor(2); 
  h_costheta_sig[0]->Scale(normfac);
  h_costheta_bac[0]-> SetFillColor(46);
  h_costheta_bac[0]->Scale(normfac);
  h_costheta_bac[1]-> SetFillColor(4);
  h_costheta_bac[1]->Scale(normfac);
  h_costheta_bac[2]-> SetFillColor(8);
  h_costheta_bac[2]->Scale(normfac);
  h_costheta_bac[3] -> SetFillColor(5);
  h_costheta_bac[3]->Scale(normfac);
  h_costheta_bac[4]-> SetFillColor(3);
  h_costheta_bac[4]->Scale(normfac);
  h_costheta_bac[5]-> SetFillColor(6);
  h_costheta_bac[5]->Scale(normfac);
  h_costheta_bac[6]-> SetFillColor(7);
  h_costheta_bac[6]->Scale(normfac);
  h_costheta_bac[7]-> SetFillColor(9);
  h_costheta_bac[7]->Scale(normfac);
  h_costheta_allsel[1]->SetFillStyle(3005);
  h_costheta_allsel[1]->SetFillColor(28);
  h_costheta_allsel[1]->Scale(scale_onoffbeam);



   hs_costheta -> Add(h_costheta_sig[0]);
   hs_costheta -> Add(h_costheta_bac[0]);
   hs_costheta -> Add(h_costheta_bac[1]);
   hs_costheta -> Add(h_costheta_bac[2]);
   hs_costheta -> Add(h_costheta_bac[3]);
   hs_costheta -> Add(h_costheta_bac[4]);
   hs_costheta -> Add(h_costheta_bac[5]);
   hs_costheta -> Add(h_costheta_bac[6]);
   hs_costheta -> Add(h_costheta_bac[7]);
   hs_costheta -> Add(h_costheta_allsel[1]);
   hs_costheta -> Draw("HIST,SAME");

   h_costheta_allsel[2]->Draw("same");
   h_costheta_allsel[0]->Draw("same");
   //h_onoff_costheta->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  legend->Draw("same");
  if(tune==3){c1->Print("figures/Tune3/BackSep/h_costheta_allsel.png");}
  else{c1->Print("figures/Tune1/BackSep/h_costheta_allsel.png");}

  //================================================================

  TH1D *h_onoff_pcostheta=(TH1D*)h_pcostheta_allsel[1]->Clone(Form("%s_on-off", h_pcostheta_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_pcostheta_allsel[0]->SetLineColor(kBlack);
  h_pcostheta_allsel[0]->SetLineWidth(2);
  h_pcostheta_allsel[0]->SetLineStyle(1);
  h_pcostheta_allsel[0]->GetXaxis()->SetTitle("cos#theta_{p}");
  h_pcostheta_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_pcostheta_allsel[0]->SetMaximum(500);
   h_pcostheta_allsel[0]->Draw();

  //h_pcostheta_allsel[1]->SetLineColor(kRed);
  //h_pcostheta_allsel[1]->SetLineWidth(2);
  //h_pcostheta_allsel[1]->SetLineStyle(1);

  h_onoff_pcostheta->Add(h_pcostheta_allsel[0],h_pcostheta_allsel[1],1,-scale_onoffbeam);
  h_onoff_pcostheta->SetLineColor(kBlack);
  h_onoff_pcostheta->SetLineWidth(2);
  h_onoff_pcostheta->SetLineStyle(1);
  h_onoff_pcostheta->SetMaximum(300);
  h_onoff_pcostheta->GetXaxis()->SetTitle("cos#theta_{p}");
  h_onoff_pcostheta->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_pcostheta->DrawCopy("");

  h_pcostheta_allsel[2]->SetLineColor(kRed);
  h_pcostheta_allsel[2]->SetLineWidth(2);
  h_pcostheta_allsel[2]->SetLineStyle(1);
  h_pcostheta_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_pcostheta = new THStack("hs_pcostheta","");

  h_pcostheta_sig[0]-> SetFillColor(2); 
  h_pcostheta_sig[0]->Scale(normfac);
  h_pcostheta_bac[0]-> SetFillColor(46);
  h_pcostheta_bac[0]->Scale(normfac);
  h_pcostheta_bac[1]-> SetFillColor(4);
  h_pcostheta_bac[1]->Scale(normfac);
  h_pcostheta_bac[2]-> SetFillColor(8);
  h_pcostheta_bac[2]->Scale(normfac);
  h_pcostheta_bac[3] -> SetFillColor(5);
  h_pcostheta_bac[3]->Scale(normfac);
  h_pcostheta_bac[4]-> SetFillColor(3);
  h_pcostheta_bac[4]->Scale(normfac);
  h_pcostheta_bac[5]-> SetFillColor(6);
  h_pcostheta_bac[5]->Scale(normfac);
  h_pcostheta_bac[6]-> SetFillColor(7);
  h_pcostheta_bac[6]->Scale(normfac);
  h_pcostheta_bac[7]-> SetFillColor(9);
  h_pcostheta_bac[7]->Scale(normfac);
  h_pcostheta_allsel[1]->SetFillStyle(3005);
  h_pcostheta_allsel[1]->SetFillColor(28);
  h_pcostheta_allsel[1]->Scale(scale_onoffbeam);



   hs_pcostheta -> Add(h_pcostheta_sig[0]);
   hs_pcostheta -> Add(h_pcostheta_bac[0]);
   hs_pcostheta -> Add(h_pcostheta_bac[1]);
   hs_pcostheta -> Add(h_pcostheta_bac[2]);
   hs_pcostheta -> Add(h_pcostheta_bac[3]);
   hs_pcostheta -> Add(h_pcostheta_bac[4]);
   hs_pcostheta -> Add(h_pcostheta_bac[5]);
   hs_pcostheta -> Add(h_pcostheta_bac[6]);
   hs_pcostheta -> Add(h_pcostheta_bac[7]);
   hs_pcostheta -> Add(h_pcostheta_allsel[1]);
   hs_pcostheta -> Draw("HIST,SAME");

   h_pcostheta_allsel[2]->Draw("same");
   h_pcostheta_allsel[0]->Draw("same");
   //h_onoff_pcostheta->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   legend->Draw("same");
  if(tune==3){c1->Print("figures/Tune3/BackSep/h_pcostheta_allsel.png");}
  else{c1->Print("figures/Tune1/BackSep/h_pcostheta_allsel.png");}

  //==================================================================================

  TH1D *h_onoff_trunmean_muon=(TH1D*)h_trunmean_muon_allsel[1]->Clone(Form("%s_on-off", h_trunmean_muon_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_trunmean_muon_allsel[0]->SetLineColor(kBlack);
  h_trunmean_muon_allsel[0]->SetLineWidth(2);
  h_trunmean_muon_allsel[0]->SetLineStyle(1);
  h_trunmean_muon_allsel[0]->GetXaxis()->SetTitle("Truncated Mean dQdx of Muon Candidate");
  h_trunmean_muon_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_trunmean_muon_allsel[0]->SetMaximum(1500);
   h_trunmean_muon_allsel[0]->Draw();

  //h_trunmean_muon_allsel[1]->SetLineColor(kRed);
  //h_trunmean_muon_allsel[1]->SetLineWidth(2);
  //h_trunmean_muon_allsel[1]->SetLineStyle(1);

  h_onoff_trunmean_muon->Add(h_trunmean_muon_allsel[0],h_trunmean_muon_allsel[1],1,-scale_onoffbeam);
  h_onoff_trunmean_muon->SetLineColor(kBlack);
  h_onoff_trunmean_muon->SetLineWidth(2);
  h_onoff_trunmean_muon->SetLineStyle(1);
  h_onoff_trunmean_muon->SetMaximum(1500);
  h_onoff_trunmean_muon->GetXaxis()->SetTitle("Truncated Mean dQdx of Muon Candidate");
  h_onoff_trunmean_muon->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_trunmean_muon->DrawCopy("");



  h_trunmean_muon_allsel[2]->SetLineColor(kRed);
  h_trunmean_muon_allsel[2]->SetLineWidth(2);
  h_trunmean_muon_allsel[2]->SetLineStyle(1);
  h_trunmean_muon_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_trunmean_muon = new THStack("hs_trunmean_muon","");

  h_trunmean_muon_sig[0]-> SetFillColor(2); 
  h_trunmean_muon_sig[0]->Scale(normfac);
  h_trunmean_muon_bac[0]-> SetFillColor(46);
  h_trunmean_muon_bac[0]->Scale(normfac);
  h_trunmean_muon_bac[1]-> SetFillColor(4);
  h_trunmean_muon_bac[1]->Scale(normfac);
  h_trunmean_muon_bac[2]-> SetFillColor(8);
  h_trunmean_muon_bac[2]->Scale(normfac);
  h_trunmean_muon_bac[3] -> SetFillColor(5);
  h_trunmean_muon_bac[3]->Scale(normfac);
  h_trunmean_muon_bac[4]-> SetFillColor(3);
  h_trunmean_muon_bac[4]->Scale(normfac);
  h_trunmean_muon_bac[5]-> SetFillColor(6);
  h_trunmean_muon_bac[5]->Scale(normfac);
  h_trunmean_muon_bac[6]-> SetFillColor(7);
  h_trunmean_muon_bac[6]->Scale(normfac);
  h_trunmean_muon_bac[7]-> SetFillColor(9);
  h_trunmean_muon_bac[7]->Scale(normfac);
  h_trunmean_muon_allsel[1]->SetFillStyle(3005);
  h_trunmean_muon_allsel[1]->SetFillColor(28);
  h_trunmean_muon_allsel[1]->Scale(scale_onoffbeam);



   hs_trunmean_muon -> Add(h_trunmean_muon_sig[0]);
   hs_trunmean_muon -> Add(h_trunmean_muon_bac[0]);
   hs_trunmean_muon -> Add(h_trunmean_muon_bac[1]);
   hs_trunmean_muon -> Add(h_trunmean_muon_bac[2]);
   hs_trunmean_muon -> Add(h_trunmean_muon_bac[3]);
   hs_trunmean_muon -> Add(h_trunmean_muon_bac[4]);
   hs_trunmean_muon -> Add(h_trunmean_muon_bac[5]);
   hs_trunmean_muon -> Add(h_trunmean_muon_bac[6]);
   hs_trunmean_muon -> Add(h_trunmean_muon_bac[7]);
   hs_trunmean_muon -> Add(h_trunmean_muon_allsel[1]);
   hs_trunmean_muon -> Draw("HIST,SAME");

   h_trunmean_muon_allsel[2]->Draw("same");
   h_trunmean_muon_allsel[0]->Draw("same");
   //h_onoff_trunmean_muon->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_trunmean_muon.png");}
   else{c1->Print("figures/Tune1/BackSep/h_trunmean_muon.png");}
  //==================================================================================

  TH1D *h_onoff_trunmean_proton=(TH1D*)h_trunmean_proton_allsel[1]->Clone(Form("%s_on-off", h_trunmean_proton_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_trunmean_proton_allsel[0]->SetLineColor(kBlack);
  h_trunmean_proton_allsel[0]->SetLineWidth(2);
  h_trunmean_proton_allsel[0]->SetLineStyle(1);
  h_trunmean_proton_allsel[0]->GetXaxis()->SetTitle("Truncated Mean dQdx of Proton Candidate");
  h_trunmean_proton_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_trunmean_proton_allsel[0]->SetMaximum(700);
   h_trunmean_proton_allsel[0]->Draw();

  //h_trunmean_proton_allsel[1]->SetLineColor(kRed);
  //h_trunmean_proton_allsel[1]->SetLineWidth(2);
  //h_trunmean_proton_allsel[1]->SetLineStyle(1);

  h_onoff_trunmean_proton->Add(h_trunmean_proton_allsel[0],h_trunmean_proton_allsel[1],1,-scale_onoffbeam);
  h_onoff_trunmean_proton->SetLineColor(kBlack);
  h_onoff_trunmean_proton->SetLineWidth(2);
  h_onoff_trunmean_proton->SetLineStyle(1);
  h_onoff_trunmean_proton->SetMaximum(700);
  h_onoff_trunmean_proton->GetXaxis()->SetTitle("Truncated Mean dQdx of Muon Candidate");
  h_onoff_trunmean_proton->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_trunmean_proton->DrawCopy("");



  h_trunmean_proton_allsel[2]->SetLineColor(kRed);
  h_trunmean_proton_allsel[2]->SetLineWidth(2);
  h_trunmean_proton_allsel[2]->SetLineStyle(1);
  h_trunmean_proton_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_trunmean_proton = new THStack("hs_trunmean_proton","");

  h_trunmean_proton_sig[0]-> SetFillColor(2); 
  h_trunmean_proton_sig[0]->Scale(normfac);
  h_trunmean_proton_bac[0]-> SetFillColor(46);
  h_trunmean_proton_bac[0]->Scale(normfac);
  h_trunmean_proton_bac[1]-> SetFillColor(4);
  h_trunmean_proton_bac[1]->Scale(normfac);
  h_trunmean_proton_bac[2]-> SetFillColor(8);
  h_trunmean_proton_bac[2]->Scale(normfac);
  h_trunmean_proton_bac[3] -> SetFillColor(5);
  h_trunmean_proton_bac[3]->Scale(normfac);
  h_trunmean_proton_bac[4]-> SetFillColor(3);
  h_trunmean_proton_bac[4]->Scale(normfac);
  h_trunmean_proton_bac[5]-> SetFillColor(6);
  h_trunmean_proton_bac[5]->Scale(normfac);
  h_trunmean_proton_bac[6]-> SetFillColor(7);
  h_trunmean_proton_bac[6]->Scale(normfac);
  h_trunmean_proton_bac[7]-> SetFillColor(9);
  h_trunmean_proton_bac[7]->Scale(normfac);
  h_trunmean_proton_allsel[1]->SetFillStyle(3005);
  h_trunmean_proton_allsel[1]->SetFillColor(28);
  h_trunmean_proton_allsel[1]->Scale(scale_onoffbeam);



   hs_trunmean_proton -> Add(h_trunmean_proton_sig[0]);
   hs_trunmean_proton -> Add(h_trunmean_proton_bac[0]);
   hs_trunmean_proton -> Add(h_trunmean_proton_bac[1]);
   hs_trunmean_proton -> Add(h_trunmean_proton_bac[2]);
   hs_trunmean_proton -> Add(h_trunmean_proton_bac[3]);
   hs_trunmean_proton -> Add(h_trunmean_proton_bac[4]);
   hs_trunmean_proton -> Add(h_trunmean_proton_bac[5]);
   hs_trunmean_proton -> Add(h_trunmean_proton_bac[6]);
   hs_trunmean_proton -> Add(h_trunmean_proton_bac[7]);
   hs_trunmean_proton -> Add(h_trunmean_proton_allsel[1]);
   hs_trunmean_proton -> Draw("HIST,SAME");

   h_trunmean_proton_allsel[2]->Draw("same");
   h_trunmean_proton_allsel[0]->Draw("same");
   //h_onoff_trunmean_proton->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_trunmean_proton.png");}
   else{c1->Print("figures/Tune1/BackSep/h_trunmean_proton.png");}
  
//===============================================

  TH1D *h_onoff_phi=(TH1D*)h_phi_allsel[1]->Clone(Form("%s_on-off", h_phi_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_phi_allsel[0]->SetLineColor(kBlack);
  h_phi_allsel[0]->SetLineWidth(2);
  h_phi_allsel[0]->SetLineStyle(1);
  h_phi_allsel[0]->GetXaxis()->SetTitle("#phi_{#mu}[Rad]");
  h_phi_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_phi_allsel[0]->SetMaximum(240);
  h_phi_allsel[0]->SetMinimum(0);
  h_phi_allsel[0]->Draw();
  
  //h_phi_allsel[1]->SetLineColor(kRed);
  //h_phi_allsel[1]->SetLineWidth(2);
  //h_phi_allsel[1]->SetLineStyle(1);

  h_onoff_phi->Add(h_phi_allsel[0],h_phi_allsel[1],1,-scale_onoffbeam);
  h_onoff_phi->SetLineColor(kBlack);
  h_onoff_phi->SetLineWidth(2);
  h_onoff_phi->SetLineStyle(1);
  h_onoff_phi->SetMaximum(80);
  h_onoff_phi->GetXaxis()->SetTitle("#phi_{#mu}[Rad]");
  h_onoff_phi->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_phi->DrawCopy("");

  h_phi_allsel[2]->SetLineColor(kRed);
  h_phi_allsel[2]->SetLineWidth(2);
  h_phi_allsel[2]->SetLineStyle(1);
  h_phi_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_phi = new THStack("hs_phi","");

  h_phi_sig[0]-> SetFillColor(2); 
  h_phi_sig[0]->Scale(normfac);
  h_phi_bac[0]-> SetFillColor(46);
  h_phi_bac[0]->Scale(normfac);
  h_phi_bac[1]-> SetFillColor(4);
  h_phi_bac[1]->Scale(normfac);
  h_phi_bac[2]-> SetFillColor(8);
  h_phi_bac[2]->Scale(normfac);
  h_phi_bac[3] -> SetFillColor(5);
  h_phi_bac[3]->Scale(normfac);
  h_phi_bac[4]-> SetFillColor(3);
  h_phi_bac[4]->Scale(normfac);
  h_phi_bac[5]-> SetFillColor(6);
  h_phi_bac[5]->Scale(normfac);
  h_phi_bac[6]-> SetFillColor(7);
  h_phi_bac[6]->Scale(normfac);
  h_phi_bac[7]-> SetFillColor(9);
  h_phi_bac[7]->Scale(normfac);
  h_phi_allsel[1]->SetFillStyle(3005);
  h_phi_allsel[1]->SetFillColor(28);
  h_phi_allsel[1]->Scale(scale_onoffbeam);



   hs_phi -> Add(h_phi_sig[0]);
   hs_phi -> Add(h_phi_bac[0]);
   hs_phi -> Add(h_phi_bac[1]);
   hs_phi -> Add(h_phi_bac[2]);
   hs_phi -> Add(h_phi_bac[3]);
   hs_phi -> Add(h_phi_bac[4]);
   hs_phi -> Add(h_phi_bac[5]);
   hs_phi -> Add(h_phi_bac[6]);
   hs_phi -> Add(h_phi_bac[7]);
   hs_phi -> Add(h_phi_allsel[1]);
   hs_phi -> Draw("HIST,SAME");

   h_phi_allsel[2]->Draw("same");
   h_phi_allsel[0]->Draw("same");
   //h_onoff_phi->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_phi_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_phi_allsel.png");}

  //=========================================================================================  
  TH1D *h_onoff_pphi=(TH1D*)h_pphi_allsel[1]->Clone(Form("%s_on-off", h_pphi_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_pphi_allsel[0]->SetLineColor(kBlack);
  h_pphi_allsel[0]->SetLineWidth(2);
  h_pphi_allsel[0]->SetLineStyle(1);
  h_pphi_allsel[0]->GetXaxis()->SetTitle("#phi_{P}[Rad]");
  h_pphi_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_pphi_allsel[0]->SetMaximum(200);
  h_pphi_allsel[0]->SetMinimum(0);
  h_pphi_allsel[0]->Draw();

  //h_pphi_allsel[1]->SetLineColor(kRed);
  //h_pphi_allsel[1]->SetLineWidth(2);
  //h_pphi_allsel[1]->SetLineStyle(1);

  h_onoff_pphi->Add(h_pphi_allsel[0],h_pphi_allsel[1],1,-scale_onoffbeam);
  h_onoff_pphi->SetLineColor(kBlack);
  h_onoff_pphi->SetLineWidth(2);
  h_onoff_pphi->SetLineStyle(1);
  h_onoff_pphi->SetMaximum(70);
  h_onoff_pphi->SetMinimum(0);
  h_onoff_pphi->GetXaxis()->SetTitle("#phi_{P}[Rad]");
  h_onoff_pphi->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_pphi->DrawCopy("");

  h_pphi_allsel[2]->SetLineColor(kRed);
  h_pphi_allsel[2]->SetLineWidth(2);
  h_pphi_allsel[2]->SetLineStyle(1);
  h_pphi_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_pphi = new THStack("hs_pphi","");

  h_pphi_sig[0]-> SetFillColor(2); 
  h_pphi_sig[0]->Scale(normfac);
  h_pphi_bac[0]-> SetFillColor(46);
  h_pphi_bac[0]->Scale(normfac);
  h_pphi_bac[1]-> SetFillColor(4);
  h_pphi_bac[1]->Scale(normfac);
  h_pphi_bac[2]-> SetFillColor(8);
  h_pphi_bac[2]->Scale(normfac);
  h_pphi_bac[3] -> SetFillColor(5);
  h_pphi_bac[3]->Scale(normfac);
  h_pphi_bac[4]-> SetFillColor(3);
  h_pphi_bac[4]->Scale(normfac);
  h_pphi_bac[5]-> SetFillColor(6);
  h_pphi_bac[5]->Scale(normfac);
  h_pphi_bac[6]-> SetFillColor(7);
  h_pphi_bac[6]->Scale(normfac);
  h_pphi_bac[7]-> SetFillColor(9);
  h_pphi_bac[7]->Scale(normfac);
  h_pphi_allsel[1]->SetFillStyle(3005);
  h_pphi_allsel[1]->SetFillColor(28);
  h_pphi_allsel[1]->Scale(scale_onoffbeam);



   hs_pphi -> Add(h_pphi_sig[0]);
   hs_pphi -> Add(h_pphi_bac[0]);
   hs_pphi -> Add(h_pphi_bac[1]);
   hs_pphi -> Add(h_pphi_bac[2]);
   hs_pphi -> Add(h_pphi_bac[3]);
   hs_pphi -> Add(h_pphi_bac[4]);
   hs_pphi -> Add(h_pphi_bac[5]);
   hs_pphi -> Add(h_pphi_bac[6]);
   hs_pphi -> Add(h_pphi_bac[7]);
   hs_pphi -> Add(h_pphi_allsel[1]);
   hs_pphi -> Draw("HIST,SAME");

   h_pphi_allsel[2]->Draw("same");
   h_pphi_allsel[0]->Draw("same");
   //h_onoff_pphi->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_pphi_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_pphi_allsel.png");}
 //================================================================================ 
  TH1D *h_onoff_ntrksp=(TH1D*)h_ntrksp_allsel[1]->Clone(Form("%s_on-off", h_ntrksp_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_ntrksp_allsel[0]->SetLineColor(kBlack);
  h_ntrksp_allsel[0]->SetLineWidth(2);
  h_ntrksp_allsel[0]->SetLineStyle(1);
  h_ntrksp_allsel[0]->GetXaxis()->SetTitle("Number of The Proton Candidates");
  h_ntrksp_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_ntrksp_allsel[0]->SetMaximum(2000);
  h_ntrksp_allsel[0]->Draw();

  //h_ntrksp_allsel[1]->SetLineColor(kRed);
  //h_ntrksp_allsel[1]->SetLineWidth(2);
  //h_ntrksp_allsel[1]->SetLineStyle(1);

  h_onoff_ntrksp->Add(h_ntrksp_allsel[0],h_ntrksp_allsel[1],1,-scale_onoffbeam);
  h_onoff_ntrksp->SetLineColor(kBlack);
  h_onoff_ntrksp->SetLineWidth(2);
  h_onoff_ntrksp->SetLineStyle(1);
  h_onoff_ntrksp->GetXaxis()->SetTitle("Number of The Proton Candidates");
  h_onoff_ntrksp->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_ntrksp->DrawCopy("");

  h_ntrksp_allsel[2]->SetLineColor(kRed);
  h_ntrksp_allsel[2]->SetLineWidth(2);
  h_ntrksp_allsel[2]->SetLineStyle(1);
  h_ntrksp_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_ntrksp = new THStack("hs_ntrksp","");

  h_ntrksp_sig[0]-> SetFillColor(2); 
  h_ntrksp_sig[0]->Scale(normfac);
  h_ntrksp_bac[0]-> SetFillColor(46);
  h_ntrksp_bac[0]->Scale(normfac);
  h_ntrksp_bac[1]-> SetFillColor(4);
  h_ntrksp_bac[1]->Scale(normfac);
  h_ntrksp_bac[2]-> SetFillColor(8);
  h_ntrksp_bac[2]->Scale(normfac);
  h_ntrksp_bac[3] -> SetFillColor(5);
  h_ntrksp_bac[3]->Scale(normfac);
  h_ntrksp_bac[4]-> SetFillColor(3);
  h_ntrksp_bac[4]->Scale(normfac);
  h_ntrksp_bac[5]-> SetFillColor(6);
  h_ntrksp_bac[5]->Scale(normfac);
  h_ntrksp_bac[6]-> SetFillColor(7);
  h_ntrksp_bac[6]->Scale(normfac);
  h_ntrksp_bac[7]-> SetFillColor(9);
  h_ntrksp_bac[7]->Scale(normfac);
  h_ntrksp_allsel[1]->SetFillStyle(3005);
  h_ntrksp_allsel[1]->SetFillColor(28);
  h_ntrksp_allsel[1]->Scale(scale_onoffbeam);



   hs_ntrksp -> Add(h_ntrksp_sig[0]);
   hs_ntrksp -> Add(h_ntrksp_bac[0]);
   hs_ntrksp -> Add(h_ntrksp_bac[1]);
   hs_ntrksp -> Add(h_ntrksp_bac[2]);
   hs_ntrksp -> Add(h_ntrksp_bac[3]);
   hs_ntrksp -> Add(h_ntrksp_bac[4]);
   hs_ntrksp -> Add(h_ntrksp_bac[5]);
   hs_ntrksp -> Add(h_ntrksp_bac[6]);
   hs_ntrksp -> Add(h_ntrksp_bac[7]);
   hs_ntrksp -> Add(h_ntrksp_allsel[1]);
   hs_ntrksp -> Draw("HIST,SAME");

   h_ntrksp_allsel[2]->Draw("same");
   h_ntrksp_allsel[0]->Draw("same");
   //h_onoff_ntrksp->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_ntrksp_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_ntrksp_allsel.png");}
//===========================================================================
  TH1D *h_onoff_plep=(TH1D*)h_plep_allsel[1]->Clone(Form("%s_on-off", h_plep_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_plep_allsel[0]->SetLineColor(kBlack);
  h_plep_allsel[0]->SetLineWidth(2);
  h_plep_allsel[0]->SetLineStyle(1);
  h_plep_allsel[0]->GetXaxis()->SetTitle("Momentum of the muon candidate[GeV]");
  h_plep_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_plep_allsel[0]->SetMaximum(600);
  h_plep_allsel[0]->Draw();  

  //h_plep_allsel[1]->SetLineColor(kRed);
  //h_plep_allsel[1]->SetLineWidth(2);
  //h_plep_allsel[1]->SetLineStyle(1);

  h_onoff_plep->Add(h_plep_allsel[0],h_plep_allsel[1],1,-scale_onoffbeam);
  h_onoff_plep->SetLineColor(kBlack);
  h_onoff_plep->SetLineWidth(2);
  h_onoff_plep->SetLineStyle(1);
  h_onoff_plep->GetXaxis()->SetTitle("Momentum of the muon candidate[GeV]");
  h_onoff_plep->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_plep->DrawCopy("");

  h_plep_allsel[2]->SetLineColor(kRed);
  h_plep_allsel[2]->SetLineWidth(2);
  h_plep_allsel[2]->SetLineStyle(1);
  h_plep_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_plep = new THStack("hs_plep","");

  h_plep_sig[0]-> SetFillColor(2); 
  h_plep_sig[0]->Scale(normfac);
  h_plep_bac[0]-> SetFillColor(46);
  h_plep_bac[0]->Scale(normfac);
  h_plep_bac[1]-> SetFillColor(4);
  h_plep_bac[1]->Scale(normfac);
  h_plep_bac[2]-> SetFillColor(8);
  h_plep_bac[2]->Scale(normfac);
  h_plep_bac[3] -> SetFillColor(5);
  h_plep_bac[3]->Scale(normfac);
  h_plep_bac[4]-> SetFillColor(3);
  h_plep_bac[4]->Scale(normfac);
  h_plep_bac[5]-> SetFillColor(6);
  h_plep_bac[5]->Scale(normfac);
  h_plep_bac[6]-> SetFillColor(7);
  h_plep_bac[6]->Scale(normfac);
  h_plep_bac[7]-> SetFillColor(9);
  h_plep_bac[7]->Scale(normfac);
  h_plep_allsel[1]->SetFillStyle(3005);
  h_plep_allsel[1]->SetFillColor(28);
  h_plep_allsel[1]->Scale(scale_onoffbeam);



   hs_plep -> Add(h_plep_sig[0]);
   hs_plep -> Add(h_plep_bac[0]);
   hs_plep -> Add(h_plep_bac[1]);
   hs_plep -> Add(h_plep_bac[2]);
   hs_plep -> Add(h_plep_bac[3]);
   hs_plep -> Add(h_plep_bac[4]);
   hs_plep -> Add(h_plep_bac[5]);
   hs_plep -> Add(h_plep_bac[6]);
   hs_plep -> Add(h_plep_bac[7]);
   hs_plep -> Add(h_plep_allsel[1]);
   hs_plep -> Draw("HIST,SAME");

   h_plep_allsel[2]->Draw("same");
   h_plep_allsel[0]->Draw("same");  
   //h_onoff_plep->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_plep_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_plep_allsel.png");}
 //==============================================================

  TH1D *h_onoff_phad=(TH1D*)h_phad_allsel[1]->Clone(Form("%s_on-off", h_phad_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_phad_allsel[0]->SetLineColor(kBlack);
  h_phad_allsel[0]->SetLineWidth(2);
  h_phad_allsel[0]->SetLineStyle(1);
  h_phad_allsel[0]->GetXaxis()->SetTitle("Momentum of the leading proton candidates[GeV]");
  h_phad_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_phad_allsel[0]->SetMaximum(600);
  h_phad_allsel[0]->Draw();  

  //h_phad_allsel[1]->SetLineColor(kRed);
  //h_phad_allsel[1]->SetLineWidth(2);
  //h_phad_allsel[1]->SetLineStyle(1);

  h_onoff_phad->Add(h_phad_allsel[0],h_phad_allsel[1],1,-scale_onoffbeam);
  h_onoff_phad->SetLineColor(kBlack);
  h_onoff_phad->SetLineWidth(2);
  h_onoff_phad->SetLineStyle(1);
  h_onoff_phad->GetXaxis()->SetTitle("Momentum of the leading proton candidates[GeV]");
  h_onoff_phad->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_phad->DrawCopy("");

  h_phad_allsel[2]->SetLineColor(kRed);
  h_phad_allsel[2]->SetLineWidth(2);
  h_phad_allsel[2]->SetLineStyle(1);
  h_phad_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_phad = new THStack("hs_phad","");

  h_phad_sig[0]-> SetFillColor(2); 
  h_phad_sig[0]->Scale(normfac);
  h_phad_bac[0]-> SetFillColor(46);
  h_phad_bac[0]->Scale(normfac);
  h_phad_bac[1]-> SetFillColor(4);
  h_phad_bac[1]->Scale(normfac);
  h_phad_bac[2]-> SetFillColor(8);
  h_phad_bac[2]->Scale(normfac);
  h_phad_bac[3] -> SetFillColor(5);
  h_phad_bac[3]->Scale(normfac);
  h_phad_bac[4]-> SetFillColor(3);
  h_phad_bac[4]->Scale(normfac);
  h_phad_bac[5]-> SetFillColor(6);
  h_phad_bac[5]->Scale(normfac);
  h_phad_bac[6]-> SetFillColor(7);
  h_phad_bac[6]->Scale(normfac);
  h_phad_bac[7]-> SetFillColor(9);
  h_phad_bac[7]->Scale(normfac);
  h_phad_allsel[1]->SetFillStyle(3005);
  h_phad_allsel[1]->SetFillColor(28);
  h_phad_allsel[1]->Scale(scale_onoffbeam);



   hs_phad -> Add(h_phad_sig[0]);
   hs_phad -> Add(h_phad_bac[0]);
   hs_phad -> Add(h_phad_bac[1]);
   hs_phad -> Add(h_phad_bac[2]);
   hs_phad -> Add(h_phad_bac[3]);
   hs_phad -> Add(h_phad_bac[4]);
   hs_phad -> Add(h_phad_bac[5]);
   hs_phad -> Add(h_phad_bac[6]);
   hs_phad -> Add(h_phad_bac[7]);
   hs_phad -> Add(h_phad_allsel[1]);
   hs_phad -> Draw("HIST,SAME");

   h_phad_allsel[2]->Draw("same");
  h_phad_allsel[0]->Draw("same");  
   //h_onoff_phad->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_phad_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_phad_allsel.png");}
   
//==============================================================

  TH1D *h_onoff_thetamup=(TH1D*)h_thetamup_allsel[1]->Clone(Form("%s_on-off", h_thetamup_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_thetamup_allsel[0]->SetLineColor(kBlack);
  h_thetamup_allsel[0]->SetLineWidth(2);
  h_thetamup_allsel[0]->SetLineStyle(1);
  h_thetamup_allsel[0]->GetXaxis()->SetTitle("#theta_{#mu p}[Rad]");
  h_thetamup_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_thetamup_allsel[0]->SetMaximum(250);
  h_thetamup_allsel[0]->Draw();
  
  //h_thetamup_allsel[1]->SetLineColor(kRed);
  //h_thetamup_allsel[1]->SetLineWidth(2);
  //h_thetamup_allsel[1]->SetLineStyle(1);

  h_onoff_thetamup->Add(h_thetamup_allsel[0],h_thetamup_allsel[1],1,-scale_onoffbeam);
  h_onoff_thetamup->SetLineColor(kBlack);
  h_onoff_thetamup->SetLineWidth(2);
  h_onoff_thetamup->SetLineStyle(1);
  h_onoff_thetamup->GetXaxis()->SetTitle("|#theta_{#mu}-#theta_{P}|[Rad]");
  h_onoff_thetamup->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_thetamup->DrawCopy("");

  h_thetamup_allsel[2]->SetLineColor(kRed);
  h_thetamup_allsel[2]->SetLineWidth(2);
  h_thetamup_allsel[2]->SetLineStyle(1);
  h_thetamup_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_thetamup = new THStack("hs_thetamup","");

  h_thetamup_sig[0]-> SetFillColor(2); 
  h_thetamup_sig[0]->Scale(normfac);
  h_thetamup_bac[0]-> SetFillColor(46);
  h_thetamup_bac[0]->Scale(normfac);
  h_thetamup_bac[1]-> SetFillColor(4);
  h_thetamup_bac[1]->Scale(normfac);
  h_thetamup_bac[2]-> SetFillColor(8);
  h_thetamup_bac[2]->Scale(normfac);
  h_thetamup_bac[3] -> SetFillColor(5);
  h_thetamup_bac[3]->Scale(normfac);
  h_thetamup_bac[4]-> SetFillColor(3);
  h_thetamup_bac[4]->Scale(normfac);
  h_thetamup_bac[5]-> SetFillColor(6);
  h_thetamup_bac[5]->Scale(normfac);
  h_thetamup_bac[6]-> SetFillColor(7);
  h_thetamup_bac[6]->Scale(normfac);
  h_thetamup_bac[7]-> SetFillColor(9);
  h_thetamup_bac[7]->Scale(normfac);
  h_thetamup_allsel[1]->SetFillStyle(3005);
  h_thetamup_allsel[1]->SetFillColor(28);
  h_thetamup_allsel[1]->Scale(scale_onoffbeam);



   hs_thetamup -> Add(h_thetamup_sig[0]);
   hs_thetamup -> Add(h_thetamup_bac[0]);
   hs_thetamup -> Add(h_thetamup_bac[1]);
   hs_thetamup -> Add(h_thetamup_bac[2]);
   hs_thetamup -> Add(h_thetamup_bac[3]);
   hs_thetamup -> Add(h_thetamup_bac[4]);
   hs_thetamup -> Add(h_thetamup_bac[5]);
   hs_thetamup -> Add(h_thetamup_bac[6]);
   hs_thetamup -> Add(h_thetamup_bac[7]);
   hs_thetamup -> Add(h_thetamup_allsel[1]);
   hs_thetamup -> Draw("HIST,SAME");

   h_thetamup_allsel[2]->Draw("same");
   h_thetamup_allsel[0]->Draw("same");
   //h_onoff_thetamup->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_thetamup_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_thetamup_allsel.png");}
//==================================================================
  TH1D *h_onoff_thetapp=(TH1D*)h_thetapp_allsel[1]->Clone(Form("%s_on-off", h_thetapp_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_thetapp_allsel[0]->SetLineColor(kBlack);
  h_thetapp_allsel[0]->SetLineWidth(2);
  h_thetapp_allsel[0]->SetLineStyle(1);
  h_thetapp_allsel[0]->GetXaxis()->SetTitle("|#theta_{P1}-#theta_{P2}|[Rad]");
  h_thetapp_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_thetapp_allsel[0]->SetMaximum(30);
  h_thetapp_allsel[0]->Draw();

  //h_thetapp_allsel[1]->SetLineColor(kRed);
  //h_thetapp_allsel[1]->SetLineWidth(2);
  //h_thetapp_allsel[1]->SetLineStyle(1);

  h_onoff_thetapp->Add(h_thetapp_allsel[0],h_thetapp_allsel[1],1,-scale_onoffbeam);
  h_onoff_thetapp->SetLineColor(kBlack);
  h_onoff_thetapp->SetLineWidth(2);
  h_onoff_thetapp->SetLineStyle(1);
  h_onoff_thetapp->GetXaxis()->SetTitle("|#theta_{P1}-#theta_{P2}|[Rad]");
  h_onoff_thetapp->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_thetapp->DrawCopy("");

  h_thetapp_allsel[2]->SetLineColor(kRed);
  h_thetapp_allsel[2]->SetLineWidth(2);
  h_thetapp_allsel[2]->SetLineStyle(1);
  h_thetapp_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_thetapp = new THStack("hs_thetapp","");

  h_thetapp_sig[0]-> SetFillColor(2); 
  h_thetapp_sig[0]->Scale(normfac);
  h_thetapp_bac[0]-> SetFillColor(46);
  h_thetapp_bac[0]->Scale(normfac);
  h_thetapp_bac[1]-> SetFillColor(4);
  h_thetapp_bac[1]->Scale(normfac);
  h_thetapp_bac[2]-> SetFillColor(8);
  h_thetapp_bac[2]->Scale(normfac);
  h_thetapp_bac[3] -> SetFillColor(5);
  h_thetapp_bac[3]->Scale(normfac);
  h_thetapp_bac[4]-> SetFillColor(3);
  h_thetapp_bac[4]->Scale(normfac);
  h_thetapp_bac[5]-> SetFillColor(6);
  h_thetapp_bac[5]->Scale(normfac);
  h_thetapp_bac[6]-> SetFillColor(7);
  h_thetapp_bac[6]->Scale(normfac);
  h_thetapp_bac[7]-> SetFillColor(9);
  h_thetapp_bac[7]->Scale(normfac);
  h_thetapp_allsel[1]->SetFillStyle(3005);
  h_thetapp_allsel[1]->SetFillColor(28);
  h_thetapp_allsel[1]->Scale(scale_onoffbeam);



   hs_thetapp -> Add(h_thetapp_sig[0]);
   hs_thetapp -> Add(h_thetapp_bac[0]);
   hs_thetapp -> Add(h_thetapp_bac[1]);
   hs_thetapp -> Add(h_thetapp_bac[2]);
   hs_thetapp -> Add(h_thetapp_bac[3]);
   hs_thetapp -> Add(h_thetapp_bac[4]);
   hs_thetapp -> Add(h_thetapp_bac[5]);
   hs_thetapp -> Add(h_thetapp_bac[6]);
   hs_thetapp -> Add(h_thetapp_bac[7]);
   hs_thetapp -> Add(h_thetapp_allsel[1]);
   hs_thetapp -> Draw("HIST,SAME");

   h_thetapp_allsel[2]->Draw("same");
   h_thetapp_allsel[0]->Draw("same");
   //h_onoff_thetapp->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_thetapp_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_thetapp_allsel.png");}


//==================================================================

  TH1D *h_onoff_phimup=(TH1D*)h_phimup_allsel[1]->Clone(Form("%s_on-off", h_phimup_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_phimup_allsel[0]->SetLineColor(kBlack);
  h_phimup_allsel[0]->SetLineWidth(2);
  h_phimup_allsel[0]->SetLineStyle(1);
  h_phimup_allsel[0]->GetXaxis()->SetTitle("|#phi_{#mu}-#phi_{P}|[Rad]");
  h_phimup_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_phimup_allsel[0]->SetMaximum(600);
  h_phimup_allsel[0]->Draw();

  //h_phimup_allsel[1]->SetLineColor(kRed);
  //h_phimup_allsel[1]->SetLineWidth(2);
  //h_phimup_allsel[1]->SetLineStyle(1);

  h_onoff_phimup->Add(h_phimup_allsel[0],h_phimup_allsel[1],1,-scale_onoffbeam);
  h_onoff_phimup->SetLineColor(kBlack);
  h_onoff_phimup->SetLineWidth(2);
  h_onoff_phimup->SetLineStyle(1);
  h_onoff_phimup->SetMaximum(300);
  h_onoff_phimup->GetXaxis()->SetTitle("|#phi_{#mu}-#phi_{P}|[Rad]");
  h_onoff_phimup->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_phimup->DrawCopy("");

  h_phimup_allsel[2]->SetLineColor(kRed);
  h_phimup_allsel[2]->SetLineWidth(2);
  h_phimup_allsel[2]->SetLineStyle(1);
  h_phimup_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_phimup = new THStack("hs_phimup","");

  h_phimup_sig[0]-> SetFillColor(2); 
  h_phimup_sig[0]->Scale(normfac);
  h_phimup_bac[0]-> SetFillColor(46);
  h_phimup_bac[0]->Scale(normfac);
  h_phimup_bac[1]-> SetFillColor(4);
  h_phimup_bac[1]->Scale(normfac);
  h_phimup_bac[2]-> SetFillColor(8);
  h_phimup_bac[2]->Scale(normfac);
  h_phimup_bac[3] -> SetFillColor(5);
  h_phimup_bac[3]->Scale(normfac);
  h_phimup_bac[4]-> SetFillColor(3);
  h_phimup_bac[4]->Scale(normfac);
  h_phimup_bac[5]-> SetFillColor(6);
  h_phimup_bac[5]->Scale(normfac);
  h_phimup_bac[6]-> SetFillColor(7);
  h_phimup_bac[6]->Scale(normfac);
  h_phimup_bac[7]-> SetFillColor(9);
  h_phimup_bac[7]->Scale(normfac);
  h_phimup_allsel[1]->SetFillStyle(3005);
  h_phimup_allsel[1]->SetFillColor(28);
  h_phimup_allsel[1]->Scale(scale_onoffbeam);



   hs_phimup -> Add(h_phimup_sig[0]);
   hs_phimup -> Add(h_phimup_bac[0]);
   hs_phimup -> Add(h_phimup_bac[1]);
   hs_phimup -> Add(h_phimup_bac[2]);
   hs_phimup -> Add(h_phimup_bac[3]);
   hs_phimup -> Add(h_phimup_bac[4]);
   hs_phimup -> Add(h_phimup_bac[5]);
   hs_phimup -> Add(h_phimup_bac[6]);
   hs_phimup -> Add(h_phimup_bac[7]);
   hs_phimup -> Add(h_phimup_allsel[1]);
   hs_phimup -> Draw("HIST,SAME");

   h_phimup_allsel[2]->Draw("same");
  h_phimup_allsel[0]->Draw("same");
   //h_onoff_phimup->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_phimup_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_phimup_allsel.png");}
//================================================================
  TH1D *h_onoff_phipp=(TH1D*)h_phipp_allsel[1]->Clone(Form("%s_on-off", h_phipp_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_phipp_allsel[0]->SetLineColor(kBlack);
  h_phipp_allsel[0]->SetLineWidth(2);
  h_phipp_allsel[0]->SetLineStyle(1);
  h_phipp_allsel[0]->GetXaxis()->SetTitle("|#phi_{P1}-#phi_{P2}|[Rad]");
  h_phipp_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_phipp_allsel[0]->SetMaximum(25);
  h_phipp_allsel[0]->Draw();

  //h_phipp_allsel[1]->SetLineColor(kRed);
  //h_phipp_allsel[1]->SetLineWidth(2);
  //h_phipp_allsel[1]->SetLineStyle(1);

  h_onoff_phipp->Add(h_phipp_allsel[0],h_phipp_allsel[1],1,-scale_onoffbeam);
  h_onoff_phipp->SetLineColor(kBlack);
  h_onoff_phipp->SetLineWidth(2);
  h_onoff_phipp->SetLineStyle(1);
  h_onoff_phipp->SetMaximum(12);
  h_onoff_phipp->GetXaxis()->SetTitle("|#phi_{P1}-#phi_{P2}|[Rad]");
  h_onoff_phipp->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_phipp->DrawCopy("");

  h_phipp_allsel[2]->SetLineColor(kRed);
  h_phipp_allsel[2]->SetLineWidth(2);
  h_phipp_allsel[2]->SetLineStyle(1);
  h_phipp_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_phipp = new THStack("hs_phipp","");

  h_phipp_sig[0]-> SetFillColor(2); 
  h_phipp_sig[0]->Scale(normfac);
  h_phipp_bac[0]-> SetFillColor(46);
  h_phipp_bac[0]->Scale(normfac);
  h_phipp_bac[1]-> SetFillColor(4);
  h_phipp_bac[1]->Scale(normfac);
  h_phipp_bac[2]-> SetFillColor(8);
  h_phipp_bac[2]->Scale(normfac);
  h_phipp_bac[3] -> SetFillColor(5);
  h_phipp_bac[3]->Scale(normfac);
  h_phipp_bac[4]-> SetFillColor(3);
  h_phipp_bac[4]->Scale(normfac);
  h_phipp_bac[5]-> SetFillColor(6);
  h_phipp_bac[5]->Scale(normfac);
  h_phipp_bac[6]-> SetFillColor(7);
  h_phipp_bac[6]->Scale(normfac);
  h_phipp_bac[7]-> SetFillColor(9);
  h_phipp_bac[7]->Scale(normfac);
  h_phipp_allsel[1]->SetFillStyle(3005);
  h_phipp_allsel[1]->SetFillColor(28);
  h_phipp_allsel[1]->Scale(scale_onoffbeam);



   hs_phipp -> Add(h_phipp_sig[0]);
   hs_phipp -> Add(h_phipp_bac[0]);
   hs_phipp -> Add(h_phipp_bac[1]);
   hs_phipp -> Add(h_phipp_bac[2]);
   hs_phipp -> Add(h_phipp_bac[3]);
   hs_phipp -> Add(h_phipp_bac[4]);
   hs_phipp -> Add(h_phipp_bac[5]);
   hs_phipp -> Add(h_phipp_bac[6]);
   hs_phipp -> Add(h_phipp_bac[7]);
   hs_phipp -> Add(h_phipp_allsel[1]);
   hs_phipp -> Draw("HIST,SAME");

   h_phipp_allsel[2]->Draw("same");
  h_phipp_allsel[0]->Draw("same");
   //h_onoff_phipp->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_phipp_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_phipp_allsel.png");}


//--------------------------------------------------------------

  TH1D *h_onoff_Nhitsmuon=(TH1D*)h_Nhitsmuon_allsel[1]->Clone(Form("%s_on-off", h_Nhitsmuon_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_Nhitsmuon_allsel[0]->SetLineColor(kBlack);
  h_Nhitsmuon_allsel[0]->SetLineWidth(2);
  h_Nhitsmuon_allsel[0]->SetLineStyle(1);
  h_Nhitsmuon_allsel[0]->GetXaxis()->SetTitle("Number of Hits of Muon Candidate");
  h_Nhitsmuon_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_Nhitsmuon_allsel[0]->SetMaximum(500);
  h_Nhitsmuon_allsel[0]->Draw();

  //h_Nhitsmuon_allsel[1]->SetLineColor(kRed);
  //h_Nhitsmuon_allsel[1]->SetLineWidth(2);
  //h_Nhitsmuon_allsel[1]->SetLineStyle(1);

  h_onoff_Nhitsmuon->Add(h_Nhitsmuon_allsel[0],h_Nhitsmuon_allsel[1],1,-scale_onoffbeam);
  h_onoff_Nhitsmuon->SetLineColor(kBlack);
  h_onoff_Nhitsmuon->SetLineWidth(2);
  h_onoff_Nhitsmuon->SetLineStyle(1);
  h_onoff_Nhitsmuon->GetXaxis()->SetTitle("Number of Hits of Muon Candidate");
  h_onoff_Nhitsmuon->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_Nhitsmuon->DrawCopy("");

  h_Nhitsmuon_allsel[2]->SetLineColor(kRed);
  h_Nhitsmuon_allsel[2]->SetLineWidth(2);
  h_Nhitsmuon_allsel[2]->SetLineStyle(1);
  h_Nhitsmuon_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_Nhitsmuon = new THStack("hs_Nhitsmuon","");

  h_Nhitsmuon_sig[0]-> SetFillColor(2); 
  h_Nhitsmuon_sig[0]->Scale(normfac);
  h_Nhitsmuon_bac[0]-> SetFillColor(46);
  h_Nhitsmuon_bac[0]->Scale(normfac);
  h_Nhitsmuon_bac[1]-> SetFillColor(4);
  h_Nhitsmuon_bac[1]->Scale(normfac);
  h_Nhitsmuon_bac[2]-> SetFillColor(8);
  h_Nhitsmuon_bac[2]->Scale(normfac);
  h_Nhitsmuon_bac[3] -> SetFillColor(5);
  h_Nhitsmuon_bac[3]->Scale(normfac);
  h_Nhitsmuon_bac[4]-> SetFillColor(3);
  h_Nhitsmuon_bac[4]->Scale(normfac);
  h_Nhitsmuon_bac[5]-> SetFillColor(6);
  h_Nhitsmuon_bac[5]->Scale(normfac);
  h_Nhitsmuon_bac[6]-> SetFillColor(7);
  h_Nhitsmuon_bac[6]->Scale(normfac);
  h_Nhitsmuon_bac[7]-> SetFillColor(9);
  h_Nhitsmuon_bac[7]->Scale(normfac);
  h_Nhitsmuon_allsel[1]->SetFillStyle(3005);
  h_Nhitsmuon_allsel[1]->SetFillColor(28);
  h_Nhitsmuon_allsel[1]->Scale(scale_onoffbeam);



   hs_Nhitsmuon -> Add(h_Nhitsmuon_sig[0]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_bac[0]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_bac[1]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_bac[2]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_bac[3]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_bac[4]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_bac[5]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_bac[6]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_bac[7]);
   hs_Nhitsmuon -> Add(h_Nhitsmuon_allsel[1]);
   hs_Nhitsmuon -> Draw("HIST,SAME");

   h_Nhitsmuon_allsel[2]->Draw("same");
   h_Nhitsmuon_allsel[0]->Draw("same");
   //h_onoff_Nhitsmuon->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_Nhitsmuon_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_Nhitsmuon_allsel.png");}

//-----------------------------------------------------------

  TH1D *h_onoff_Nhitsproton=(TH1D*)h_Nhitsproton_allsel[1]->Clone(Form("%s_on-off", h_Nhitsproton_allsel[1]->GetName()));
  cout<<"get all the histograms!"<<endl;


  h_Nhitsproton_allsel[0]->SetLineColor(kBlack);
  h_Nhitsproton_allsel[0]->SetLineWidth(2);
  h_Nhitsproton_allsel[0]->SetLineStyle(1);
  h_Nhitsproton_allsel[0]->GetXaxis()->SetTitle("Number of Hits of The Proton Candidate");
  h_Nhitsproton_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_Nhitsproton_allsel[0]->SetMaximum(600);
  h_Nhitsproton_allsel[0]->Draw();

  //h_Nhitsproton_allsel[1]->SetLineColor(kRed);
  //h_Nhitsproton_allsel[1]->SetLineWidth(2);
  //h_Nhitsproton_allsel[1]->SetLineStyle(1);

  h_onoff_Nhitsproton->Add(h_Nhitsproton_allsel[0],h_Nhitsproton_allsel[1],1,-scale_onoffbeam);
  h_onoff_Nhitsproton->SetLineColor(kBlack);
  h_onoff_Nhitsproton->SetLineWidth(2);
  h_onoff_Nhitsproton->SetLineStyle(1);
  h_onoff_Nhitsproton->GetXaxis()->SetTitle("Number of Hits of The Proton Candidate");
  h_onoff_Nhitsproton->GetYaxis()->SetTitle("No. of Tracks");
  //h_onoff_Nhitsproton->DrawCopy("");

  h_Nhitsproton_allsel[2]->SetLineColor(kRed);
  h_Nhitsproton_allsel[2]->SetLineWidth(2);
  h_Nhitsproton_allsel[2]->SetLineStyle(1);
  h_Nhitsproton_allsel[2]->Scale(normfac);
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_Nhitsproton = new THStack("hs_Nhitsproton","");

  h_Nhitsproton_sig[0]-> SetFillColor(2); 
  h_Nhitsproton_sig[0]->Scale(normfac);
  h_Nhitsproton_bac[0]-> SetFillColor(46);
  h_Nhitsproton_bac[0]->Scale(normfac);
  h_Nhitsproton_bac[1]-> SetFillColor(4);
  h_Nhitsproton_bac[1]->Scale(normfac);
  h_Nhitsproton_bac[2]-> SetFillColor(8);
  h_Nhitsproton_bac[2]->Scale(normfac);
  h_Nhitsproton_bac[3] -> SetFillColor(5);
  h_Nhitsproton_bac[3]->Scale(normfac);
  h_Nhitsproton_bac[4]-> SetFillColor(3);
  h_Nhitsproton_bac[4]->Scale(normfac);
  h_Nhitsproton_bac[5]-> SetFillColor(6);
  h_Nhitsproton_bac[5]->Scale(normfac);
  h_Nhitsproton_bac[6]-> SetFillColor(7);
  h_Nhitsproton_bac[6]->Scale(normfac);
  h_Nhitsproton_bac[7]-> SetFillColor(9);
  h_Nhitsproton_bac[7]->Scale(normfac);
  h_Nhitsproton_allsel[1]->SetFillStyle(3005);
  h_Nhitsproton_allsel[1]->SetFillColor(28);
  h_Nhitsproton_allsel[1]->Scale(scale_onoffbeam);



   hs_Nhitsproton -> Add(h_Nhitsproton_sig[0]);
   hs_Nhitsproton -> Add(h_Nhitsproton_bac[0]);
   hs_Nhitsproton -> Add(h_Nhitsproton_bac[1]);
   hs_Nhitsproton -> Add(h_Nhitsproton_bac[2]);
   hs_Nhitsproton -> Add(h_Nhitsproton_bac[3]);
   hs_Nhitsproton -> Add(h_Nhitsproton_bac[4]);
   hs_Nhitsproton -> Add(h_Nhitsproton_bac[5]);
   hs_Nhitsproton -> Add(h_Nhitsproton_bac[6]);
   hs_Nhitsproton -> Add(h_Nhitsproton_bac[7]);
   hs_Nhitsproton -> Add(h_Nhitsproton_allsel[1]);
   hs_Nhitsproton -> Draw("HIST,SAME");

   h_Nhitsproton_allsel[2]->Draw("same");
   h_Nhitsproton_allsel[0]->Draw("same");
   //h_onoff_Nhitsproton->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_Nhitsproton_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_Nhitsproton_allsel.png");}
   //=====================================================================
   TH1D *h_onoff_vtxx=(TH1D*)h_vtxx_allsel[1]->Clone(Form("%s_on-off", h_vtxx_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_vtxx_allsel[0]->SetLineColor(kBlack);
   h_vtxx_allsel[0]->SetLineWidth(2);
   h_vtxx_allsel[0]->SetLineStyle(1);
   h_vtxx_allsel[0]->GetXaxis()->SetTitle("Vertex Position X[cm]");
   h_vtxx_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_vtxx_allsel[0]->SetMaximum(200); 
   h_vtxx_allsel[0]->Draw();

   //h_vtxx_allsel[1]->SetLineColor(kRed);
   //h_vtxx_allsel[1]->SetLineWidth(2);
   //h_vtxx_allsel[1]->SetLineStyle(1);

   h_onoff_vtxx->Add(h_vtxx_allsel[0],h_vtxx_allsel[1],1,-scale_onoffbeam);
   h_onoff_vtxx->SetLineColor(kBlack);
   h_onoff_vtxx->SetLineWidth(2);
   h_onoff_vtxx->SetLineStyle(1);
   h_onoff_vtxx->SetMaximum(120);
   h_onoff_vtxx->GetXaxis()->SetTitle("Vertex Position X[cm]");
   h_onoff_vtxx->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_vtxx->DrawCopy("");

   h_vtxx_allsel[2]->SetLineColor(kRed);
   h_vtxx_allsel[2]->SetLineWidth(2);
   h_vtxx_allsel[2]->SetLineStyle(1);
   h_vtxx_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_vtxx = new THStack("hs_vtxx","");

   h_vtxx_sig[0]-> SetFillColor(2); 
   h_vtxx_sig[0]->Scale(normfac);
   h_vtxx_bac[0]-> SetFillColor(46);
   h_vtxx_bac[0]->Scale(normfac);
   h_vtxx_bac[1]-> SetFillColor(4);
   h_vtxx_bac[1]->Scale(normfac);
   h_vtxx_bac[2]-> SetFillColor(8);
   h_vtxx_bac[2]->Scale(normfac);
   h_vtxx_bac[3] -> SetFillColor(5);
   h_vtxx_bac[3]->Scale(normfac);
   h_vtxx_bac[4]-> SetFillColor(3);
   h_vtxx_bac[4]->Scale(normfac);
   h_vtxx_bac[5]-> SetFillColor(6);
   h_vtxx_bac[5]->Scale(normfac);
   h_vtxx_bac[6]-> SetFillColor(7);
   h_vtxx_bac[6]->Scale(normfac);
   h_vtxx_bac[7]-> SetFillColor(9);
   h_vtxx_bac[7]->Scale(normfac);
  h_vtxx_allsel[1]->SetFillStyle(3005);
  h_vtxx_allsel[1]->SetFillColor(28);
  h_vtxx_allsel[1]->Scale(scale_onoffbeam);



   hs_vtxx -> Add(h_vtxx_sig[0]);
   hs_vtxx -> Add(h_vtxx_bac[0]);
   hs_vtxx -> Add(h_vtxx_bac[1]);
   hs_vtxx -> Add(h_vtxx_bac[2]);
   hs_vtxx -> Add(h_vtxx_bac[3]);
   hs_vtxx -> Add(h_vtxx_bac[4]);
   hs_vtxx -> Add(h_vtxx_bac[5]);
   hs_vtxx -> Add(h_vtxx_bac[6]);
   hs_vtxx -> Add(h_vtxx_bac[7]);
   hs_vtxx -> Add(h_vtxx_allsel[1]);
   hs_vtxx -> Draw("HIST,SAME");

   h_vtxx_allsel[2]->Draw("same");
   h_vtxx_allsel[0]->Draw("same");
   //h_onoff_vtxx->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_vtxx_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_vtxx_allsel.png");}
   //============================================================
   TH1D *h_onoff_vtxy=(TH1D*)h_vtxy_allsel[1]->Clone(Form("%s_on-off", h_vtxy_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_vtxy_allsel[0]->SetLineColor(kBlack);
   h_vtxy_allsel[0]->SetLineWidth(2);
   h_vtxy_allsel[0]->SetLineStyle(1);
   h_vtxy_allsel[0]->GetXaxis()->SetTitle("Vertex Position Y[cm]");
   h_vtxy_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_vtxy_allsel[0]->SetMaximum(200);
   h_vtxy_allsel[0]->Draw();

   //h_vtxy_allsel[1]->SetLineColor(kRed);
   //h_vtxy_allsel[1]->SetLineWidth(2);
   //h_vtxy_allsel[1]->SetLineStyle(1);

   h_onoff_vtxy->Add(h_vtxy_allsel[0],h_vtxy_allsel[1],1,-scale_onoffbeam);
   h_onoff_vtxy->SetLineColor(kBlack);
   h_onoff_vtxy->SetLineWidth(2);
   h_onoff_vtxy->SetLineStyle(1);
   h_onoff_vtxy->SetMaximum(120);
   h_onoff_vtxy->GetXaxis()->SetTitle("Vertex Position Y[cm]");
   h_onoff_vtxy->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_vtxy->DrawCopy("");

   h_vtxy_allsel[2]->SetLineColor(kRed);
   h_vtxy_allsel[2]->SetLineWidth(2);
   h_vtxy_allsel[2]->SetLineStyle(1);
   h_vtxy_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_vtxy = new THStack("hs_vtxy","");

   h_vtxy_sig[0]-> SetFillColor(2); 
   h_vtxy_sig[0]->Scale(normfac);
   h_vtxy_bac[0]-> SetFillColor(46);
   h_vtxy_bac[0]->Scale(normfac);
   h_vtxy_bac[1]-> SetFillColor(4);
   h_vtxy_bac[1]->Scale(normfac);
   h_vtxy_bac[2]-> SetFillColor(8);
   h_vtxy_bac[2]->Scale(normfac);
   h_vtxy_bac[3] -> SetFillColor(5);
   h_vtxy_bac[3]->Scale(normfac);
   h_vtxy_bac[4]-> SetFillColor(3);
   h_vtxy_bac[4]->Scale(normfac);
   h_vtxy_bac[5]-> SetFillColor(6);
   h_vtxy_bac[5]->Scale(normfac);
   h_vtxy_bac[6]-> SetFillColor(7);
   h_vtxy_bac[6]->Scale(normfac);
   h_vtxy_bac[7]-> SetFillColor(9);
   h_vtxy_bac[7]->Scale(normfac);
  h_vtxy_allsel[1]->SetFillStyle(3005);
  h_vtxy_allsel[1]->SetFillColor(28);
  h_vtxy_allsel[1]->Scale(scale_onoffbeam);



   hs_vtxy -> Add(h_vtxy_sig[0]);
   hs_vtxy -> Add(h_vtxy_bac[0]);
   hs_vtxy -> Add(h_vtxy_bac[1]);
   hs_vtxy -> Add(h_vtxy_bac[2]);
   hs_vtxy -> Add(h_vtxy_bac[3]);
   hs_vtxy -> Add(h_vtxy_bac[4]);
   hs_vtxy -> Add(h_vtxy_bac[5]);
   hs_vtxy -> Add(h_vtxy_bac[6]);
   hs_vtxy -> Add(h_vtxy_bac[7]);
   hs_vtxy -> Add(h_vtxy_allsel[1]);
   hs_vtxy -> Draw("HIST,SAME");

   h_vtxy_allsel[2]->Draw("same");
   h_vtxy_allsel[0]->Draw("same");
   //h_onoff_vtxy->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_vtxy_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_vtxy_allsel.png");}
   //================================================================================ 
   TH1D *h_onoff_vtxz=(TH1D*)h_vtxz_allsel[1]->Clone(Form("%s_on-off", h_vtxz_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_vtxz_allsel[0]->SetLineColor(kBlack);
   h_vtxz_allsel[0]->SetLineWidth(2);
   h_vtxz_allsel[0]->SetLineStyle(1);
   h_vtxz_allsel[0]->GetXaxis()->SetTitle("Vertex Position Z[cm]");
   h_vtxz_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_vtxz_allsel[0]->SetMaximum(200);
   h_vtxz_allsel[0]->Draw();
  
   //h_vtxz_allsel[1]->SetLineColor(kRed);
   //h_vtxz_allsel[1]->SetLineWidth(2);
   //h_vtxz_allsel[1]->SetLineStyle(1);

   h_onoff_vtxz->Add(h_vtxz_allsel[0],h_vtxz_allsel[1],1,-scale_onoffbeam);
   h_onoff_vtxz->SetLineColor(kBlack);
   h_onoff_vtxz->SetLineWidth(2);
   h_onoff_vtxz->SetLineStyle(1);
   h_onoff_vtxz->SetMaximum(120);
   h_onoff_vtxz->GetXaxis()->SetTitle("Vertex Position Z[cm]");
   h_onoff_vtxz->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_vtxz->DrawCopy("");

   h_vtxz_allsel[2]->SetLineColor(kRed);
   h_vtxz_allsel[2]->SetLineWidth(2);
   h_vtxz_allsel[2]->SetLineStyle(1);
   h_vtxz_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_vtxz = new THStack("hs_vtxz","");

   h_vtxz_sig[0]-> SetFillColor(2); 
   h_vtxz_sig[0]->Scale(normfac);
   h_vtxz_bac[0]-> SetFillColor(46);
   h_vtxz_bac[0]->Scale(normfac);
   h_vtxz_bac[1]-> SetFillColor(4);
   h_vtxz_bac[1]->Scale(normfac);
   h_vtxz_bac[2]-> SetFillColor(8);
   h_vtxz_bac[2]->Scale(normfac);
   h_vtxz_bac[3] -> SetFillColor(5);
   h_vtxz_bac[3]->Scale(normfac);
   h_vtxz_bac[4]-> SetFillColor(3);
   h_vtxz_bac[4]->Scale(normfac);
   h_vtxz_bac[5]-> SetFillColor(6);
   h_vtxz_bac[5]->Scale(normfac);
   h_vtxz_bac[6]-> SetFillColor(7);
   h_vtxz_bac[6]->Scale(normfac);
   h_vtxz_bac[7]-> SetFillColor(9);
   h_vtxz_bac[7]->Scale(normfac);
  h_vtxz_allsel[1]->SetFillStyle(3005);
  h_vtxz_allsel[1]->SetFillColor(28);
  h_vtxz_allsel[1]->Scale(scale_onoffbeam);



   hs_vtxz -> Add(h_vtxz_sig[0]);
   hs_vtxz -> Add(h_vtxz_bac[0]);
   hs_vtxz -> Add(h_vtxz_bac[1]);
   hs_vtxz -> Add(h_vtxz_bac[2]);
   hs_vtxz -> Add(h_vtxz_bac[3]);
   hs_vtxz -> Add(h_vtxz_bac[4]);
   hs_vtxz -> Add(h_vtxz_bac[5]);
   hs_vtxz -> Add(h_vtxz_bac[6]);
   hs_vtxz -> Add(h_vtxz_bac[7]);
   hs_vtxz -> Add(h_vtxz_allsel[1]);
   hs_vtxz -> Draw("HIST,SAME");

   h_vtxz_allsel[2]->Draw("same");
   h_vtxz_allsel[0]->Draw("same");
   //h_onoff_vtxz->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_vtxz_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_vtxz_allsel.png");}
   //==================================================    
   TH1D *h_onoff_mustartx=(TH1D*)h_mustartx_allsel[1]->Clone(Form("%s_on-off", h_mustartx_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_mustartx_allsel[0]->SetLineColor(kBlack);
   h_mustartx_allsel[0]->SetLineWidth(2);
   h_mustartx_allsel[0]->SetLineStyle(1);
   h_mustartx_allsel[0]->GetXaxis()->SetTitle("Start Position X of Muon Candidate[cm]");
   h_mustartx_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_mustartx_allsel[0]->SetMaximum(160);
   h_mustartx_allsel[0]->Draw();

   //h_mustartx_allsel[1]->SetLineColor(kRed);
   //h_mustartx_allsel[1]->SetLineWidth(2);
   //h_mustartx_allsel[1]->SetLineStyle(1);

   h_onoff_mustartx->Add(h_mustartx_allsel[0],h_mustartx_allsel[1],1,-scale_onoffbeam);
   h_onoff_mustartx->SetLineColor(kBlack);
   h_onoff_mustartx->SetLineWidth(2);
   h_onoff_mustartx->SetLineStyle(1);
   h_onoff_mustartx->SetMaximum(120);
   h_onoff_mustartx->GetXaxis()->SetTitle("Start Position X of Muon Candidate[cm]");
   h_onoff_mustartx->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_mustartx->DrawCopy("");

   h_mustartx_allsel[2]->SetLineColor(kRed);
   h_mustartx_allsel[2]->SetLineWidth(2);
   h_mustartx_allsel[2]->SetLineStyle(1);
   h_mustartx_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_mustartx = new THStack("hs_mustartx","");

   h_mustartx_sig[0]-> SetFillColor(2); 
   h_mustartx_sig[0]->Scale(normfac);
   h_mustartx_bac[0]-> SetFillColor(46);
   h_mustartx_bac[0]->Scale(normfac);
   h_mustartx_bac[1]-> SetFillColor(4);
   h_mustartx_bac[1]->Scale(normfac);
   h_mustartx_bac[2]-> SetFillColor(8);
   h_mustartx_bac[2]->Scale(normfac);
   h_mustartx_bac[3] -> SetFillColor(5);
   h_mustartx_bac[3]->Scale(normfac);
   h_mustartx_bac[4]-> SetFillColor(3);
   h_mustartx_bac[4]->Scale(normfac);
   h_mustartx_bac[5]-> SetFillColor(6);
   h_mustartx_bac[5]->Scale(normfac);
   h_mustartx_bac[6]-> SetFillColor(7);
   h_mustartx_bac[6]->Scale(normfac);
   h_mustartx_bac[7]-> SetFillColor(9);
   h_mustartx_bac[7]->Scale(normfac);
  h_mustartx_allsel[1]->SetFillStyle(3005);
  h_mustartx_allsel[1]->SetFillColor(28);
  h_mustartx_allsel[1]->Scale(scale_onoffbeam);



   hs_mustartx -> Add(h_mustartx_sig[0]);
   hs_mustartx -> Add(h_mustartx_bac[0]);
   hs_mustartx -> Add(h_mustartx_bac[1]);
   hs_mustartx -> Add(h_mustartx_bac[2]);
   hs_mustartx -> Add(h_mustartx_bac[3]);
   hs_mustartx -> Add(h_mustartx_bac[4]);
   hs_mustartx -> Add(h_mustartx_bac[5]);
   hs_mustartx -> Add(h_mustartx_bac[6]);
   hs_mustartx -> Add(h_mustartx_bac[7]);
   hs_mustartx -> Add(h_mustartx_allsel[1]);
   hs_mustartx -> Draw("HIST,SAME");

   h_mustartx_allsel[2]->Draw("same");
   h_mustartx_allsel[0]->Draw("same");
   //h_onoff_mustartx->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_mustartx_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_mustartx_allsel.png");}
 
   //==========================================================
   TH1D *h_onoff_mustarty=(TH1D*)h_mustarty_allsel[1]->Clone(Form("%s_on-off", h_mustarty_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_mustarty_allsel[0]->SetLineColor(kBlack);
   h_mustarty_allsel[0]->SetLineWidth(2);
   h_mustarty_allsel[0]->SetLineStyle(1);
   h_mustarty_allsel[0]->GetXaxis()->SetTitle("Start Position Y of Muon Candidate[cm]");
   h_mustarty_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_mustarty_allsel[0]->SetMaximum(200);
    h_mustarty_allsel[0]->Draw();
  
   //h_mustarty_allsel[1]->SetLineColor(kRed);
   //h_mustarty_allsel[1]->SetLineWidth(2);
   //h_mustarty_allsel[1]->SetLineStyle(1);

   h_onoff_mustarty->Add(h_mustarty_allsel[0],h_mustarty_allsel[1],1,-scale_onoffbeam);
   h_onoff_mustarty->SetLineColor(kBlack);
   h_onoff_mustarty->SetLineWidth(2);
   h_onoff_mustarty->SetLineStyle(1);
   h_onoff_mustarty->SetMaximum(120);
   h_onoff_mustarty->GetXaxis()->SetTitle("Start Position Y of Muon Candidate[cm]");
   h_onoff_mustarty->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_mustarty->DrawCopy("");

   h_mustarty_allsel[2]->SetLineColor(kRed);
   h_mustarty_allsel[2]->SetLineWidth(2);
   h_mustarty_allsel[2]->SetLineStyle(1);
   h_mustarty_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_mustarty = new THStack("hs_mustarty","");

   h_mustarty_sig[0]-> SetFillColor(2); 
   h_mustarty_sig[0]->Scale(normfac);
   h_mustarty_bac[0]-> SetFillColor(46);
   h_mustarty_bac[0]->Scale(normfac);
   h_mustarty_bac[1]-> SetFillColor(4);
   h_mustarty_bac[1]->Scale(normfac);
   h_mustarty_bac[2]-> SetFillColor(8);
   h_mustarty_bac[2]->Scale(normfac);
   h_mustarty_bac[3] -> SetFillColor(5);
   h_mustarty_bac[3]->Scale(normfac);
   h_mustarty_bac[4]-> SetFillColor(3);
   h_mustarty_bac[4]->Scale(normfac);
   h_mustarty_bac[5]-> SetFillColor(6);
   h_mustarty_bac[5]->Scale(normfac);
   h_mustarty_bac[6]-> SetFillColor(7);
   h_mustarty_bac[6]->Scale(normfac);
   h_mustarty_bac[7]-> SetFillColor(9);
   h_mustarty_bac[7]->Scale(normfac);
  h_mustarty_allsel[1]->SetFillStyle(3005);
  h_mustarty_allsel[1]->SetFillColor(28);
  h_mustarty_allsel[1]->Scale(scale_onoffbeam);



   hs_mustarty -> Add(h_mustarty_sig[0]);
   hs_mustarty -> Add(h_mustarty_bac[0]);
   hs_mustarty -> Add(h_mustarty_bac[1]);
   hs_mustarty -> Add(h_mustarty_bac[2]);
   hs_mustarty -> Add(h_mustarty_bac[3]);
   hs_mustarty -> Add(h_mustarty_bac[4]);
   hs_mustarty -> Add(h_mustarty_bac[5]);
   hs_mustarty -> Add(h_mustarty_bac[6]);
   hs_mustarty -> Add(h_mustarty_bac[7]);
   hs_mustarty -> Add(h_mustarty_allsel[1]);
   hs_mustarty -> Draw("HIST,SAME");

   h_mustarty_allsel[2]->Draw("same");
    h_mustarty_allsel[0]->Draw("same");
   //h_onoff_mustarty->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_mustarty_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_mustarty_allsel.png");}
  //==================================================================== 
   TH1D *h_onoff_mustartz=(TH1D*)h_mustartz_allsel[1]->Clone(Form("%s_on-off", h_mustartz_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_mustartz_allsel[0]->SetLineColor(kBlack);
   h_mustartz_allsel[0]->SetLineWidth(2);
   h_mustartz_allsel[0]->SetLineStyle(1);
   h_mustartz_allsel[0]->GetXaxis()->SetTitle("Start Position Z of Muon Candidate[cm]");
   h_mustartz_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_mustartz_allsel[0]->SetMaximum(200);
    h_mustartz_allsel[0]->Draw();

   //h_mustartz_allsel[1]->SetLineColor(kRed);
   //h_mustartz_allsel[1]->SetLineWidth(2);
   //h_mustartz_allsel[1]->SetLineStyle(1);

   h_onoff_mustartz->Add(h_mustartz_allsel[0],h_mustartz_allsel[1],1,-scale_onoffbeam);
   h_onoff_mustartz->SetLineColor(kBlack);
   h_onoff_mustartz->SetLineWidth(2);
   h_onoff_mustartz->SetLineStyle(1);
   h_onoff_mustartz->SetMaximum(120);
   h_onoff_mustartz->GetXaxis()->SetTitle("Start Position Z of Muon Candidate[cm]");
   h_onoff_mustartz->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_mustartz->DrawCopy("");

   h_mustartz_allsel[2]->SetLineColor(kRed);
   h_mustartz_allsel[2]->SetLineWidth(2);
   h_mustartz_allsel[2]->SetLineStyle(1);
   h_mustartz_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_mustartz = new THStack("hs_mustartz","");

   h_mustartz_sig[0]-> SetFillColor(2); 
   h_mustartz_sig[0]->Scale(normfac);
   h_mustartz_bac[0]-> SetFillColor(46);
   h_mustartz_bac[0]->Scale(normfac);
   h_mustartz_bac[1]-> SetFillColor(4);
   h_mustartz_bac[1]->Scale(normfac);
   h_mustartz_bac[2]-> SetFillColor(8);
   h_mustartz_bac[2]->Scale(normfac);
   h_mustartz_bac[3] -> SetFillColor(5);
   h_mustartz_bac[3]->Scale(normfac);
   h_mustartz_bac[4]-> SetFillColor(3);
   h_mustartz_bac[4]->Scale(normfac);
   h_mustartz_bac[5]-> SetFillColor(6);
   h_mustartz_bac[5]->Scale(normfac);
   h_mustartz_bac[6]-> SetFillColor(7);
   h_mustartz_bac[6]->Scale(normfac);
   h_mustartz_bac[7]-> SetFillColor(9);
   h_mustartz_bac[7]->Scale(normfac);
  h_mustartz_allsel[1]->SetFillStyle(3005);
  h_mustartz_allsel[1]->SetFillColor(28);
  h_mustartz_allsel[1]->Scale(scale_onoffbeam);



   hs_mustartz -> Add(h_mustartz_sig[0]);
   hs_mustartz -> Add(h_mustartz_bac[0]);
   hs_mustartz -> Add(h_mustartz_bac[1]);
   hs_mustartz -> Add(h_mustartz_bac[2]);
   hs_mustartz -> Add(h_mustartz_bac[3]);
   hs_mustartz -> Add(h_mustartz_bac[4]);
   hs_mustartz -> Add(h_mustartz_bac[5]);
   hs_mustartz -> Add(h_mustartz_bac[6]);
   hs_mustartz -> Add(h_mustartz_bac[7]);
   hs_mustartz -> Add(h_mustartz_allsel[1]);
   hs_mustartz -> Draw("HIST,SAME");

   h_mustartz_allsel[2]->Draw("same");
    h_mustartz_allsel[0]->Draw("same");
   //h_onoff_mustartz->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_mustartz_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_mustartz_allsel.png");}
  //================================================================= 
   TH1D *h_onoff_muendx=(TH1D*)h_muendx_allsel[1]->Clone(Form("%s_on-off", h_muendx_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_muendx_allsel[0]->SetLineColor(kBlack);
   h_muendx_allsel[0]->SetLineWidth(2);
   h_muendx_allsel[0]->SetLineStyle(1);
   h_muendx_allsel[0]->GetXaxis()->SetTitle("End Position X of Muon Candidate[cm]");
   h_muendx_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_muendx_allsel[0]->SetMaximum(500);
   h_muendx_allsel[0]->Draw();

   //h_muendx_allsel[1]->SetLineColor(kRed);
   //h_muendx_allsel[1]->SetLineWidth(2);
   //h_muendx_allsel[1]->SetLineStyle(1);

   h_onoff_muendx->Add(h_muendx_allsel[0],h_muendx_allsel[1],1,-scale_onoffbeam);
   h_onoff_muendx->SetLineColor(kBlack);
   h_onoff_muendx->SetLineWidth(2);
   h_onoff_muendx->SetLineStyle(1);
   h_onoff_muendx->SetMaximum(300);
   h_onoff_muendx->GetXaxis()->SetTitle("End Position X of Muon Candidate[cm]");
   h_onoff_muendx->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_muendx->DrawCopy("");

   h_muendx_allsel[2]->SetLineColor(kRed);
   h_muendx_allsel[2]->SetLineWidth(2);
   h_muendx_allsel[2]->SetLineStyle(1);
   h_muendx_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_muendx = new THStack("hs_muendx","");

   h_muendx_sig[0]-> SetFillColor(2); 
   h_muendx_sig[0]->Scale(normfac);
   h_muendx_bac[0]-> SetFillColor(46);
   h_muendx_bac[0]->Scale(normfac);
   h_muendx_bac[1]-> SetFillColor(4);
   h_muendx_bac[1]->Scale(normfac);
   h_muendx_bac[2]-> SetFillColor(8);
   h_muendx_bac[2]->Scale(normfac);
   h_muendx_bac[3] -> SetFillColor(5);
   h_muendx_bac[3]->Scale(normfac);
   h_muendx_bac[4]-> SetFillColor(3);
   h_muendx_bac[4]->Scale(normfac);
   h_muendx_bac[5]-> SetFillColor(6);
   h_muendx_bac[5]->Scale(normfac);
   h_muendx_bac[6]-> SetFillColor(7);
   h_muendx_bac[6]->Scale(normfac);
   h_muendx_bac[7]-> SetFillColor(9);
   h_muendx_bac[7]->Scale(normfac);
  h_muendx_allsel[1]->SetFillStyle(3005);
  h_muendx_allsel[1]->SetFillColor(28);
  h_muendx_allsel[1]->Scale(scale_onoffbeam);



   hs_muendx -> Add(h_muendx_sig[0]);
   hs_muendx -> Add(h_muendx_bac[0]);
   hs_muendx -> Add(h_muendx_bac[1]);
   hs_muendx -> Add(h_muendx_bac[2]);
   hs_muendx -> Add(h_muendx_bac[3]);
   hs_muendx -> Add(h_muendx_bac[4]);
   hs_muendx -> Add(h_muendx_bac[5]);
   hs_muendx -> Add(h_muendx_bac[6]);
   hs_muendx -> Add(h_muendx_bac[7]);
   hs_muendx -> Add(h_muendx_allsel[1]);
   hs_muendx -> Draw("HIST,SAME");

   h_muendx_allsel[2]->Draw("same");
   h_muendx_allsel[0]->Draw("same");
   //h_onoff_muendx->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_muendx_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_muendx_allsel.png");}
 
   //==========================================================
   TH1D *h_onoff_muendy=(TH1D*)h_muendy_allsel[1]->Clone(Form("%s_on-off", h_muendy_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_muendy_allsel[0]->SetLineColor(kBlack);
   h_muendy_allsel[0]->SetLineWidth(2);
   h_muendy_allsel[0]->SetLineStyle(1);
   h_muendy_allsel[0]->GetXaxis()->SetTitle("End Position Y of Muon Candidate[cm]");
   h_muendy_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_muendy_allsel[0]->SetMaximum(700);
   h_muendy_allsel[0]->Draw();

   //h_muendy_allsel[1]->SetLineColor(kRed);
   //h_muendy_allsel[1]->SetLineWidth(2);
   //h_muendy_allsel[1]->SetLineStyle(1);

   h_onoff_muendy->Add(h_muendy_allsel[0],h_muendy_allsel[1],1,-scale_onoffbeam);
   h_onoff_muendy->SetLineColor(kBlack);
   h_onoff_muendy->SetLineWidth(2);
   h_onoff_muendy->SetLineStyle(1);
   h_onoff_muendy->SetMaximum(300);
   h_onoff_muendy->GetXaxis()->SetTitle("End Position Y of Muon Candidate[cm]");
   h_onoff_muendy->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_muendy->DrawCopy("");

   h_muendy_allsel[2]->SetLineColor(kRed);
   h_muendy_allsel[2]->SetLineWidth(2);
   h_muendy_allsel[2]->SetLineStyle(1);
   h_muendy_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_muendy = new THStack("hs_muendy","");

   h_muendy_sig[0]-> SetFillColor(2); 
   h_muendy_sig[0]->Scale(normfac);
   h_muendy_bac[0]-> SetFillColor(46);
   h_muendy_bac[0]->Scale(normfac);
   h_muendy_bac[1]-> SetFillColor(4);
   h_muendy_bac[1]->Scale(normfac);
   h_muendy_bac[2]-> SetFillColor(8);
   h_muendy_bac[2]->Scale(normfac);
   h_muendy_bac[3] -> SetFillColor(5);
   h_muendy_bac[3]->Scale(normfac);
   h_muendy_bac[4]-> SetFillColor(3);
   h_muendy_bac[4]->Scale(normfac);
   h_muendy_bac[5]-> SetFillColor(6);
   h_muendy_bac[5]->Scale(normfac);
   h_muendy_bac[6]-> SetFillColor(7);
   h_muendy_bac[6]->Scale(normfac);
   h_muendy_bac[7]-> SetFillColor(9);
   h_muendy_bac[7]->Scale(normfac);
  h_muendy_allsel[1]->SetFillStyle(3005);
  h_muendy_allsel[1]->SetFillColor(28);
  h_muendy_allsel[1]->Scale(scale_onoffbeam);



   hs_muendy -> Add(h_muendy_sig[0]);
   hs_muendy -> Add(h_muendy_bac[0]);
   hs_muendy -> Add(h_muendy_bac[1]);
   hs_muendy -> Add(h_muendy_bac[2]);
   hs_muendy -> Add(h_muendy_bac[3]);
   hs_muendy -> Add(h_muendy_bac[4]);
   hs_muendy -> Add(h_muendy_bac[5]);
   hs_muendy -> Add(h_muendy_bac[6]);
   hs_muendy -> Add(h_muendy_bac[7]);
   hs_muendy -> Add(h_muendy_allsel[1]);
   hs_muendy -> Draw("HIST,SAME");

   h_muendy_allsel[2]->Draw("same");
   h_muendy_allsel[0]->Draw("same");
   //h_onoff_muendy->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_muendy_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_muendy_allsel.png");}
  //==================================================================== 
   TH1D *h_onoff_muendz=(TH1D*)h_muendz_allsel[1]->Clone(Form("%s_on-off", h_muendz_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_muendz_allsel[0]->SetLineColor(kBlack);
   h_muendz_allsel[0]->SetLineWidth(2);
   h_muendz_allsel[0]->SetLineStyle(1);
   h_muendz_allsel[0]->GetXaxis()->SetTitle("End Position Z of Muon Candidate[cm]");
   h_muendz_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_muendz_allsel[0]->SetMaximum(500);
   h_muendz_allsel[0]->Draw();

   //h_muendz_allsel[1]->SetLineColor(kRed);
   //h_muendz_allsel[1]->SetLineWidth(2);
   //h_muendz_allsel[1]->SetLineStyle(1);

   h_onoff_muendz->Add(h_muendz_allsel[0],h_muendz_allsel[1],1,-scale_onoffbeam);
   h_onoff_muendz->SetLineColor(kBlack);
   h_onoff_muendz->SetLineWidth(2);
   h_onoff_muendz->SetLineStyle(1);
   h_onoff_muendz->SetMaximum(300);
   h_onoff_muendz->GetXaxis()->SetTitle("End Position Z of Muon Candidate[cm]");
   h_onoff_muendz->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_muendz->DrawCopy("");

   h_muendz_allsel[2]->SetLineColor(kRed);
   h_muendz_allsel[2]->SetLineWidth(2);
   h_muendz_allsel[2]->SetLineStyle(1);
   h_muendz_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_muendz = new THStack("hs_muendz","");

   h_muendz_sig[0]-> SetFillColor(2); 
   h_muendz_sig[0]->Scale(normfac);
   h_muendz_bac[0]-> SetFillColor(46);
   h_muendz_bac[0]->Scale(normfac);
   h_muendz_bac[1]-> SetFillColor(4);
   h_muendz_bac[1]->Scale(normfac);
   h_muendz_bac[2]-> SetFillColor(8);
   h_muendz_bac[2]->Scale(normfac);
   h_muendz_bac[3] -> SetFillColor(5);
   h_muendz_bac[3]->Scale(normfac);
   h_muendz_bac[4]-> SetFillColor(3);
   h_muendz_bac[4]->Scale(normfac);
   h_muendz_bac[5]-> SetFillColor(6);
   h_muendz_bac[5]->Scale(normfac);
   h_muendz_bac[6]-> SetFillColor(7);
   h_muendz_bac[6]->Scale(normfac);
   h_muendz_bac[7]-> SetFillColor(9);
   h_muendz_bac[7]->Scale(normfac);
  h_muendz_allsel[1]->SetFillStyle(3005);
  h_muendz_allsel[1]->SetFillColor(28);
  h_muendz_allsel[1]->Scale(scale_onoffbeam);



   hs_muendz -> Add(h_muendz_sig[0]);
   hs_muendz -> Add(h_muendz_bac[0]);
   hs_muendz -> Add(h_muendz_bac[1]);
   hs_muendz -> Add(h_muendz_bac[2]);
   hs_muendz -> Add(h_muendz_bac[3]);
   hs_muendz -> Add(h_muendz_bac[4]);
   hs_muendz -> Add(h_muendz_bac[5]);
   hs_muendz -> Add(h_muendz_bac[6]);
   hs_muendz -> Add(h_muendz_bac[7]);
   hs_muendz -> Add(h_muendz_allsel[1]);
   hs_muendz -> Draw("HIST,SAME");

   h_muendz_allsel[2]->Draw("same");
   h_muendz_allsel[0]->Draw("same");
   //h_onoff_muendz->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_muendz_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_muendz_allsel.png");}
   //==================================================================   
   TH1D *h_onoff_pstartx=(TH1D*)h_pstartx_allsel[1]->Clone(Form("%s_on-off", h_pstartx_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_pstartx_allsel[0]->SetLineColor(kBlack);
   h_pstartx_allsel[0]->SetLineWidth(2);
   h_pstartx_allsel[0]->SetLineStyle(1);
   h_pstartx_allsel[0]->GetXaxis()->SetTitle("Start Position X of Proton Candidate[cm]");
   h_pstartx_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_pstartx_allsel[0]->SetMaximum(250);
   h_pstartx_allsel[0]->Draw();

   //h_pstartx_allsel[1]->SetLineColor(kRed);
   //h_pstartx_allsel[1]->SetLineWidth(2);
   //h_pstartx_allsel[1]->SetLineStyle(1);

   h_onoff_pstartx->Add(h_pstartx_allsel[0],h_pstartx_allsel[1],1,-scale_onoffbeam);
   h_onoff_pstartx->SetLineColor(kBlack);
   h_onoff_pstartx->SetLineWidth(2);
   h_onoff_pstartx->SetLineStyle(1);
   h_onoff_pstartx->SetMaximum(150);
   h_onoff_pstartx->GetXaxis()->SetTitle("Start Position X of Proton Candidate[cm]");
   h_onoff_pstartx->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_pstartx->DrawCopy("");

   h_pstartx_allsel[2]->SetLineColor(kRed);
   h_pstartx_allsel[2]->SetLineWidth(2);
   h_pstartx_allsel[2]->SetLineStyle(1);
   h_pstartx_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_pstartx = new THStack("hs_pstartx","");

   h_pstartx_sig[0]-> SetFillColor(2); 
   h_pstartx_sig[0]->Scale(normfac);
   h_pstartx_bac[0]-> SetFillColor(46);
   h_pstartx_bac[0]->Scale(normfac);
   h_pstartx_bac[1]-> SetFillColor(4);
   h_pstartx_bac[1]->Scale(normfac);
   h_pstartx_bac[2]-> SetFillColor(8);
   h_pstartx_bac[2]->Scale(normfac);
   h_pstartx_bac[3] -> SetFillColor(5);
   h_pstartx_bac[3]->Scale(normfac);
   h_pstartx_bac[4]-> SetFillColor(3);
   h_pstartx_bac[4]->Scale(normfac);
   h_pstartx_bac[5]-> SetFillColor(6);
   h_pstartx_bac[5]->Scale(normfac);
   h_pstartx_bac[6]-> SetFillColor(7);
   h_pstartx_bac[6]->Scale(normfac);
   h_pstartx_bac[7]-> SetFillColor(9);
   h_pstartx_bac[7]->Scale(normfac);
  h_pstartx_allsel[1]->SetFillStyle(3005);
  h_pstartx_allsel[1]->SetFillColor(28);
  h_pstartx_allsel[1]->Scale(scale_onoffbeam);



   hs_pstartx -> Add(h_pstartx_sig[0]);
   hs_pstartx -> Add(h_pstartx_bac[0]);
   hs_pstartx -> Add(h_pstartx_bac[1]);
   hs_pstartx -> Add(h_pstartx_bac[2]);
   hs_pstartx -> Add(h_pstartx_bac[3]);
   hs_pstartx -> Add(h_pstartx_bac[4]);
   hs_pstartx -> Add(h_pstartx_bac[5]);
   hs_pstartx -> Add(h_pstartx_bac[6]);
   hs_pstartx -> Add(h_pstartx_bac[7]);
   hs_pstartx -> Add(h_pstartx_allsel[1]);
   hs_pstartx -> Draw("HIST,SAME");

   h_pstartx_allsel[2]->Draw("same");
   h_pstartx_allsel[0]->Draw("same");
   //h_onoff_pstartx->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_pstartx_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_pstartx_allsel.png");}
 
   //==========================================================
   TH1D *h_onoff_pstarty=(TH1D*)h_pstarty_allsel[1]->Clone(Form("%s_on-off", h_pstarty_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_pstarty_allsel[0]->SetLineColor(kBlack);
   h_pstarty_allsel[0]->SetLineWidth(2);
   h_pstarty_allsel[0]->SetLineStyle(1);
   h_pstarty_allsel[0]->GetXaxis()->SetTitle("Start Position Y of Proton Candidate[cm]");
   h_pstarty_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_pstarty_allsel[0]->SetMaximum(250);
    h_pstarty_allsel[0]->Draw();
  
   //h_pstarty_allsel[1]->SetLineColor(kRed);
   //h_pstarty_allsel[1]->SetLineWidth(2);
   //h_pstarty_allsel[1]->SetLineStyle(1);

   h_onoff_pstarty->Add(h_pstarty_allsel[0],h_pstarty_allsel[1],1,-scale_onoffbeam);
   h_onoff_pstarty->SetLineColor(kBlack);
   h_onoff_pstarty->SetLineWidth(2);
   h_onoff_pstarty->SetLineStyle(1);
   h_onoff_pstarty->SetMaximum(150);
   h_onoff_pstarty->GetXaxis()->SetTitle("Start Position Y of Proton Candidate[cm]");
   h_onoff_pstarty->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_pstarty->DrawCopy("");

   h_pstarty_allsel[2]->SetLineColor(kRed);
   h_pstarty_allsel[2]->SetLineWidth(2);
   h_pstarty_allsel[2]->SetLineStyle(1);
   h_pstarty_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_pstarty = new THStack("hs_pstarty","");

   h_pstarty_sig[0]-> SetFillColor(2); 
   h_pstarty_sig[0]->Scale(normfac);
   h_pstarty_bac[0]-> SetFillColor(46);
   h_pstarty_bac[0]->Scale(normfac);
   h_pstarty_bac[1]-> SetFillColor(4);
   h_pstarty_bac[1]->Scale(normfac);
   h_pstarty_bac[2]-> SetFillColor(8);
   h_pstarty_bac[2]->Scale(normfac);
   h_pstarty_bac[3] -> SetFillColor(5);
   h_pstarty_bac[3]->Scale(normfac);
   h_pstarty_bac[4]-> SetFillColor(3);
   h_pstarty_bac[4]->Scale(normfac);
   h_pstarty_bac[5]-> SetFillColor(6);
   h_pstarty_bac[5]->Scale(normfac);
   h_pstarty_bac[6]-> SetFillColor(7);
   h_pstarty_bac[6]->Scale(normfac);
   h_pstarty_bac[7]-> SetFillColor(9);
   h_pstarty_bac[7]->Scale(normfac);
  h_pstarty_allsel[1]->SetFillStyle(3005);
  h_pstarty_allsel[1]->SetFillColor(28);
  h_pstarty_allsel[1]->Scale(scale_onoffbeam);



   hs_pstarty -> Add(h_pstarty_sig[0]);
   hs_pstarty -> Add(h_pstarty_bac[0]);
   hs_pstarty -> Add(h_pstarty_bac[1]);
   hs_pstarty -> Add(h_pstarty_bac[2]);
   hs_pstarty -> Add(h_pstarty_bac[3]);
   hs_pstarty -> Add(h_pstarty_bac[4]);
   hs_pstarty -> Add(h_pstarty_bac[5]);
   hs_pstarty -> Add(h_pstarty_bac[6]);
   hs_pstarty -> Add(h_pstarty_bac[7]);
   hs_pstarty -> Add(h_pstarty_allsel[1]);
   hs_pstarty -> Draw("HIST,SAME");

   h_pstarty_allsel[2]->Draw("same");
    h_pstarty_allsel[0]->Draw("same");
   //h_onoff_pstarty->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_pstarty_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_pstarty_allsel.png");}
  //==================================================================== 
   TH1D *h_onoff_pstartz=(TH1D*)h_pstartz_allsel[1]->Clone(Form("%s_on-off", h_pstartz_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_pstartz_allsel[0]->SetLineColor(kBlack);
   h_pstartz_allsel[0]->SetLineWidth(2);
   h_pstartz_allsel[0]->SetLineStyle(1);
   h_pstartz_allsel[0]->GetXaxis()->SetTitle("Start Position Z of Proton Candidate[cm]");
   h_pstartz_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_pstartz_allsel[0]->SetMaximum(250);
   h_pstartz_allsel[0]->Draw();

   //h_pstartz_allsel[1]->SetLineColor(kRed);
   //h_pstartz_allsel[1]->SetLineWidth(2);
   //h_pstartz_allsel[1]->SetLineStyle(1);

   h_onoff_pstartz->Add(h_pstartz_allsel[0],h_pstartz_allsel[1],1,-scale_onoffbeam);
   h_onoff_pstartz->SetLineColor(kBlack);
   h_onoff_pstartz->SetLineWidth(2);
   h_onoff_pstartz->SetLineStyle(1);
   h_onoff_pstartz->SetMaximum(150);
   h_onoff_pstartz->GetXaxis()->SetTitle("Start Position Z of Proton Candidate[cm]");
   h_onoff_pstartz->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_pstartz->DrawCopy("");

   h_pstartz_allsel[2]->SetLineColor(kRed);
   h_pstartz_allsel[2]->SetLineWidth(2);
   h_pstartz_allsel[2]->SetLineStyle(1);
   h_pstartz_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_pstartz = new THStack("hs_pstartz","");

   h_pstartz_sig[0]-> SetFillColor(2); 
   h_pstartz_sig[0]->Scale(normfac);
   h_pstartz_bac[0]-> SetFillColor(46);
   h_pstartz_bac[0]->Scale(normfac);
   h_pstartz_bac[1]-> SetFillColor(4);
   h_pstartz_bac[1]->Scale(normfac);
   h_pstartz_bac[2]-> SetFillColor(8);
   h_pstartz_bac[2]->Scale(normfac);
   h_pstartz_bac[3] -> SetFillColor(5);
   h_pstartz_bac[3]->Scale(normfac);
   h_pstartz_bac[4]-> SetFillColor(3);
   h_pstartz_bac[4]->Scale(normfac);
   h_pstartz_bac[5]-> SetFillColor(6);
   h_pstartz_bac[5]->Scale(normfac);
   h_pstartz_bac[6]-> SetFillColor(7);
   h_pstartz_bac[6]->Scale(normfac);
   h_pstartz_bac[7]-> SetFillColor(9);
   h_pstartz_bac[7]->Scale(normfac);
  h_pstartz_allsel[1]->SetFillStyle(3005);
  h_pstartz_allsel[1]->SetFillColor(28);
  h_pstartz_allsel[1]->Scale(scale_onoffbeam);



   hs_pstartz -> Add(h_pstartz_sig[0]);
   hs_pstartz -> Add(h_pstartz_bac[0]);
   hs_pstartz -> Add(h_pstartz_bac[1]);
   hs_pstartz -> Add(h_pstartz_bac[2]);
   hs_pstartz -> Add(h_pstartz_bac[3]);
   hs_pstartz -> Add(h_pstartz_bac[4]);
   hs_pstartz -> Add(h_pstartz_bac[5]);
   hs_pstartz -> Add(h_pstartz_bac[6]);
   hs_pstartz -> Add(h_pstartz_bac[7]);
   hs_pstartz -> Add(h_pstartz_allsel[1]);
   hs_pstartz -> Draw("HIST,SAME");

   h_pstartz_allsel[2]->Draw("same");
   h_pstartz_allsel[0]->Draw("same");
//h_onoff_pstartz->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_pstartz_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_pstartz_allsel.png");}
  //================================================================= 
   TH1D *h_onoff_pendx=(TH1D*)h_pendx_allsel[1]->Clone(Form("%s_on-off", h_pendx_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_pendx_allsel[0]->SetLineColor(kBlack);
   h_pendx_allsel[0]->SetLineWidth(2);
   h_pendx_allsel[0]->SetLineStyle(1);
   h_pendx_allsel[0]->GetXaxis()->SetTitle("End Position X of Proton Candidate[cm]");
   h_pendx_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_pendx_allsel[0]->SetMaximum(250);
   h_pendx_allsel[0]->Draw();

   //h_pendx_allsel[1]->SetLineColor(kRed);
   //h_pendx_allsel[1]->SetLineWidth(2);
   //h_pendx_allsel[1]->SetLineStyle(1);

   h_onoff_pendx->Add(h_pendx_allsel[0],h_pendx_allsel[1],1,-scale_onoffbeam);
   h_onoff_pendx->SetLineColor(kBlack);
   h_onoff_pendx->SetLineWidth(2);
   h_onoff_pendx->SetLineStyle(1);
   h_onoff_pendx->SetMaximum(150);
   h_onoff_pendx->GetXaxis()->SetTitle("End Position X of Proton Candidate[cm]");
   h_onoff_pendx->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_pendx->DrawCopy("");

   h_pendx_allsel[2]->SetLineColor(kRed);
   h_pendx_allsel[2]->SetLineWidth(2);
   h_pendx_allsel[2]->SetLineStyle(1);
   h_pendx_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_pendx = new THStack("hs_pendx","");

   h_pendx_sig[0]-> SetFillColor(2); 
   h_pendx_sig[0]->Scale(normfac);
   h_pendx_bac[0]-> SetFillColor(46);
   h_pendx_bac[0]->Scale(normfac);
   h_pendx_bac[1]-> SetFillColor(4);
   h_pendx_bac[1]->Scale(normfac);
   h_pendx_bac[2]-> SetFillColor(8);
   h_pendx_bac[2]->Scale(normfac);
   h_pendx_bac[3] -> SetFillColor(5);
   h_pendx_bac[3]->Scale(normfac);
   h_pendx_bac[4]-> SetFillColor(3);
   h_pendx_bac[4]->Scale(normfac);
   h_pendx_bac[5]-> SetFillColor(6);
   h_pendx_bac[5]->Scale(normfac);
   h_pendx_bac[6]-> SetFillColor(7);
   h_pendx_bac[6]->Scale(normfac);
   h_pendx_bac[7]-> SetFillColor(9);
   h_pendx_bac[7]->Scale(normfac);
  h_pendx_allsel[1]->SetFillStyle(3005);
  h_pendx_allsel[1]->SetFillColor(28);
  h_pendx_allsel[1]->Scale(scale_onoffbeam);



   hs_pendx -> Add(h_pendx_sig[0]);
   hs_pendx -> Add(h_pendx_bac[0]);
   hs_pendx -> Add(h_pendx_bac[1]);
   hs_pendx -> Add(h_pendx_bac[2]);
   hs_pendx -> Add(h_pendx_bac[3]);
   hs_pendx -> Add(h_pendx_bac[4]);
   hs_pendx -> Add(h_pendx_bac[5]);
   hs_pendx -> Add(h_pendx_bac[6]);
   hs_pendx -> Add(h_pendx_bac[7]);
   hs_pendx -> Add(h_pendx_allsel[1]);
   hs_pendx -> Draw("HIST,SAME");

   h_pendx_allsel[2]->Draw("same");
   h_pendx_allsel[0]->Draw("same");
//h_onoff_pendx->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_pendx_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_pendx_allsel.png");}
 
   //==========================================================
   TH1D *h_onoff_pendy=(TH1D*)h_pendy_allsel[1]->Clone(Form("%s_on-off", h_pendy_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_pendy_allsel[0]->SetLineColor(kBlack);
   h_pendy_allsel[0]->SetLineWidth(2);
   h_pendy_allsel[0]->SetLineStyle(1);
   h_pendy_allsel[0]->GetXaxis()->SetTitle("End Position Y of Proton Candidate[cm]");
   h_pendy_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_pendy_allsel[0]->SetMaximum(250);
    h_pendy_allsel[0]->Draw();

   //h_pendy_allsel[1]->SetLineColor(kRed);
   //h_pendy_allsel[1]->SetLineWidth(2);
   //h_pendy_allsel[1]->SetLineStyle(1);

   h_onoff_pendy->Add(h_pendy_allsel[0],h_pendy_allsel[1],1,-scale_onoffbeam);
   h_onoff_pendy->SetLineColor(kBlack);
   h_onoff_pendy->SetLineWidth(2);
   h_onoff_pendy->SetLineStyle(1);
   h_onoff_pendy->SetMaximum(150);
   h_onoff_pendy->GetXaxis()->SetTitle("End Position Y of Proton Candidate[cm]");
   h_onoff_pendy->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_pendy->DrawCopy("");

   h_pendy_allsel[2]->SetLineColor(kRed);
   h_pendy_allsel[2]->SetLineWidth(2);
   h_pendy_allsel[2]->SetLineStyle(1);
   h_pendy_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_pendy = new THStack("hs_pendy","");

   h_pendy_sig[0]-> SetFillColor(2); 
   h_pendy_sig[0]->Scale(normfac);
   h_pendy_bac[0]-> SetFillColor(46);
   h_pendy_bac[0]->Scale(normfac);
   h_pendy_bac[1]-> SetFillColor(4);
   h_pendy_bac[1]->Scale(normfac);
   h_pendy_bac[2]-> SetFillColor(8);
   h_pendy_bac[2]->Scale(normfac);
   h_pendy_bac[3] -> SetFillColor(5);
   h_pendy_bac[3]->Scale(normfac);
   h_pendy_bac[4]-> SetFillColor(3);
   h_pendy_bac[4]->Scale(normfac);
   h_pendy_bac[5]-> SetFillColor(6);
   h_pendy_bac[5]->Scale(normfac);
   h_pendy_bac[6]-> SetFillColor(7);
   h_pendy_bac[6]->Scale(normfac);
   h_pendy_bac[7]-> SetFillColor(9);
   h_pendy_bac[7]->Scale(normfac);
  h_pendy_allsel[1]->SetFillStyle(3005);
  h_pendy_allsel[1]->SetFillColor(28);
  h_pendy_allsel[1]->Scale(scale_onoffbeam);



   hs_pendy -> Add(h_pendy_sig[0]);
   hs_pendy -> Add(h_pendy_bac[0]);
   hs_pendy -> Add(h_pendy_bac[1]);
   hs_pendy -> Add(h_pendy_bac[2]);
   hs_pendy -> Add(h_pendy_bac[3]);
   hs_pendy -> Add(h_pendy_bac[4]);
   hs_pendy -> Add(h_pendy_bac[5]);
   hs_pendy -> Add(h_pendy_bac[6]);
   hs_pendy -> Add(h_pendy_bac[7]);
   hs_pendy -> Add(h_pendy_allsel[1]);
   hs_pendy -> Draw("HIST,SAME");

   h_pendy_allsel[2]->Draw("same");
  h_pendy_allsel[0]->Draw("same");
  //h_onoff_pendy->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_pendy_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_pendy_allsel.png");}
  //==================================================================== 
   TH1D *h_onoff_pendz=(TH1D*)h_pendz_allsel[1]->Clone(Form("%s_on-off", h_pendz_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_pendz_allsel[0]->SetLineColor(kBlack);
   h_pendz_allsel[0]->SetLineWidth(2);
   h_pendz_allsel[0]->SetLineStyle(1);
   h_pendz_allsel[0]->GetXaxis()->SetTitle("End Position Z of Proton Candidate[cm]");
   h_pendz_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_pendz_allsel[0]->SetMaximum(250);
   h_pendz_allsel[0]->Draw();
  
   //h_pendz_allsel[1]->SetLineColor(kRed);
   //h_pendz_allsel[1]->SetLineWidth(2);
   //h_pendz_allsel[1]->SetLineStyle(1);

   h_onoff_pendz->Add(h_pendz_allsel[0],h_pendz_allsel[1],1,-scale_onoffbeam);
   h_onoff_pendz->SetLineColor(kBlack);
   h_onoff_pendz->SetLineWidth(2);
   h_onoff_pendz->SetLineStyle(1);
   h_onoff_pendz->SetMaximum(150);
   h_onoff_pendz->GetXaxis()->SetTitle("End Position Z of Proton Candidate[cm]");
   h_onoff_pendz->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_pendz->DrawCopy("");

   h_pendz_allsel[2]->SetLineColor(kRed);
   h_pendz_allsel[2]->SetLineWidth(2);
   h_pendz_allsel[2]->SetLineStyle(1);
   h_pendz_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_pendz = new THStack("hs_pendz","");

   h_pendz_sig[0]-> SetFillColor(2); 
   h_pendz_sig[0]->Scale(normfac);
   h_pendz_bac[0]-> SetFillColor(46);
   h_pendz_bac[0]->Scale(normfac);
   h_pendz_bac[1]-> SetFillColor(4);
   h_pendz_bac[1]->Scale(normfac);
   h_pendz_bac[2]-> SetFillColor(8);
   h_pendz_bac[2]->Scale(normfac);
   h_pendz_bac[3] -> SetFillColor(5);
   h_pendz_bac[3]->Scale(normfac);
   h_pendz_bac[4]-> SetFillColor(3);
   h_pendz_bac[4]->Scale(normfac);
   h_pendz_bac[5]-> SetFillColor(6);
   h_pendz_bac[5]->Scale(normfac);
   h_pendz_bac[6]-> SetFillColor(7);
   h_pendz_bac[6]->Scale(normfac);
   h_pendz_bac[7]-> SetFillColor(9);
   h_pendz_bac[7]->Scale(normfac);
  h_pendz_allsel[1]->SetFillStyle(3005);
  h_pendz_allsel[1]->SetFillColor(28);
  h_pendz_allsel[1]->Scale(scale_onoffbeam);



   hs_pendz -> Add(h_pendz_sig[0]);
   hs_pendz -> Add(h_pendz_bac[0]);
   hs_pendz -> Add(h_pendz_bac[1]);
   hs_pendz -> Add(h_pendz_bac[2]);
   hs_pendz -> Add(h_pendz_bac[3]);
   hs_pendz -> Add(h_pendz_bac[4]);
   hs_pendz -> Add(h_pendz_bac[5]);
   hs_pendz -> Add(h_pendz_bac[6]);
   hs_pendz -> Add(h_pendz_bac[7]);
   hs_pendz -> Add(h_pendz_allsel[1]);
   hs_pendz -> Draw("HIST,SAME");

   h_pendz_allsel[2]->Draw("same");
   h_pendz_allsel[0]->Draw("same");
//h_onoff_pendz->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_pendz_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_pendz_allsel.png");}
  //==================================================
   TH1D *h_onoff_Evis=(TH1D*)h_Evis_allsel[1]->Clone(Form("%s_on-off", h_Evis_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_Evis_allsel[0]->SetLineColor(kBlack);
   h_Evis_allsel[0]->SetLineWidth(2);
   h_Evis_allsel[0]->SetLineStyle(1);
   h_Evis_allsel[0]->GetXaxis()->SetTitle("Visible Energy [GeV]");
   h_Evis_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_Evis_allsel[0]->SetMaximum(500);
   h_Evis_allsel[0]->Draw();

   //h_Evis_allsel[1]->SetLineColor(kRed);
   //h_Evis_allsel[1]->SetLineWidth(2);
   //h_Evis_allsel[1]->SetLineStyle(1);

   h_onoff_Evis->Add(h_Evis_allsel[0],h_Evis_allsel[1],1,-scale_onoffbeam);
   h_onoff_Evis->SetLineColor(kBlack);
   h_onoff_Evis->SetLineWidth(2);
   h_onoff_Evis->SetLineStyle(1);
   h_onoff_Evis->SetMaximum(180);
   h_onoff_Evis->GetXaxis()->SetTitle("Visible Energy [GeV]");
   h_onoff_Evis->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_Evis->DrawCopy("");

   h_Evis_allsel[2]->SetLineColor(kRed);
   h_Evis_allsel[2]->SetLineWidth(2);
   h_Evis_allsel[2]->SetLineStyle(1);
   h_Evis_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_Evis = new THStack("hs_Evis","");

   h_Evis_sig[0]-> SetFillColor(2); 
   h_Evis_sig[0]->Scale(normfac);
   h_Evis_bac[0]-> SetFillColor(46);
   h_Evis_bac[0]->Scale(normfac);
   h_Evis_bac[1]-> SetFillColor(4);
   h_Evis_bac[1]->Scale(normfac);
   h_Evis_bac[2]-> SetFillColor(8);
   h_Evis_bac[2]->Scale(normfac);
   h_Evis_bac[3] -> SetFillColor(5);
   h_Evis_bac[3]->Scale(normfac);
   h_Evis_bac[4]-> SetFillColor(3);
   h_Evis_bac[4]->Scale(normfac);
   h_Evis_bac[5]-> SetFillColor(6);
   h_Evis_bac[5]->Scale(normfac);
   h_Evis_bac[6]-> SetFillColor(7);
   h_Evis_bac[6]->Scale(normfac);
   h_Evis_bac[7]-> SetFillColor(9);
   h_Evis_bac[7]->Scale(normfac);
  h_Evis_allsel[1]->SetFillStyle(3005);
  h_Evis_allsel[1]->SetFillColor(28);
  h_Evis_allsel[1]->Scale(scale_onoffbeam);



   hs_Evis -> Add(h_Evis_sig[0]);
   hs_Evis -> Add(h_Evis_bac[0]);
   hs_Evis -> Add(h_Evis_bac[1]);
   hs_Evis -> Add(h_Evis_bac[2]);
   hs_Evis -> Add(h_Evis_bac[3]);
   hs_Evis -> Add(h_Evis_bac[4]);
   hs_Evis -> Add(h_Evis_bac[5]);
   hs_Evis -> Add(h_Evis_bac[6]);
   hs_Evis -> Add(h_Evis_bac[7]);
   hs_Evis -> Add(h_Evis_allsel[1]);
   hs_Evis -> Draw("HIST,SAME");

   h_Evis_allsel[2]->Draw("same");
    h_Evis_allsel[0]->Draw("same");
//h_onoff_Evis->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_Evis_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_Evis_allsel.png");}
  //============================================================
   TH1D *h_onoff_Q2cal=(TH1D*)h_Q2cal_allsel[1]->Clone(Form("%s_on-off", h_Q2cal_allsel[1]->GetName()));
   cout<<"get all the histograms!"<<endl;


   h_Q2cal_allsel[0]->SetLineColor(kBlack);
   h_Q2cal_allsel[0]->SetLineWidth(2);
   h_Q2cal_allsel[0]->SetLineStyle(1);
   h_Q2cal_allsel[0]->GetXaxis()->SetTitle("Momentum Transfer Q^{2} [GeV^{2}]");
   h_Q2cal_allsel[0]->GetYaxis()->SetTitle("No. of Events");
   h_Q2cal_allsel[0]->SetMaximum(600);
   h_Q2cal_allsel[0]->Draw();

   //h_Q2cal_allsel[1]->SetLineColor(kRed);
   //h_Q2cal_allsel[1]->SetLineWidth(2);
   //h_Q2cal_allsel[1]->SetLineStyle(1);

   h_onoff_Q2cal->Add(h_Q2cal_allsel[0],h_Q2cal_allsel[1],1,-scale_onoffbeam);
   h_onoff_Q2cal->SetLineColor(kBlack);
   h_onoff_Q2cal->SetLineWidth(2);
   h_onoff_Q2cal->SetLineStyle(1);
   h_onoff_Q2cal->SetMaximum(240);
   h_onoff_Q2cal->GetXaxis()->SetTitle("Momentum Transfer Q^{2} [GeV^{2}]");
   h_onoff_Q2cal->GetYaxis()->SetTitle("No. of Events");
   //h_onoff_Q2cal->DrawCopy("");

   h_Q2cal_allsel[2]->SetLineColor(kRed);
   h_Q2cal_allsel[2]->SetLineWidth(2);
   h_Q2cal_allsel[2]->SetLineStyle(1);
   h_Q2cal_allsel[2]->Scale(normfac);
   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   THStack *hs_Q2cal = new THStack("hs_Q2cal","");

   h_Q2cal_sig[0]-> SetFillColor(2); 
   h_Q2cal_sig[0]->Scale(normfac);
   h_Q2cal_bac[0]-> SetFillColor(46);
   h_Q2cal_bac[0]->Scale(normfac);
   h_Q2cal_bac[1]-> SetFillColor(4);
   h_Q2cal_bac[1]->Scale(normfac);
   h_Q2cal_bac[2]-> SetFillColor(8);
   h_Q2cal_bac[2]->Scale(normfac);
   h_Q2cal_bac[3] -> SetFillColor(5);
   h_Q2cal_bac[3]->Scale(normfac);
   h_Q2cal_bac[4]-> SetFillColor(3);
   h_Q2cal_bac[4]->Scale(normfac);
   h_Q2cal_bac[5]-> SetFillColor(6);
   h_Q2cal_bac[5]->Scale(normfac);
   h_Q2cal_bac[6]-> SetFillColor(7);
   h_Q2cal_bac[6]->Scale(normfac);
   h_Q2cal_bac[7]-> SetFillColor(9);
   h_Q2cal_bac[7]->Scale(normfac);
   h_Q2cal_allsel[1]->SetFillStyle(3005);
   h_Q2cal_allsel[1]->SetFillColor(28);
   h_Q2cal_allsel[1]->Scale(scale_onoffbeam);



   hs_Q2cal -> Add(h_Q2cal_sig[0]);
   hs_Q2cal -> Add(h_Q2cal_bac[0]);
   hs_Q2cal -> Add(h_Q2cal_bac[1]);
   hs_Q2cal -> Add(h_Q2cal_bac[2]);
   hs_Q2cal -> Add(h_Q2cal_bac[3]);
   hs_Q2cal -> Add(h_Q2cal_bac[4]);
   hs_Q2cal -> Add(h_Q2cal_bac[5]);
   hs_Q2cal -> Add(h_Q2cal_bac[6]);
   hs_Q2cal -> Add(h_Q2cal_bac[7]);
   hs_Q2cal -> Add(h_Q2cal_allsel[1]);
   hs_Q2cal -> Draw("HIST,SAME");

   h_Q2cal_allsel[2]->Draw("same");
    h_Q2cal_allsel[0]->Draw("same");
   //h_onoff_Q2cal->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
   legend->Draw("same");
   if(tune==3){c1->Print("figures/Tune3/BackSep/h_Q2cal_allsel.png");}
   else{c1->Print("figures/Tune1/BackSep/h_Q2cal_allsel.png");}




}
