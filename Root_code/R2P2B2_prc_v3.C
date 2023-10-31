/**
 * @Author: Prof. B. K. Nandi, Prof. P. R. Pujahari, Prof. Claude A Pruneau & Baidyanath Sahoo
 * @Date:   2017-06-18T23:15:25+05:30
 * @Filename: R2P2B2.C
 * @Extracted From: CrossSecScale.C
 */
#include <iostream>
#include <fstream>
#include <TH2D.h>
#include <TH1.h>
#include <string>
#include <TMath.h>
#include <vector>

#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TBuffer.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
using namespace std;

void Normalize2D(TH2D *hh, Double_t scale, Double_t dphi, Double_t deta){
	double norm = scale / (deta * dphi);
	hh->Scale(norm);
}



//----------------------R2----------------------------------
Double_t Rho1Rho1EtaPhi(const TH2D *h_1, const TH2D * h_2, TH2D * h_12);

void CalculateR2(const TH2D * n2_12, const TH2D * n1n1_12, TH2D * r2_12);
void EtaPhiTodEtadPhi(const TH2D * source, TH2D * target, Int_t nEtaBins, Int_t nPhiBins);
void SymmetrizeDEtaDPhi(TH2D *hin);
void shiftY(const TH2D & source, TH2D & target, Int_t nbins);

//--------------DpTDpT-------------------------------

bool sameDimensions2D(const TH2D* h1, const TH2D* h2);
void calculateDptDpt(const TH2D * spp, const TH2D * spn, const TH2D * snp, const TH2D * snn, const TH1D * avgpt1, const TH1D * avgpt2, const Float_t kMinPt, const Float_t kMaxPt, TH2D * p2DptDpt,  TH2D * dptdpt, Int_t nEta, Int_t nPhi);
//-------------calculate CI & CD -------------------
void calculate_CI_CD(const TH2D *hpm, const TH2D *hpp, const TH2D *hmm,  TH2D *hci,  TH2D *hcd);
//void smoothFunXRebinY(const TH2D *source, const TH2D *target);
void smoothFunXRebinY(TH2D *source, TH2D *target);

//-------------calculate B2 -------------------
void calculate_B2( const TH2D *hp, const TH2D *hm, const TH2D *hin,  TH2D *hout);

//typedef enum {kLowPt, kMidPt, kHighPt} switchPtCut_t;

//void R2P2B2_prc_v3(Int_t switchPtCut = kLowPt)
void R2P2B2_prc_v3(/*const TString ptCase, const TString strMult, const TString strPt*/)
{
  cout <<"\n ========== R2P2B2_prc_v3 starts ======= \n" <<endl;
//  cout <<"R2P2B2_prc_v3: "<<ptCase;<<"\t"<<strMult<<"\t"<< strPt<<endl;
  
  //--------------Initialization Of Varibles-----------------

  //Double_t kminpt = 0.0, kmaxpt = 0.0;

  Float_t kminpt  = 0, kmaxpt = 0; 
  int ptCut = 3;//atoi(ptCase); //convert string to int
  switch (ptCut)
    {
    case 0:
      kminpt  = 0.2, kmaxpt = 2.0;
      break;
    case 1:
      kminpt  = 2.0, kmaxpt = 5.0;
      break;
    case 2:
      kminpt = 5.0, kmaxpt = 30.0;
      break;
    case 3:
      kminpt = 1.0, kmaxpt = 3.0;
      break;
    }
    kminpt=0.2;
    kmaxpt=2.0;
  
  cout <<"\nPt Range: "<<kminpt<<"\t"<< kmaxpt<<endl;


  
  //----------Getter-----------------------------------------
  //TFile *ftpcorr      = TFile::Open("TwoParticleCorrScaled.root");
  TFile *ftpcorr      = TFile::Open(/*"TwoParticleCorrScaled_"+strMult+"_"+strPt+"*/"AnalysisResults.root", "read");
  //TFile *ftpcorr      = 
  //TFile::Open("ScaledPythia_qq2qq_HighPt.root");
  if (!ftpcorr) return;

  /*TH1D *h1i_n1_multPM = (TH1D*)ftpcorr->Get("h1i_n1_multPM");
  TH1D *h1d_n1_etaPM  = (TH1D*)ftpcorr->Get("h1d_n1_etaPM");
  TH1D *h1d_n1_phiPM  = (TH1D*)ftpcorr->Get("h1d_n1_phiPM");
  TH1D *h1d_n1_etaP   = (TH1D*)ftpcorr->Get("h1d_n1_etaP");
  TH1D *h1d_n1_phiP   = (TH1D*)ftpcorr->Get("h1d_n1_phiP");
  TH1D *h1d_n1_ptP    = (TH1D*)ftpcorr->Get("h1d_n1_ptP");
  TH1D *h1d_n1_etaM   = (TH1D*)ftpcorr->Get("h1d_n1_etaM");
  TH1D *h1d_n1_phiM   = (TH1D*)ftpcorr->Get("h1d_n1_phiM");
  TH1D *h1d_n1_ptM    = (TH1D*)ftpcorr->Get("h1d_n1_ptM");
  */
  TH1D *h1d_n1_ptP    = (TH1D*)ftpcorr->Get("h1d_n1_ptP");
  TH1D *h1d_n1_ptM    = (TH1D*)ftpcorr->Get("h1d_n1_ptM");

  //----------R2----------------------------
  TH2D *h2d_n1_etaPhiP             = (TH2D*)ftpcorr->Get("h2d_n1_etaPhiP");
  TH2D *h2d_n1_etaPhiM             = (TH2D*)ftpcorr->Get("h2d_n1_etaPhiM");
  TH2D *h2d_n2_eta1Phi1Eta2Phi2PM  = (TH2D*)ftpcorr->Get("h2d_n2_eta1Phi1Eta2Phi2PM");
  TH2D *h2d_n2_eta1Phi1Eta2Phi2PP  = (TH2D*)ftpcorr->Get("h2d_n2_eta1Phi1Eta2Phi2PP");
  TH2D *h2d_n2_eta1Phi1Eta2Phi2MM  = (TH2D*)ftpcorr->Get("h2d_n2_eta1Phi1Eta2Phi2MM");

  //--------------------DpTDpT---------------------
  TH2D *h2d_pt_etaPhiP               = (TH2D*)ftpcorr->Get("h2d_pt_etaPhiP");
  TH2D *h2d_pt_etaPhiM               = (TH2D*)ftpcorr->Get("h2d_pt_etaPhiM");
  TH2D *h2d_ptpt_eta1Phi1Eta2Phi2PM  = (TH2D*)ftpcorr->Get("h2d_ptpt_eta1Phi1Eta2Phi2PM");
  TH2D *h2d_ptpt_eta1Phi1Eta2Phi2PP  = (TH2D*)ftpcorr->Get("h2d_ptpt_eta1Phi1Eta2Phi2PP");
  TH2D *h2d_ptpt_eta1Phi1Eta2Phi2MM  = (TH2D*)ftpcorr->Get("h2d_ptpt_eta1Phi1Eta2Phi2MM");
  TH2D *h2d_ptn_eta1Phi1Eta2Phi2PM   = (TH2D*)ftpcorr->Get("h2d_ptn_eta1Phi1Eta2Phi2PM");
  TH2D *h2d_ptn_eta1Phi1Eta2Phi2PP   = (TH2D*)ftpcorr->Get("h2d_ptn_eta1Phi1Eta2Phi2PP");
  TH2D *h2d_ptn_eta1Phi1Eta2Phi2MM   = (TH2D*)ftpcorr->Get("h2d_ptn_eta1Phi1Eta2Phi2MM");
  TH2D *h2d_npt_eta1Phi1Eta2Phi2PM   = (TH2D*)ftpcorr->Get("h2d_npt_eta1Phi1Eta2Phi2PM");
  TH2D *h2d_npt_eta1Phi1Eta2Phi2PP   = (TH2D*)ftpcorr->Get("h2d_npt_eta1Phi1Eta2Phi2PP");
  TH2D *h2d_npt_eta1Phi1Eta2Phi2MM   = (TH2D*)ftpcorr->Get("h2d_npt_eta1Phi1Eta2Phi2MM");




  const Double_t kmineta    = -0.8;
  const Double_t kmaxeta    =  0.8;
  const Double_t kmindeta   = 2.0 * kmineta;
  const Double_t kmaxdeta   = 2.0 * kmaxeta;

  const Double_t kminphi    = 0.0;
  const Double_t kmaxphi    = 2.*TMath::Pi();




  const Int_t knetabins     = h2d_n1_etaPhiP->GetNbinsX() ;
  const Int_t knphibins     = h2d_n1_etaPhiP->GetNbinsY() ;

  const Int_t nBins_ShiftY  = knphibins / 4.0;
  //cout << "\v \v nBins_Dphi_shft \v \v" << nBins_ShiftY << endl;


  Double_t deta             = h2d_n1_etaPhiP->GetXaxis()->GetBinWidth(0);
  Double_t deta2            = TMath::Power(deta, 2);
  Double_t dphi             = h2d_n1_etaPhiP->GetYaxis()->GetBinWidth(0);
  Double_t dphi2            = TMath::Power(dphi, 2);

  const Int_t knetaphibinsx = knphibins * knetabins;
  const Int_t knetaphibinsy = knetaphibinsx;

  const Double_t kminetaphi = 0.;
  const Double_t kmaxetaphi = (Double_t) knetaphibinsx;
  const Int_t nBins_Deta    = 2.*knetabins - 1;

  //-------------------Histogram Booking----------------------
  //-----------------R2-----------------------------------------
  TH2D *h2d_n1n1_eta1Phi1Eta2Phi2PM
    = new TH2D("h2d_n1n1_eta1Phi1Eta2Phi2PM",
               "h2d_n1n1_eta1Phi1Eta2Phi2PM",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_n1n1_eta1Phi1Eta2Phi2PP
    = new TH2D("h2d_n1n1_eta1Phi1Eta2Phi2PP",
               "h2d_n1n1_eta1Phi1Eta2Phi2PP",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_n1n1_eta1Phi1Eta2Phi2MM
    = new TH2D("h2d_n1n1_eta1Phi1Eta2Phi2MM",
               "h2d_n1n1_eta1Phi1Eta2Phi2MM",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);

  TH2D *h2d_r2_eta1Phi1Eta2Phi2PM
    = new TH2D("h2d_r2_eta1Phi1Eta2Phi2PM",
               "h2d_r2_eta1Phi1Eta2Phi2PM",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_r2_eta1Phi1Eta2Phi2PP
    = new TH2D("h2d_r2_eta1Phi1Eta2Phi2PP",
               "h2d_r2_eta1Phi1Eta2Phi2PP",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_r2_eta1Phi1Eta2Phi2MM
    = new TH2D("h2d_r2_eta1Phi1Eta2Phi2MM",
               "h2d_r2_eta1Phi1Eta2Phi2MM",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);


  TH2D *h2d_n1n1_DetaDphiPM
    = new TH2D("h2d_n1n1_DetaDphiPM",
               "h2d_n1n1_DetaDphiPM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_n1n1_DetaDphiPP
    = new TH2D("h2d_n1n1_DetaDphiPP",
               "h2d_n1n1_DetaDphiPP",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_n1n1_DetaDphiMM
    = new TH2D("h2d_n1n1_DetaDphiMM",
               "h2d_n1n1_DetaDphiMM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_n2_DetaDphiPM
    = new TH2D("h2d_n2_DetaDphiPM",
               "h2d_n2_DetaDphiPM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_n2_DetaDphiPP
    = new TH2D("h2d_n2_DetaDphiPP",
               "h2d_n2_DetaDphiPP",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_n2_DetaDphiMM
    = new TH2D("h2d_n2_DetaDphiMM",
               "h2d_n2_DetaDphiMM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_r2_DetaDphiPM
    = new TH2D("h2d_r2_DetaDphiPM",
               "h2d_r2_DetaDphiPM;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_r2_DetaDphiPP
    = new TH2D("h2d_r2_DetaDphiPP",
               "h2d_r2_DetaDphiPP;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_r2_DetaDphiMM
    = new TH2D("h2d_r2_DetaDphiMM",
               "h2d_r2_DetaDphiMM;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_n1n1_DetaDphiPM_ShiftY
    = new TH2D("h2d_n1n1_DetaDphiPM_ShiftY",
               "h2d_n1n1_DetaDphiPM_ShiftY",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_n1n1_DetaDphiPP_ShiftY
    = new TH2D("h2d_n1n1_DetaDphiPP_ShiftY",
               "h2d_n1n1_DetaDphiPP_ShiftY",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_n1n1_DetaDphiMM_ShiftY
    = new TH2D("h2d_n1n1_DetaDphiMM_ShiftY",
               "h2d_n1n1_DetaDphiMM_ShiftY",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_n2_DetaDphiPM_ShiftY
    = new TH2D("h2d_n2_DetaDphiPM_ShiftY",
               "h2d_n2_DetaDphiPM_ShiftY",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_n2_DetaDphiPP_ShiftY =
    new TH2D("h2d_n2_DetaDphiPP_ShiftY",
             "h2d_n2_DetaDphiPP_ShiftY",
             nBins_Deta, kmindeta, kmaxdeta,
             knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_n2_DetaDphiMM_ShiftY
    = new TH2D("h2d_n2_DetaDphiMM_ShiftY",
               "h2d_n2_DetaDphiMM_ShiftY",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_r2_DetaDphiPM_ShiftY
    = new TH2D("h2d_r2_DetaDphiPM_ShiftY",
               "h2d_r2_DetaDphiPM_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_r2_DetaDphiPP_ShiftY
    = new TH2D("h2d_r2_DetaDphiPP_ShiftY",
               "h2d_r2_DetaDphiPP_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_r2_DetaDphiMM_ShiftY
    = new TH2D("h2d_r2_DetaDphiMM_ShiftY",
               "h2d_r2_DetaDphiMM_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_r2_DetaDphiLS
    = new TH2D("h2d_r2_DetaDphiLS",
               "h2d_r2_DetaDphiLS;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());


  TH2D *h2d_r2CI_DetaDphi
    = new TH2D("h2d_r2CI_DetaDphi",
               "h2d_r2CI_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_r2CD_DetaDphi
    = new TH2D("h2d_r2CD_DetaDphi",
               "h2d_r2CD_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_B2_DetaDphi
    = new TH2D("h2d_B2_DetaDphi",
               "h2d_B2_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());


  //-----------------DpTDpT----------------------------

  TH2D *h2d_pt1pt1_eta1Phi1Eta2Phi2PM
    = new TH2D("h2d_pt1pt1_eta1Phi1Eta2Phi2PM",
               "h2d_pt1pt1_eta1Phi1Eta2Phi2PM",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_pt1pt1_eta1Phi1Eta2Phi2PP
    = new TH2D("h2d_pt1pt1_eta1Phi1Eta2Phi2PP",
               "h2d_pt1pt1_eta1Phi1Eta2Phi2PP",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_pt1pt1_eta1Phi1Eta2Phi2MM
    = new TH2D("h2d_pt1pt1_eta1Phi1Eta2Phi2MM",
               "h2d_pt1pt1_eta1Phi1Eta2Phi2MM",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);

  TH2D *h2d_DptDpt_eta1Phi1Eta2Phi2PM
    = new TH2D("h2d_DptDpt_eta1Phi1Eta2Phi2PM",
               "h2d_DptDpt_eta1Phi1Eta2Phi2PM",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_DptDpt_eta1Phi1Eta2Phi2PP
    = new TH2D("h2d_DptDpt_eta1Phi1Eta2Phi2PP",
               "h2d_DptDpt_eta1Phi1Eta2Phi2PP",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_DptDpt_eta1Phi1Eta2Phi2MM
    = new TH2D("h2d_DptDpt_eta1Phi1Eta2Phi2MM",
               "h2d_DptDpt_eta1Phi1Eta2Phi2MM",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);

  TH2D *h2d_p2DptDpt_eta1Phi1Eta2Phi2PM
    = new TH2D("h2d_p2DptDpt_eta1Phi1Eta2Phi2PM",
               "h2d_p2DptDpt_eta1Phi1Eta2Phi2PM;deta;dphi",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_p2DptDpt_eta1Phi1Eta2Phi2PP
    = new TH2D("h2d_p2DptDpt_eta1Phi1Eta2Phi2PP",
               "h2d_p2DptDpt_eta1Phi1Eta2Phi2PP;deta;dphi",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);
  TH2D *h2d_p2DptDpt_eta1Phi1Eta2Phi2MM
    = new TH2D("h2d_p2DptDpt_eta1Phi1Eta2Phi2MM",
               "h2d_p2DptDpt_eta1Phi1Eta2Phi2MM;deta;dphi",
               knetaphibinsx, kminetaphi, kmaxetaphi,
               knetaphibinsy, kminetaphi, kmaxetaphi);

  TH2D *h2d_pt1pt1_DetaDphiPM
    = new TH2D("h2d_pt1pt1_DetaDphiPM",
               "h2d_pt1pt1_DetaDphiPM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_pt1pt1_DetaDphiPP
    = new TH2D("h2d_pt1pt1_DetaDphiPP",
               "h2d_pt1pt1_DetaDphiPP",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_pt1pt1_DetaDphiMM
    = new TH2D("h2d_pt1pt1_DetaDphiMM",
               "h2d_pt1pt1_DetaDphiMM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_ptpt_DetaDphiPM
    = new TH2D("h2d_ptpt_DetaDphiPM",
               "h2d_ptpt_DetaDphiPM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_ptpt_DetaDphiPP
    = new TH2D("h2d_ptpt_DetaDphiPP",
               "h2d_ptpt_DetaDphiPP",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_ptpt_DetaDphiMM
    = new TH2D("h2d_ptpt_DetaDphiMM",
               "h2d_ptpt_DetaDphiMM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_ptn_DetaDphiPM
    = new TH2D("h2d_ptn_DetaDphiPM",
               "h2d_ptn_DetaDphiPM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_ptn_DetaDphiPP
    = new TH2D("h2d_ptn_DetaDphiPP",
               "h2d_ptn_DetaDphiPP",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_ptn_DetaDphiMM
    = new TH2D("h2d_ptn_DetaDphiMM",
               "h2d_ptn_DetaDphiMM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_npt_DetaDphiPM
    = new TH2D("h2d_npt_DetaDphiPM",
               "h2d_npt_DetaDphiPM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_npt_DetaDphiPP
    = new TH2D("h2d_npt_DetaDphiPP",
               "h2d_npt_DetaDphiPP",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_npt_DetaDphiMM
    = new TH2D("h2d_npt_DetaDphiMM",
               "h2d_npt_DetaDphiMM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_DptDpt_DetaDphiPM
    = new TH2D("h2d_DptDpt_DetaDphiPM",
               "h2d_DptDpt_DetaDphiPM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_DptDpt_DetaDphiPP
    = new TH2D("h2d_DptDpt_DetaDphiPP",
               "h2d_DptDpt_DetaDphiPP",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());
  TH2D *h2d_DptDpt_DetaDphiMM
    = new TH2D("h2d_DptDpt_DetaDphiMM",
               "h2d_DptDpt_DetaDphiMM",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_DptDpt_DetaDphiPM_ShiftY
    = new TH2D("h2d_DptDpt_DetaDphiPM_ShiftY",
               "h2d_DptDpt_DetaDphiPM_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_DptDpt_DetaDphiPP_ShiftY
    = new TH2D("h2d_DptDpt_DetaDphiPP_ShiftY",
               "h2d_DptDpt_DetaDphiPP_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_DptDpt_DetaDphiMM_ShiftY
    = new TH2D("h2d_DptDpt_DetaDphiMM_ShiftY",
               "h2d_DptDpt_DetaDphiMM_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_DptDptCI_DetaDphi
    = new TH2D("h2d_DptDptCI_DetaDphi",
               "h2d_DptDptCI_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_DptDptCD_DetaDphi
    = new TH2D("h2d_DptDptCD_DetaDphi",
               "h2d_DptDptCD_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_B2DptDpt_DetaDphi
    = new TH2D("h2d_B2DptDpt_DetaDphi",
               "h2d_B2DptDpt_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_p2DptDpt_DetaDphiPM
    = new TH2D("h2d_p2DptDpt_DetaDphiPM",
               "h2d_p2DptDpt_DetaDphiPM;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_p2DptDpt_DetaDphiPP
    = new TH2D("h2d_p2DptDpt_DetaDphiPP",
               "h2d_p2DptDpt_DetaDphiPP;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_p2DptDpt_DetaDphiMM
    = new TH2D("h2d_p2DptDpt_DetaDphiMM",
               "h2d_p2DptDpt_DetaDphiMM;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, 0., 2.*TMath::Pi());

  TH2D *h2d_p2DptDpt_DetaDphiPM_ShiftY
    = new TH2D("h2d_p2DptDpt_DetaDphiPM_ShiftY",
               "h2d_p2DptDpt_DetaDphiPM_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_p2DptDpt_DetaDphiPP_ShiftY
    = new TH2D("h2d_p2DptDpt_DetaDphiPP_ShiftY",
               "h2d_p2DptDpt_DetaDphiPP_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_p2DptDpt_DetaDphiMM_ShiftY
    = new TH2D("h2d_p2DptDpt_DetaDphiMM_ShiftY",
               "h2d_p2DptDpt_DetaDphiMM_ShiftY;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_p2DptDpt_DetaDphiLS
    = new TH2D("h2d_p2DptDpt_DetaDphiLS",
               "h2d_p2DptDpt_DetaDphiLS;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_p2DptDptCI_DetaDphi
    = new TH2D("h2d_p2DptDptCI_DetaDphi",
               "h2d_p2DptDptCI_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D *h2d_p2DptDptCD_DetaDphi
    = new TH2D("h2d_p2DptDptCD_DetaDphi",
               "h2d_p2DptDptCD_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
TH2D *h2d_p2DptDptCD_DetaDphiRebined
    = new TH2D("h2d_p2DptDptCD_DetaDphiRebined",
               "h2d_p2DptDptCD_DetaDphiRebined;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  TH2D *h2d_B2p2DptDpt_DetaDphi
    = new TH2D("h2d_B2p2DptDpt_DetaDphi",
               "h2d_B2p2DptDpt_DetaDphi;deta;dphi",
               nBins_Deta, kmindeta, kmaxdeta,
               knphibins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
//-------------------------------------------------------------------------------------------------------
TH1I *h1i_n1_multPM = (TH1I*)ftpcorr->Get("h1i_n1_multPM");
	double sigmaFrac = 1.0;
	double nevt = h1i_n1_multPM->GetEntries();
	double scale = sigmaFrac / nevt;
	Normalize2D(h2d_n1_etaPhiP, scale, dphi, deta);
	Normalize2D(h2d_n1_etaPhiM, scale, dphi, deta);
	Normalize2D(h2d_n2_eta1Phi1Eta2Phi2PM, scale, dphi2, deta2);
	Normalize2D(h2d_n2_eta1Phi1Eta2Phi2PP, scale, dphi2, deta2);
	Normalize2D(h2d_n2_eta1Phi1Eta2Phi2MM, scale, dphi2, deta2);
	Normalize2D(h2d_pt_etaPhiP, scale, dphi, deta);
	Normalize2D(h2d_pt_etaPhiM, scale, dphi, deta);
	Normalize2D(h2d_npt_eta1Phi1Eta2Phi2PM, scale, dphi2, deta2);
	Normalize2D(h2d_npt_eta1Phi1Eta2Phi2PP, scale, dphi2, deta2);
	Normalize2D(h2d_npt_eta1Phi1Eta2Phi2MM, scale, dphi2, deta2);
	Normalize2D(h2d_ptpt_eta1Phi1Eta2Phi2PM, scale, dphi2, deta2);
	Normalize2D(h2d_ptpt_eta1Phi1Eta2Phi2PP, scale, dphi2, deta2);
	Normalize2D(h2d_ptpt_eta1Phi1Eta2Phi2MM, scale, dphi2, deta2);
	Normalize2D(h2d_ptn_eta1Phi1Eta2Phi2PM, scale, dphi2, deta2);
	Normalize2D(h2d_ptn_eta1Phi1Eta2Phi2PP, scale, dphi2, deta2);
	Normalize2D(h2d_ptn_eta1Phi1Eta2Phi2MM, scale, dphi2, deta2);

//-----------------_______________________________--------------------------------------___________________--------------------*/

  //--------R2---------------------------------

  h2d_n1n1_eta1Phi1Eta2Phi2PM->Sumw2();
  h2d_n1n1_eta1Phi1Eta2Phi2PP->Sumw2();
  h2d_n1n1_eta1Phi1Eta2Phi2MM->Sumw2();
  h2d_n1n1_DetaDphiPM->Sumw2();
  h2d_n1n1_DetaDphiPP->Sumw2();
  h2d_n1n1_DetaDphiMM->Sumw2();
  h2d_n2_DetaDphiPM->Sumw2();
  h2d_n2_DetaDphiPP->Sumw2();
  h2d_n2_DetaDphiMM->Sumw2();
  h2d_r2_eta1Phi1Eta2Phi2PM->Sumw2();
  h2d_r2_eta1Phi1Eta2Phi2PP->Sumw2();
  h2d_r2_eta1Phi1Eta2Phi2MM->Sumw2();
  //h2d_n1n1_DetaDphiPM_ShiftY->Sumw2();
  //h2d_n1n1_DetaDphiPP_ShiftY->Sumw2();
  //h2d_n1n1_DetaDphiMM_ShiftY->Sumw2();
  //h2d_n2_DetaDphiPM_ShiftY->Sumw2();
  //h2d_n2_DetaDphiPP_ShiftY->Sumw2();
  //h2d_n2_DetaDphiMM_ShiftY->Sumw2();
  //h2d_r2_DetaDphiPM_ShiftY->Sumw2();
  //h2d_r2_DetaDphiPP_ShiftY->Sumw2();
  //h2d_r2_DetaDphiMM_ShiftY->Sumw2();
  h2d_r2_DetaDphiLS->Sumw2();
  h2d_r2CI_DetaDphi->Sumw2();
  h2d_r2CD_DetaDphi->Sumw2();
  //h2d_B2_DetaDphi->Sumw2();

  //-----------------DpTDpT----------------------
  h2d_pt1pt1_eta1Phi1Eta2Phi2PM->Sumw2();
  h2d_pt1pt1_eta1Phi1Eta2Phi2PP->Sumw2();
  h2d_pt1pt1_eta1Phi1Eta2Phi2MM->Sumw2();
  h2d_DptDpt_eta1Phi1Eta2Phi2PM->Sumw2();
  h2d_DptDpt_eta1Phi1Eta2Phi2PP->Sumw2();
  h2d_DptDpt_eta1Phi1Eta2Phi2MM->Sumw2();
  h2d_p2DptDpt_eta1Phi1Eta2Phi2PM->Sumw2();
  h2d_p2DptDpt_eta1Phi1Eta2Phi2PP->Sumw2();
  h2d_p2DptDpt_eta1Phi1Eta2Phi2MM->Sumw2();
  h2d_pt1pt1_DetaDphiPM->Sumw2();
  h2d_pt1pt1_DetaDphiPP->Sumw2();
  h2d_pt1pt1_DetaDphiMM->Sumw2();
  h2d_ptn_DetaDphiPM->Sumw2();
  h2d_ptn_DetaDphiPP->Sumw2();
  h2d_ptn_DetaDphiMM->Sumw2();
  h2d_npt_DetaDphiPM->Sumw2();
  h2d_npt_DetaDphiPP->Sumw2();
  h2d_npt_DetaDphiMM->Sumw2();
  h2d_DptDpt_DetaDphiPM->Sumw2();
  h2d_DptDpt_DetaDphiPP->Sumw2();
  h2d_DptDpt_DetaDphiMM->Sumw2();
  h2d_DptDpt_DetaDphiPM_ShiftY->Sumw2();
  h2d_DptDpt_DetaDphiPP_ShiftY->Sumw2();
  h2d_DptDpt_DetaDphiMM_ShiftY->Sumw2();
  h2d_p2DptDpt_DetaDphiPM->Sumw2();
  h2d_p2DptDpt_DetaDphiPP->Sumw2();
  h2d_p2DptDpt_DetaDphiMM->Sumw2();
  //h2d_p2DptDpt_DetaDphiPM_ShiftY->Sumw2();
  //h2d_p2DptDpt_DetaDphiPP_ShiftY->Sumw2();
  //h2d_p2DptDpt_DetaDphiMM_ShiftY->Sumw2();
  h2d_p2DptDpt_DetaDphiLS->Sumw2();
  h2d_DptDptCI_DetaDphi->Sumw2();
  h2d_DptDptCD_DetaDphi->Sumw2();
  //h2d_B2DptDpt_DetaDphi->Sumw2();
  h2d_p2DptDptCI_DetaDphi->Sumw2();
  h2d_p2DptDptCD_DetaDphi->Sumw2();
  h2d_p2DptDptCD_DetaDphiRebined->Sumw2();
  //h2d_B2p2DptDpt_DetaDphi->Sumw2();
 //product of single particle density
  //---------------R2---------------------------------
  Rho1Rho1EtaPhi(h2d_n1_etaPhiP, h2d_n1_etaPhiM, h2d_n1n1_eta1Phi1Eta2Phi2PM);
  Rho1Rho1EtaPhi(h2d_n1_etaPhiP, h2d_n1_etaPhiP, h2d_n1n1_eta1Phi1Eta2Phi2PP);
  Rho1Rho1EtaPhi(h2d_n1_etaPhiM, h2d_n1_etaPhiM, h2d_n1n1_eta1Phi1Eta2Phi2MM);

  //-------------------DpTDpT---------------------------------
  Rho1Rho1EtaPhi(h2d_pt_etaPhiP, h2d_pt_etaPhiM, h2d_pt1pt1_eta1Phi1Eta2Phi2PM);
  Rho1Rho1EtaPhi(h2d_pt_etaPhiP, h2d_pt_etaPhiP, h2d_pt1pt1_eta1Phi1Eta2Phi2PP);
  Rho1Rho1EtaPhi(h2d_pt_etaPhiM, h2d_pt_etaPhiM, h2d_pt1pt1_eta1Phi1Eta2Phi2MM);

  //---------------calculation of r2-----------------------
  CalculateR2(h2d_n2_eta1Phi1Eta2Phi2PM, h2d_n1n1_eta1Phi1Eta2Phi2PM,  h2d_r2_eta1Phi1Eta2Phi2PM);
  CalculateR2(h2d_n2_eta1Phi1Eta2Phi2PP, h2d_n1n1_eta1Phi1Eta2Phi2PP,  h2d_r2_eta1Phi1Eta2Phi2PP);
  CalculateR2(h2d_n2_eta1Phi1Eta2Phi2MM, h2d_n1n1_eta1Phi1Eta2Phi2MM , h2d_r2_eta1Phi1Eta2Phi2MM);


  //--------------calculation of DpTDpT---------------------
  calculateDptDpt(h2d_ptpt_eta1Phi1Eta2Phi2PM, h2d_ptn_eta1Phi1Eta2Phi2PM,
                  h2d_npt_eta1Phi1Eta2Phi2PM, h2d_n2_eta1Phi1Eta2Phi2PM,
                  h1d_n1_ptP,
                  h1d_n1_ptM,
                  kminpt,
                  kmaxpt,
                  h2d_p2DptDpt_eta1Phi1Eta2Phi2PM, h2d_DptDpt_eta1Phi1Eta2Phi2PM,
                  knetabins,  knphibins);

  calculateDptDpt(h2d_ptpt_eta1Phi1Eta2Phi2PP, h2d_ptn_eta1Phi1Eta2Phi2PP,
                  h2d_npt_eta1Phi1Eta2Phi2PP, h2d_n2_eta1Phi1Eta2Phi2PP,
                  h1d_n1_ptP,
                  h1d_n1_ptP,
                  kminpt,
                  kmaxpt,
                  h2d_p2DptDpt_eta1Phi1Eta2Phi2PP, h2d_DptDpt_eta1Phi1Eta2Phi2PP,
                  knetabins,  knphibins);

  calculateDptDpt(h2d_ptpt_eta1Phi1Eta2Phi2MM, h2d_ptn_eta1Phi1Eta2Phi2MM,
                  h2d_npt_eta1Phi1Eta2Phi2MM, h2d_n2_eta1Phi1Eta2Phi2MM,
                  h1d_n1_ptM,
                  h1d_n1_ptM,
                  kminpt,
                  kmaxpt,
                  h2d_p2DptDpt_eta1Phi1Eta2Phi2MM, h2d_DptDpt_eta1Phi1Eta2Phi2MM,
                  knetabins,  knphibins);


  //reduce  etaphi to deltaEta deltaPhi
  //---------------R2-----------------------------
  EtaPhiTodEtadPhi(h2d_n1n1_eta1Phi1Eta2Phi2PM, h2d_n1n1_DetaDphiPM,
                   knetabins,  knphibins);

  EtaPhiTodEtadPhi(h2d_n1n1_eta1Phi1Eta2Phi2PP, h2d_n1n1_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_n1n1_eta1Phi1Eta2Phi2MM, h2d_n1n1_DetaDphiMM,
                   knetabins,  knphibins);

  EtaPhiTodEtadPhi(h2d_n2_eta1Phi1Eta2Phi2PM, h2d_n2_DetaDphiPM,
                   knetabins,  knphibins);

  EtaPhiTodEtadPhi(h2d_n2_eta1Phi1Eta2Phi2PP, h2d_n2_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_n2_eta1Phi1Eta2Phi2MM, h2d_n2_DetaDphiMM,
                   knetabins,  knphibins);

  EtaPhiTodEtadPhi(h2d_r2_eta1Phi1Eta2Phi2PM, h2d_r2_DetaDphiPM,
                   knetabins,  knphibins);

  EtaPhiTodEtadPhi(h2d_r2_eta1Phi1Eta2Phi2PP, h2d_r2_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_r2_eta1Phi1Eta2Phi2MM, h2d_r2_DetaDphiMM,
                   knetabins,  knphibins);

  //------------------DpTDpT-----------------------
  EtaPhiTodEtadPhi(h2d_ptpt_eta1Phi1Eta2Phi2PM, h2d_ptpt_DetaDphiPM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_ptpt_eta1Phi1Eta2Phi2PP, h2d_ptpt_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_ptpt_eta1Phi1Eta2Phi2MM, h2d_ptpt_DetaDphiMM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_ptn_eta1Phi1Eta2Phi2PM, h2d_ptn_DetaDphiPM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_ptn_eta1Phi1Eta2Phi2PP, h2d_ptn_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_ptn_eta1Phi1Eta2Phi2MM, h2d_ptn_DetaDphiMM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_npt_eta1Phi1Eta2Phi2PM, h2d_npt_DetaDphiPM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_npt_eta1Phi1Eta2Phi2PP, h2d_npt_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_npt_eta1Phi1Eta2Phi2MM, h2d_npt_DetaDphiMM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_pt1pt1_eta1Phi1Eta2Phi2PM, h2d_pt1pt1_DetaDphiPM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_pt1pt1_eta1Phi1Eta2Phi2PP, h2d_pt1pt1_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_pt1pt1_eta1Phi1Eta2Phi2MM, h2d_pt1pt1_DetaDphiMM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_DptDpt_eta1Phi1Eta2Phi2PM, h2d_DptDpt_DetaDphiPM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_DptDpt_eta1Phi1Eta2Phi2PP, h2d_DptDpt_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_DptDpt_eta1Phi1Eta2Phi2MM, h2d_DptDpt_DetaDphiMM,
                   knetabins,  knphibins);

  EtaPhiTodEtadPhi(h2d_p2DptDpt_eta1Phi1Eta2Phi2PM, h2d_p2DptDpt_DetaDphiPM,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_p2DptDpt_eta1Phi1Eta2Phi2PP, h2d_p2DptDpt_DetaDphiPP,
                   knetabins,  knphibins);
  EtaPhiTodEtadPhi(h2d_p2DptDpt_eta1Phi1Eta2Phi2MM, h2d_p2DptDpt_DetaDphiMM,
                   knetabins,  knphibins);


  //Symmetrization starts
  //--------------R2------------------------
  SymmetrizeDEtaDPhi(h2d_n1n1_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_n1n1_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_n1n1_DetaDphiMM);
  SymmetrizeDEtaDPhi(h2d_n2_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_n2_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_n2_DetaDphiMM);
  SymmetrizeDEtaDPhi(h2d_r2_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_r2_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_r2_DetaDphiMM);

  //---------------DpTDpT-------------------
  SymmetrizeDEtaDPhi(h2d_ptpt_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_ptpt_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_ptpt_DetaDphiMM);
  SymmetrizeDEtaDPhi(h2d_ptn_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_ptn_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_ptn_DetaDphiMM);
  SymmetrizeDEtaDPhi(h2d_npt_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_npt_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_npt_DetaDphiMM);
  SymmetrizeDEtaDPhi(h2d_pt1pt1_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_pt1pt1_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_pt1pt1_DetaDphiMM);
  SymmetrizeDEtaDPhi(h2d_DptDpt_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_DptDpt_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_DptDpt_DetaDphiMM);
  SymmetrizeDEtaDPhi(h2d_p2DptDpt_DetaDphiPM);
  SymmetrizeDEtaDPhi(h2d_p2DptDpt_DetaDphiPP);
  SymmetrizeDEtaDPhi(h2d_p2DptDpt_DetaDphiMM);

  //Shiftting starts
  //----------------R2---------------------------------
  shiftY(*h2d_r2_DetaDphiPM,   *h2d_r2_DetaDphiPM_ShiftY, nBins_ShiftY);
  shiftY(*h2d_r2_DetaDphiPP,   *h2d_r2_DetaDphiPP_ShiftY, nBins_ShiftY);
  shiftY(*h2d_r2_DetaDphiMM,   *h2d_r2_DetaDphiMM_ShiftY, nBins_ShiftY);

  //addition of PP & MM to get LS
  h2d_r2_DetaDphiLS->Add(h2d_r2_DetaDphiPP_ShiftY, h2d_r2_DetaDphiMM_ShiftY, .5, .5);

  //-------------------DpTDpT----------------------

  shiftY(*h2d_DptDpt_DetaDphiPM,   *h2d_DptDpt_DetaDphiPM_ShiftY, nBins_ShiftY);
  shiftY(*h2d_DptDpt_DetaDphiPP,   *h2d_DptDpt_DetaDphiPP_ShiftY, nBins_ShiftY);
  shiftY(*h2d_DptDpt_DetaDphiMM,   *h2d_DptDpt_DetaDphiMM_ShiftY, nBins_ShiftY);
  shiftY(*h2d_p2DptDpt_DetaDphiPM,   *h2d_p2DptDpt_DetaDphiPM_ShiftY, nBins_ShiftY);
  shiftY(*h2d_p2DptDpt_DetaDphiPP,   *h2d_p2DptDpt_DetaDphiPP_ShiftY, nBins_ShiftY);
  shiftY(*h2d_p2DptDpt_DetaDphiMM,   *h2d_p2DptDpt_DetaDphiMM_ShiftY, nBins_ShiftY);



  //addition of PP & MM to get LS
  h2d_p2DptDpt_DetaDphiLS->Add(h2d_p2DptDpt_DetaDphiPP_ShiftY, h2d_p2DptDpt_DetaDphiMM_ShiftY, .5, .5);

  //-------------calculate CI  -------------------
  //R2ci = (2*R2pm + R2pp + R2mm)/4.; or (0.5*R2pm + 0.25*R2pp + 0.25*R2mm)
  //-----------R2------------------------
  calculate_CI_CD(h2d_r2_DetaDphiPM_ShiftY, h2d_r2_DetaDphiPP_ShiftY, h2d_r2_DetaDphiMM_ShiftY, h2d_r2CI_DetaDphi, h2d_r2CD_DetaDphi);

  //---------------DpTDpT----------------

  calculate_CI_CD(h2d_DptDpt_DetaDphiPM_ShiftY, h2d_DptDpt_DetaDphiPP_ShiftY,
                  h2d_DptDpt_DetaDphiMM_ShiftY, h2d_DptDptCI_DetaDphi, h2d_DptDptCD_DetaDphi);
  //---------------p2DpTDpT----------------

  calculate_CI_CD(h2d_p2DptDpt_DetaDphiPM_ShiftY, h2d_p2DptDpt_DetaDphiPP_ShiftY,
                  h2d_p2DptDpt_DetaDphiMM_ShiftY, h2d_p2DptDptCI_DetaDphi, h2d_p2DptDptCD_DetaDphi);

smoothFunXRebinY(h2d_p2DptDptCD_DetaDphi, h2d_p2DptDptCD_DetaDphiRebined);

  //-------------calculate B2 -------------------
  //-----------R2------------------------
  //cout << "h2d_B2_DetaDphi Starts " << endl;
  calculate_B2(h2d_n1_etaPhiP, h2d_n1_etaPhiM,
               h2d_r2CD_DetaDphi, h2d_B2_DetaDphi);
  //cout << "h2d_B2_DetaDphi Ends " << endl;
  //---------------DpTDpT----------------
  calculate_B2(h2d_pt_etaPhiP, h2d_pt_etaPhiM,
               h2d_DptDptCD_DetaDphi, h2d_B2DptDpt_DetaDphi);

  //---------------p2DpTDpT----------------
  calculate_B2(h2d_pt_etaPhiP, h2d_pt_etaPhiM,
               h2d_p2DptDptCD_DetaDphi, h2d_B2p2DptDpt_DetaDphi);


  // create root file
  TFile *f = new TFile("R2P2B2_prc_.root"/*+strMult+"_"+strPt+".root"*/, "RECREATE");
f->cd();
 TH2D *h2d[20] = {h2d_r2_DetaDphiPM_ShiftY, h2d_r2_DetaDphiPP_ShiftY, h2d_r2_DetaDphiMM_ShiftY, h2d_r2_DetaDphiLS, h2d_r2CI_DetaDphi, h2d_r2CD_DetaDphi, h2d_B2_DetaDphi,  h2d_DptDpt_DetaDphiPM_ShiftY, h2d_DptDpt_DetaDphiPP_ShiftY, h2d_DptDpt_DetaDphiMM_ShiftY, h2d_DptDptCI_DetaDphi, h2d_DptDptCD_DetaDphi,  h2d_p2DptDpt_DetaDphiPM_ShiftY, h2d_p2DptDpt_DetaDphiPP_ShiftY, h2d_p2DptDpt_DetaDphiMM_ShiftY, h2d_p2DptDpt_DetaDphiLS, h2d_p2DptDptCI_DetaDphi, h2d_p2DptDptCD_DetaDphi, h2d_p2DptDptCD_DetaDphiRebined,h2d_B2p2DptDpt_DetaDphi};

  for (Int_t iHist2d = 0; iHist2d < 20; ++iHist2d)
  {
    h2d[iHist2d]->Write();
    delete h2d[iHist2d];
  }
  /*h1i_n1_multPM->Write();
    h1d_n1_etaPM->Write();
    h1d_n1_phiPM->Write();
    h1d_n1_ptP->Write();
    h1d_n1_etaP->Write();
    h1d_n1_phiP->Write();
    h1d_n1_ptM->Write();
    h1d_n1_etaM->Write();
    h1d_n1_phiM->Write();
    h1d_n1_ptPM_Tot->Write();
    //-----------R2------------------------
    h2d_n1_etaPhiP->Write();
    h2d_n1_etaPhiM->Write();

    h2d_n1n1_eta1Phi1Eta2Phi2PM->Write();
    h2d_n1n1_eta1Phi1Eta2Phi2PP->Write();
    h2d_n1n1_eta1Phi1Eta2Phi2MM->Write();
    h2d_n1n1_DetaDphiPM->Write();
    h2d_n1n1_DetaDphiPP->Write();
    h2d_n1n1_DetaDphiMM->Write();
    h2d_n2_eta1Phi1Eta2Phi2PM->Write();
    h2d_n2_eta1Phi1Eta2Phi2PP->Write();
    h2d_n2_eta1Phi1Eta2Phi2MM->Write();
    h2d_n2_DetaDphiPM->Write();
    h2d_n2_DetaDphiPP->Write();
    h2d_n2_DetaDphiMM->Write();
    h2d_r2_eta1Phi1Eta2Phi2PM->Write();
    h2d_r2_eta1Phi1Eta2Phi2PP->Write();
    h2d_r2_eta1Phi1Eta2Phi2MM->Write();
    h2d_n1n1_DetaDphiPM_ShiftY->Write();
    h2d_n1n1_DetaDphiPP_ShiftY->Write();
    h2d_n1n1_DetaDphiMM_ShiftY->Write();
    h2d_n2_DetaDphiPM_ShiftY->Write();
    h2d_n2_DetaDphiPP_ShiftY->Write();
    h2d_n2_DetaDphiMM_ShiftY->Write();
    h2d_r2_DetaDphiPP_ShiftY->Write();
    h2d_r2_DetaDphiMM_ShiftY->Write();*/

  //---------------DpTDpT----------------
  /*h2d_pt_etaPhiP->Write();
    h2d_pt_etaPhiM->Write();

    h2d_ptpt_eta1Phi1Eta2Phi2PM->Write();
    h2d_ptpt_eta1Phi1Eta2Phi2PP->Write();
    h2d_ptpt_eta1Phi1Eta2Phi2MM->Write();
    h2d_ptn_eta1Phi1Eta2Phi2PM->Write();
    h2d_ptn_eta1Phi1Eta2Phi2PP->Write();
    h2d_ptn_eta1Phi1Eta2Phi2MM->Write();
    h2d_npt_eta1Phi1Eta2Phi2PM->Write();
    h2d_npt_eta1Phi1Eta2Phi2PP->Write();
    h2d_npt_eta1Phi1Eta2Phi2MM->Write();

    h2d_pt1pt1_eta1Phi1Eta2Phi2PM->Write();
    h2d_pt1pt1_eta1Phi1Eta2Phi2PP->Write();
    h2d_pt1pt1_eta1Phi1Eta2Phi2MM->Write();
    h2d_ptn_eta1Phi1Eta2Phi2PM->Write();
    h2d_ptn_eta1Phi1Eta2Phi2PP->Write();
    h2d_ptn_eta1Phi1Eta2Phi2MM->Write();
    h2d_npt_eta1Phi1Eta2Phi2PM->Write();
    h2d_npt_eta1Phi1Eta2Phi2PP->Write();
    h2d_npt_eta1Phi1Eta2Phi2MM->Write();
    h2d_DptDpt_eta1Phi1Eta2Phi2PM->Write();
    h2d_DptDpt_eta1Phi1Eta2Phi2PP->Write();
    h2d_DptDpt_eta1Phi1Eta2Phi2MM->Write();
    h2d_p2DptDpt_eta1Phi1Eta2Phi2PM->Write();
    h2d_p2DptDpt_eta1Phi1Eta2Phi2PP->Write();
    h2d_p2DptDpt_eta1Phi1Eta2Phi2MM->Write();
    h2d_pt1pt1_DetaDphiPM->Write();
    h2d_pt1pt1_DetaDphiPP->Write();
    h2d_pt1pt1_DetaDphiMM->Write();
    h2d_ptn_DetaDphiPM->Write();
    h2d_ptn_DetaDphiPP->Write();
    h2d_ptn_DetaDphiMM->Write();
    h2d_npt_DetaDphiPM->Write();
    h2d_npt_DetaDphiPP->Write();
    h2d_npt_DetaDphiMM->Write();
    h2d_DptDpt_DetaDphiPM->Write();
    h2d_DptDpt_DetaDphiPP->Write();
    h2d_DptDpt_DetaDphiMM->Write();
    h2d_DptDpt_DetaDphiPM_ShiftY->Write();
    h2d_DptDpt_DetaDphiPP_ShiftY->Write();
    h2d_DptDpt_DetaDphiMM_ShiftY->Write();
    h2d_DptDptCI_DetaDphi->Write();
    h2d_DptDptCD_DetaDphi->Write();
    h2d_B2DptDpt_DetaDphi->Write();

    h2d_p2DptDpt_DetaDphiPM->Write();
    h2d_p2DptDpt_DetaDphiPP->Write();
    h2d_p2DptDpt_DetaDphiMM->Write();
    h2d_p2DptDpt_DetaDphiPP_ShiftY->Write();
    h2d_p2DptDpt_DetaDphiMM_ShiftY->Write();*/

  //-------------pt1pt2dphi Corr. ---------------------
  /*h2d_n1_ptPhiP->Write();
    h2d_n1_ptPhiM->Write();
    h3d_n2_pt1Pt2DphiPM->Write();
    h3d_n2_pt1Pt2DphiPP->Write();
    h3d_n2_pt1Pt2DphiMM->Write();

    h3d_n1n1_pt1Pt2DphiPM->Write();
    h3d_n1n1_pt1Pt2DphiPP->Write();
    h3d_n1n1_pt1Pt2DphiMM->Write();
    h3d_r2_pt1Pt2DphiPM->Write();
    h3d_r2_pt1Pt2DphiPP->Write();
    h3d_r2_pt1Pt2DphiMM->Write();
    h2d_r2_pt1Pt2PM_NearSide->Write();
    h2d_r2_pt1Pt2PP_NearSide->Write();
    h2d_r2_pt1Pt2MM_NearSide->Write();
    h2d_r2_pt1Pt2PM_AwaySide->Write();
    h2d_r2_pt1Pt2PP_AwaySide->Write();
    h2d_r2_pt1Pt2MM_AwaySide->Write();
    h2d_r2_pt1Pt2PM_NearSide_Symm->Write();
    h2d_r2_pt1Pt2PP_NearSide_Symm->Write();
    h2d_r2_pt1Pt2MM_NearSide_Symm->Write();
    h2d_r2_pt1Pt2PM_AwaySide_Symm->Write();
    h2d_r2_pt1Pt2PP_AwaySide_Symm->Write();
    h2d_r2_pt1Pt2MM_AwaySide_Symm->Write();
    h2d_r2CI_pt1Pt2_NearSide->Write();
    h2d_r2CI_pt1Pt2_AwaySide->Write();
    h2d_r2CD_pt1Pt2_NearSide->Write();
    h2d_r2CD_pt1Pt2_AwaySide->Write();
    h2d_B2_pt1Pt2_NearSide->Write();
    h2d_B2_pt1Pt2_AwaySide->Write();
  */
  f->Write();
  delete f;

  /*delete h1i_n1_multPM;
  delete h1d_n1_etaPM;
  delete h1d_n1_phiPM;
  delete h1d_n1_ptP;
  delete h1d_n1_etaP;
  delete h1d_n1_phiP;
  delete h1d_n1_ptM;
  delete h1d_n1_etaM;
  delete h1d_n1_phiM;
  //delete h1d_n1_ptPM_Tot;
  //-----------R2------------------------
  delete h2d_n1_etaPhiP;
  delete h2d_n1_etaPhiM;
  delete h2d_n1n1_eta1Phi1Eta2Phi2PM;
  delete h2d_n1n1_eta1Phi1Eta2Phi2PP;
  delete h2d_n1n1_eta1Phi1Eta2Phi2MM;
  delete h2d_n1n1_DetaDphiPM;
  delete h2d_n2_eta1Phi1Eta2Phi2PM;
  delete h2d_n2_DetaDphiPM;
  delete h2d_r2_eta1Phi1Eta2Phi2PM;
  delete h2d_n1n1_DetaDphiPM_ShiftY;
  delete h2d_n1n1_DetaDphiPP_ShiftY;
  delete h2d_n1n1_DetaDphiMM_ShiftY;
  delete h2d_n2_DetaDphiPM_ShiftY;
  delete h2d_n2_DetaDphiPP_ShiftY;
  delete h2d_n2_DetaDphiMM_ShiftY;
  delete h2d_r2_DetaDphiPM_ShiftY;
  delete h2d_r2_DetaDphiPP_ShiftY;
  delete h2d_r2_DetaDphiMM_ShiftY;
  delete h2d_r2_DetaDphiLS;
  delete h2d_r2CI_DetaDphi;
  delete h2d_r2CD_DetaDphi;
  delete h2d_B2_DetaDphi;
  //hdPhi;
  //hdEta;

  //---------------DpTDpT----------------
  delete h2d_pt_etaPhiP;
  delete h2d_pt_etaPhiM;
  delete h2d_ptpt_eta1Phi1Eta2Phi2PM;
  delete h2d_ptpt_eta1Phi1Eta2Phi2PP;
  delete h2d_ptpt_eta1Phi1Eta2Phi2MM;
  delete h2d_ptn_eta1Phi1Eta2Phi2PM;
  delete h2d_ptn_eta1Phi1Eta2Phi2PP;
  delete h2d_ptn_eta1Phi1Eta2Phi2MM;
  delete h2d_npt_eta1Phi1Eta2Phi2PM;
  delete h2d_npt_eta1Phi1Eta2Phi2PP;
  delete h2d_npt_eta1Phi1Eta2Phi2MM;

  delete h2d_pt1pt1_eta1Phi1Eta2Phi2PM;
  delete h2d_pt1pt1_eta1Phi1Eta2Phi2PP;
  delete h2d_pt1pt1_eta1Phi1Eta2Phi2MM;
  delete h2d_ptn_eta1Phi1Eta2Phi2PM;
  delete h2d_ptn_eta1Phi1Eta2Phi2PP;
  delete h2d_ptn_eta1Phi1Eta2Phi2MM;
  delete h2d_npt_eta1Phi1Eta2Phi2PM;
  delete h2d_npt_eta1Phi1Eta2Phi2PP;
  delete h2d_npt_eta1Phi1Eta2Phi2MM;
  delete h2d_DptDpt_eta1Phi1Eta2Phi2PM;
  delete h2d_DptDpt_eta1Phi1Eta2Phi2PP;
  delete h2d_DptDpt_eta1Phi1Eta2Phi2MM;
  delete h2d_p2DptDpt_eta1Phi1Eta2Phi2PM;
  delete h2d_p2DptDpt_eta1Phi1Eta2Phi2PP;
  delete h2d_p2DptDpt_eta1Phi1Eta2Phi2MM;
  delete h2d_pt1pt1_DetaDphiPM;
  delete h2d_pt1pt1_DetaDphiPP;
  delete h2d_pt1pt1_DetaDphiMM;
  delete h2d_ptn_DetaDphiPM;
  delete h2d_ptn_DetaDphiPP;
  delete h2d_ptn_DetaDphiMM;
  delete h2d_npt_DetaDphiPM;
  delete h2d_npt_DetaDphiPP;
  delete h2d_npt_DetaDphiMM;
  delete h2d_DptDpt_DetaDphiPM;
  delete h2d_DptDpt_DetaDphiPP;
  delete h2d_DptDpt_DetaDphiMM;
  delete h2d_DptDpt_DetaDphiPM_ShiftY;
  delete h2d_DptDpt_DetaDphiPP_ShiftY;
  delete h2d_DptDpt_DetaDphiMM_ShiftY;
  delete h2d_p2DptDpt_DetaDphiPM;
  delete h2d_p2DptDpt_DetaDphiPP;
  delete h2d_p2DptDpt_DetaDphiMM;
  delete h2d_p2DptDpt_DetaDphiPM_ShiftY;
  delete h2d_p2DptDpt_DetaDphiPP_ShiftY;
  delete h2d_p2DptDpt_DetaDphiMM_ShiftY;
  delete h2d_p2DptDpt_DetaDphiLS;
  delete h2d_DptDptCI_DetaDphi;
  delete h2d_DptDptCD_DetaDphi;
  delete h2d_B2DptDpt_DetaDphi;
  delete h2d_p2DptDptCI_DetaDphi;
  delete h2d_p2DptDptCD_DetaDphi;
  delete h2d_B2p2DptDpt_DetaDphi;*/
  cout <<"\n ========== R2P2B2_prc_v3 ends ======= \n" <<endl;
}

//---------Methods Starts-------------------------------

//------Calculate External Product n1_1 x n1_2 and store into n1n1_12-----
Double_t Rho1Rho1EtaPhi(const TH2D *h_1, const TH2D * h_2, TH2D * h_12)
{
  // This function calculates dN/(deta1*dphi1*deta2*dphi2)
  if (!h_1  || !h_2 || !h_12)
  {
    cout << "-DEBUG- -E- calculateN1N1_H2H2H2(...) Null pointers as arguments" << endl;
    cout << "-DEBUG- ABORT!!!!!" << endl;
    return -1.;
  }
  Int_t n1x = h_1->GetNbinsX();
  Int_t n1y = h_1->GetNbinsY();
  Int_t n2x = h_2->GetNbinsX();
  Int_t n2y = h_2->GetNbinsY();
  Int_t n3x = h_12->GetNbinsX();
  Int_t n3y = h_12->GetNbinsY();
  if (n3x != (n1x * n1y) || n3y != (n2x * n2y) )
  {
    cout << "-DEBUG- -E- calculateN1N1_H2H2H2(...) Incompatible histo dimensions" << endl;
    cout << "-DEBUG- H1: " << h_1->GetName()  << " nBins_x:" << n1x << " nBins_y:" << n1y << endl;
    cout << "-DEBUG- H2: " << h_2->GetName()  << " nBins_x:" << n2x << " nBins_y:" << n2y << endl;
    cout << "-DEBUG- H3: " << h_12->GetName() << " nBins_x:" << n3x << " nBins_y:" << n3y << endl;
    cout << "-DEBUG- ABORT!!!!!" << endl;
    return -1.;
  }

  Double_t v1, ev1, v2, ev2, v, ev, r1, r2;
  Double_t sum  = 0.;
  Double_t norm = 0.;
  for (Int_t i1x = 0; i1x < n1x; ++i1x)
  {
    for (Int_t i1y = 0; i1y < n1y; ++i1y)
    {
      v1  = h_1->GetBinContent(i1x + 1, i1y + 1);
      ev1 = h_1->GetBinError(i1x + 1, i1y + 1);
      for (Int_t i2x = 0; i2x < n2x; ++i2x)
      {
        for (Int_t i2y = 0; i2y < n2y; ++i2y)
        {
          v2  = h_2->GetBinContent(i2x + 1, i2y + 1);
          ev2 = h_2->GetBinError(i2x + 1, i2y + 1);
          v = v1 * v2;
          if (v > 0)
          {
            r1 = ev1 / v1;
            r2 = ev2 / v2;
            ev = v * TMath::Sqrt(r1 * r1 + r2 * r2);
          }
          else
          {
            v = 0.;
            ev = 0.;
          }
          Int_t i3x = i1x * n1y + i1y + 1;
          Int_t i3y = i2x * n2y + i2y + 1;
          h_12->SetBinContent(i3x, i3y, v);
          h_12->SetBinError(i3x, i3y, ev);
          sum += v;
          norm = 1.;
        }
      }
    }
  }
  //return average across bins
  return sum / norm;
}

//---------Calculate R2 = N2/N1/N1 - 1--------------------
void CalculateR2(const TH2D * n2_12, const TH2D * n1n1_12, TH2D * r2_12)
{
  //cout << "r2DetaDphi Starts :" << r2_12->GetName() << endl;
  if (!n2_12  || !n1n1_12 || !r2_12)
  {
    cout << "-DEBUG- -E- calculateR2_H2H2H2(...) Null pointers as arguments" << endl;
    cout << "-DEBUG- ABORT!!!!!" << endl;
    return;
  }
  Int_t n2_12_n_x    = n2_12->GetNbinsX();
  Int_t n2_12_n_y    = n2_12->GetNbinsY();
  Int_t n1n1_12_n_x  = n1n1_12->GetNbinsX();
  Int_t n1n1_12_n_y  = n1n1_12->GetNbinsY();
  Int_t r2_12_n_x    = r2_12->GetNbinsX();
  Int_t r2_12_n_y    = r2_12->GetNbinsY();


  if (n2_12_n_x != n1n1_12_n_x || n2_12_n_x != r2_12_n_x || n2_12_n_y != n1n1_12_n_y || n2_12_n_y != r2_12_n_y)
  {
    cout << "-DEBUG- -E- calculateR2_H2H2H2(...) Incompatible histo dimensions" << endl;
    cout << "-DEBUG- H1: " << n2_12->GetName()   << " n_x:" << n2_12_n_x   << " n_y:" << n2_12_n_y   << endl;
    cout << "-DEBUG- H2: " << n1n1_12->GetName() << " n_x:" << n1n1_12_n_x << " n_y:" << n1n1_12_n_y << endl;
    cout << "-DEBUG- H3: " << r2_12->GetName()   << " n_x:" << r2_12_n_x   << " n_y:" << r2_12_n_y   << endl;
    cout << "-DEBUG- ABORT!!!!!" << endl;
    return;
  }

  // cout << "*****  STEP 2 -DEBUG- H1: " << n2_12->GetName()   << " n_x:" << n2_12_n_x   << " n_y:" << n2_12_n_y   << endl;



  Double_t v1, ev1, v2, ev2, v, ev, re1, re2;
  Int_t i_x = 1;
  Int_t i_y = 1;

  for ( i_x = 1; i_x <= n2_12_n_x; ++i_x)
  {
    for ( i_y = 1; i_y <= n2_12_n_y; ++i_y)
    {
      v1  = n2_12->GetBinContent(i_x, i_y);
      ev1 = n2_12->GetBinError(i_x, i_y);
      v2  = n1n1_12->GetBinContent(i_x, i_y);
      ev2 = n1n1_12->GetBinError(i_x, i_y);
      // cout << " *** INSIDE **** " << i_x << " " << i_y << " " << v1 << " " << v2 << endl;
      if (v1 > 0 && v2 > 0) //   && ev1/v1<0.5  && ev2/v2<0.5)
      {
        v   = v1 / v2; // all pairs counted - no need to multiply by 2
        re1 = ev1 / v1;
        re2 = ev2 / v2;
        ev  = v * TMath::Sqrt(re1 * re1 + re2 * re2);//std deviation
        v   -= 1.;
        // cout << " Hello v :" << v << endl;
      }
      //cout << i_x << "\t"<< i_y << "\t"<<  v << endl;

      else
      {
        v = 0.;
        ev = 0;
      }
      r2_12->SetBinContent(i_x, i_y, v);
      r2_12->SetBinError(i_x, i_y, ev);


    }

  }

  //cout << "r2DetaDphi ends :" << r2_12->GetName() << endl;
}

void EtaPhiTodEtadPhi(const TH2D * source, TH2D * target, Int_t nEtaBins, Int_t nPhiBins)
{
  //
  // Input is eta phi
  // Output  is deltaeta deltaphi
  //
  // cout << "-DEBUG- reduce_n2xEtaPhi_n2DetaDphi() ==============  New Version From TH2D" << endl;
  Double_t v1, v2, ev1;
  Int_t dPhi, dEta, iPhi, iEta, jPhi, jEta, i, j;
  Int_t nBins = nEtaBins * nPhiBins;
  Int_t nWrk  = nPhiBins * (2 * nEtaBins - 1);
  Int_t index;
  Double_t * numerator    = new Double_t[nWrk];
  Double_t * numeratorErr = new Double_t[nWrk];
  Double_t * denominator  = new Double_t[nWrk];
  for (Int_t k = 0; k < nWrk; ++k)
  {
    numerator[k]    = 0;
    numeratorErr[k] = 0;
    denominator[k]  = 0;
  }

  TString name = target->GetName();

  i = 1;
  for (iEta = 0; iEta < nEtaBins; ++iEta)
  {
    for (iPhi = 0; iPhi < nPhiBins; ++iPhi)
    {
      j = 1;
      for (jEta = 0; jEta < nEtaBins; ++jEta)
      {
        for (jPhi = 0; jPhi < nPhiBins; ++jPhi)
        {
          dPhi = iPhi - jPhi; if (dPhi < 0) dPhi += nPhiBins; dPhi += 1;
          dEta = iEta - jEta + nEtaBins;
          v1   = source->GetBinContent(i, j);
          ev1  = source->GetBinError(  i, j);
          index = (dEta - 1) * nPhiBins + dPhi - 1;
          numerator[index]    += v1;
          numeratorErr[index] += ev1 * ev1;
          denominator[index]  += 1.;
          //cout << "-DEBUG-  " << name << "  iEta:" << iEta << " iPhi:" << iPhi << " jEta:" << jEta << " jPhi:" << jPhi << " v1:" << v1 << " ev1:" << ev1 << endl;
          ++j;
        }
      }
      ++i;
    }
  }
  for (dEta = 0; dEta < 2 * nEtaBins - 1; ++dEta)
  {
    for (dPhi = 0; dPhi < nPhiBins; ++dPhi)
    {
      v1   = target->GetBinContent(dEta + 1, dPhi + 1);
      ev1  = target->GetBinError(dEta + 1, dPhi + 1);
      index = dEta * nPhiBins + dPhi;
      v1    = numerator[index];
      ev1   = numeratorErr[index];
      v2    = denominator[index];
      if (v2 <= 0)  cout << "-DEBUG- HistogramCollection::reduce_n2xEtaPhi_n2DetaDphi() Elements of denominator are negative." << endl;
      target->SetBinContent(dEta + 1, dPhi + 1, v1 / v2);
      target->SetBinError(  dEta + 1, dPhi + 1, TMath::Sqrt(ev1) / v2);
    }
  }

  delete [] numerator;
  delete [] numeratorErr;
  delete [] denominator;
}

// h must be writtable otherwise this does not work...
// x axis is DeltaEta. It must have an odd number of bins
// y axis is DeltaPhi. It must have an even number of bins

void SymmetrizeDEtaDPhi(TH2D *h)
{

  Double_t v1, v2, v3, v4;
  Double_t ev1, ev2, ev3, ev4;
  Double_t sv, esv;
  Int_t nEta = h->GetNbinsX(); //DeltaEta
  Int_t nPhi = h->GetNbinsY(); //DeltaPhi

  Int_t nEtaHalf = (nEta - 1) / 2.;
  Int_t nPhiHalf = (nPhi - 2) / 2.;
  Int_t iEta, iPhi, iPhi1, iEta1;
  Double_t * v = new Double_t[nEta * nPhi];
  Double_t * ev = new Double_t[nEta * nPhi];
  // cout << "-DEBUG- symmetrizeDeltaEtaDeltaPhi(TH2D * h) Arrays created" << endl;
  for (Int_t iPhi = 0; iPhi < nPhi; iPhi++)
  {
    for (Int_t iEta = 0; iEta < nEta; iEta++)
    {
      iPhi1 = iPhi + 1;
      iEta1 = iEta + 1;
      v[ iEta + iPhi * nEta]  = h->GetBinContent(iEta1, iPhi1);
      ev[iEta + iPhi * nEta]  = h->GetBinError(  iEta1, iPhi1);
    }
  }
  // cout << "-DEBUG- symmetrizeDeltaEtaDeltaPhi(TH2D * h) Arrays copied" << endl;
  for (iEta = 0; iEta < nEtaHalf; iEta++)
  {
    iEta1 = iEta + 1;
    for (iPhi = 0; iPhi < nPhiHalf; iPhi++)
    {
      iPhi1 = iPhi + 1;
      v1 = v[  nEta - iEta1 + (nPhi - iPhi1) * nEta];
      v2 = v[  nEta - iEta1 + (     iPhi1) * nEta];
      v3 = v[        iEta + (nPhi - iPhi1) * nEta];
      v4 = v[        iEta + (     iPhi1) * nEta];
      ev1 = ev[nEta - iEta1 + (nPhi - iPhi1) * nEta];
      ev2 = ev[nEta - iEta1 + (     iPhi1) * nEta];
      ev3 = ev[      iEta + (nPhi - iPhi1) * nEta];
      ev4 = ev[      iEta + (     iPhi1) * nEta];

      sv = (v1 + v2 + v3 + v4) / 4.;
      esv = TMath::Sqrt(ev1 * ev1 + ev2 * ev2 + ev3 * ev3 + ev4 * ev4) / 4.;

      h->SetBinContent( nEta - iEta, nPhi - iPhi, sv);
      h->SetBinContent( nEta - iEta,   iPhi1 + 1, sv);
      h->SetBinContent(     iEta1, nPhi - iPhi, sv);
      h->SetBinContent(     iEta1,   iPhi1 + 1, sv);
      h->SetBinError(   nEta - iEta, nPhi - iPhi, esv);
      h->SetBinError(   nEta - iEta,   iPhi1 + 1, esv);
      h->SetBinError(       iEta1, nPhi - iPhi, esv);
      h->SetBinError(       iEta1,   iPhi1 + 1, esv);
    }
  }
  // cout << "-DEBUG- symmetrizeDeltaEtaDeltaPhi(TH2D * h) Part 1 Done" << endl;
  iEta  = nEtaHalf;
  iEta1 = iEta + 1;
  for (iPhi = 0; iPhi < nPhiHalf; iPhi++) // iEta center bin
  {
    iPhi1 = iPhi + 1;
    v3 = v[        iEta + (nPhi - iPhi1) * nEta];
    v4 = v[        iEta + (     iPhi1) * nEta];
    ev3 = ev[      iEta + (nPhi - iPhi1) * nEta];
    ev4 = ev[      iEta + (     iPhi1) * nEta];

    sv = (v3 + v4) / 2.;
    esv = TMath::Sqrt(ev3 * ev3 + ev4 * ev4) / 2.;

    h->SetBinContent(     iEta1, nPhi - iPhi, sv);
    h->SetBinContent(     iEta1,   iPhi1 + 1, sv);
    h->SetBinError(       iEta1, nPhi - iPhi, esv);
    h->SetBinError(       iEta1,   iPhi1 + 1, esv);
  }
  // cout << "-DEBUG- symmetrizeDeltaEtaDeltaPhi(TH2D * h) Part 2 Done" << endl;

  iPhi = 0;
  iPhi1 = iPhi + 1;
  for (iEta = 0; iEta < nEtaHalf; iEta++)
  {
    iEta1 = iEta + 1;
    v1 = v[  nEta - iEta1];
    v3 = v[        iEta];
    ev1 = ev[nEta - iEta1];
    ev3 = ev[      iEta];

    sv = (v1 + v3) / 2.;
    esv = TMath::Sqrt(ev1 * ev1 + ev3 * ev3) / 2.;

    h->SetBinContent( nEta - iEta, 1, sv);
    h->SetBinContent(     iEta1, 1, sv);
    h->SetBinError(   nEta - iEta, 1, esv);
    h->SetBinError(       iEta1, 1, esv);

    iPhi  = nPhi / 2;
    iPhi1 = iPhi + 1;
    v1 = v[  nEta - iEta1 + iPhi * nEta];
    v3 = v[        iEta + iPhi * nEta];
    ev1 = ev[nEta - iEta1 + iPhi * nEta];
    ev3 = ev[      iEta + iPhi * nEta];

    sv = (v1 + v3) / 2.;
    esv = TMath::Sqrt(ev1 * ev1 + ev3 * ev3) / 2.;

    h->SetBinContent( nEta - iEta, iPhi1, sv);
    h->SetBinContent(     iEta1, iPhi1, sv);
    h->SetBinError(   nEta - iEta, iPhi1, esv);
    h->SetBinError(       iEta1, iPhi1, esv);

  }
  delete[] v;

}


///shift the given source to the target vertically by nbins
void shiftY(const TH2D & source, TH2D & target, Int_t nbins)
{
  Int_t i_x, i_y;
  Int_t n_x = source.GetNbinsX();
  Int_t n_y = source.GetNbinsY();

  //cout << "nx" << n_x << "ny" << n_y << "nBins" << nbins << "n_y-nbins" << n_y-nbins << "n_y-nbins+1" <<n_y-nbins+1 << endl;
  //shift the 1st area
  for (i_x = 1; i_x <= n_x; ++i_x)
  {
    for (i_y = 1; i_y <= n_y - nbins; ++i_y)
    {
      Double_t v  = source.GetBinContent(i_x, i_y);
      Double_t ev = source.GetBinError(i_x, i_y);
      target.SetBinContent(i_x,  i_y + nbins, v);
      target.SetBinError(i_x,    i_y + nbins, ev);
      //cout << "i_y :" << i_y << "i_y + nbins :" << i_y + nbins << "\t";
    }
    for (i_y = n_y - nbins + 1; i_y <= n_y; ++i_y)
    {
      Double_t v  = source.GetBinContent(i_x, i_y);
      Double_t ev = source.GetBinError(i_x, i_y);
      target.SetBinContent(i_x, i_y - (n_y - nbins), v);
      target.SetBinError(i_x,   i_y - (n_y - nbins), ev);
      //cout << "i_y :" << i_y << "i_y - (n_y - nbins) :" << i_y - (n_y - nbins) << endl;
    }
  }
}

//--------------DpTDpT-------------------------------

bool sameDimensions2D(const TH2D* h1, const TH2D* h2)
{
  Int_t n1x = h1->GetNbinsX();
  Int_t n1y = h1->GetNbinsY();
  Int_t n2x = h2->GetNbinsX();
  Int_t n2y = h2->GetNbinsY();
  if (n1x == n2x && n1y == n2y)
  {
    return true;
  }
  else
  {
    cout << "-DEBUG- -ERROR- Histograms titled " << h1->GetTitle() << " and " << h2->GetTitle() << " have incompatible dimensions" << endl;
    return false;
  }

}

void calculateDptDpt(const TH2D * spp, const TH2D * spn, const TH2D * snp, const TH2D * snn, const TH1D * avgpt1, const TH1D * avgpt2, const Float_t kMinPt, const Float_t kMaxPt, TH2D * p2DptDpt,  TH2D * dptdpt, Int_t nEta, Int_t nPhi) {
  if (!sameDimensions2D(spp, spn)) return;
  if (!sameDimensions2D(spp, snp)) return;
  if (!sameDimensions2D(spp, snn)) return;
  //Int_t nx = spp->GetNbinsX();
  //Int_t ny = spp->GetNbinsY();

  Double_t v1, ev1, v2, ev2, v3, ev3, v4, ev4, v5, ev5;
  Double_t v6, ev6, v7, ev7, p1, p2;
  Int_t k, k1, k2;

  Int_t lowBinPt1 = avgpt1->GetXaxis()->FindBin(kMinPt +0.01);
  Int_t highBinPt1 = avgpt1->GetXaxis()->FindBin(kMaxPt - 0.01);

  /*  cout << "lowPt1: " << kMinPt << "\t"
       << "highPt1: " << kMaxPt <<  endl;

  cout << "lowBinPt1: " << lowBinPt1 << "\t"
       << "highBinPt1: " << highBinPt1 <<  endl;
  */
  v1 = 0.0, v2 = 0.0;
  Double_t  x1 = 0.0, sum1 = 0.0, xSum1 = 0.0;
  Double_t  x2 = 0.0, sum2 = 0.0, xSum2 = 0.0;

  for (Int_t ipt1 = lowBinPt1; ipt1 <= highBinPt1; ++ipt1)
  {
    v1 = avgpt1->GetBinContent(ipt1);
    x1 = avgpt1->GetXaxis()->GetBinCenter(ipt1);
    sum1     += v1;
    xSum1    += v1 * x1;

    v2 = avgpt2->GetBinContent(ipt1);
    x2 = avgpt2->GetXaxis()->GetBinCenter(ipt1);
    sum2     += v2;
    xSum2    += v2 * x2;
  }

  p1 = xSum1 / sum1;
  p2 = xSum2 / sum2;			//THese are avg momentum of particle 1 & 2...

  /*  cout << "What is my average pT?" << endl;
  cout << avgpt1->GetName() << ":: p1: " << p1 << endl;
  cout << avgpt2->GetName() << ":: p2: " << p2 << endl;
  */

  for (Int_t iEta1 = 0; iEta1 < nEta; ++iEta1)
  {
    //p1  = ptAvg1[iEta1];
    for (Int_t iPhi1 = 0; iPhi1 < nPhi; ++iPhi1)
    {
      k1 = 1 + iEta1 * nPhi + iPhi1;
      //c1 = avgp1->GetBinContent(k1+1)/p1;
      for (Int_t iEta2 = 0; iEta2 < nEta; ++iEta2)
      {
        //p2  = ptAvg2[iEta2];
        for (Int_t iPhi2 = 0; iPhi2 < nPhi; ++iPhi2)
        {
          k2 = 1 + iEta2 * nPhi + iPhi2;
          v1  = spp->GetBinContent(k1, k2); ev1 = spp->GetBinError(k1, k2);
          v2  = spn->GetBinContent(k1, k2); ev2 = spn->GetBinError(k1, k2);
          v3  = snp->GetBinContent(k1, k2); ev3 = snp->GetBinError(k1, k2);
          v4  = snn->GetBinContent(k1, k2); ev4 = snn->GetBinError(k1, k2);
          //c2 = avgp2->GetBinContent(k2+1)/p2;

          if (v4 != 0 ) // && ev4/v4<0.5)
          {
            v5  = v1 - v2 * p2 - p1 * v3 + p1 * p2 * v4;
            ev5 = v5 * ev4 / v4;
            v6  = v5 / v4;
            ev6 = v6 * ev4 / v4;

            v7 = (v6 / (p1 * p2));			//Ok so doing this does give the desired plot...although is it legal.
            ev7 = ev6 / (p1 * p2);
          }
          else
          {
            v5 = v6 = v7 = ev5 = ev6 = ev7 = 0;
          }

          dptdpt->SetBinContent(  k1, k2, v6);
          dptdpt->SetBinError(  k1, k2, ev6);

          p2DptDpt->SetBinContent(k1, k2, v7);
          p2DptDpt->SetBinError(k1, k2, ev7);
        }
      }
    }
  }

}


void calculate_CI_CD(const TH2D *hpm, const TH2D *hpp, const TH2D *hmm,  TH2D *hci,  TH2D *hcd)
{

  hci->Add(hpp, hmm, 0.25, 0.25);
  hci->Add(hpm, 0.5);

  hcd->Add(hpp, hmm, -0.25, -0.25);
  hcd->Add(hpm, 0.5);
}

///rebin by 2 in dphi-direction
//use smoothing function for deta-dir as it has odd-bin
//T(i) =  2S(i) + S(i+1)/3  for first bin
//T(i) = S(i-1) + 2S(i)/3   for last bin
//T(i) = S(i-1) + 2S(i) + S(i+1)/4 for  other bin except end-bin
//void smoothFunXRebinY(const TH2D *source, const TH2D *target){
void smoothFunXRebinY( TH2D *source,  TH2D *target){
  int i_x, i_y;
  int n_x = source->GetNbinsX();
  int n_y = source->GetNbinsY();
  double vLTnx = 0., evLTnx = 0., vGT1 = 0., evGT1 = 0.; 
  for (i_x = 1; i_x <= n_x; ++i_x)
  {
    for (i_y = 1; i_y <= n_y; ++i_y)
    {
//---------------Getter-----------------
      //for S(i) bin
      double v  = source->GetBinContent(i_x, i_y);
      double ev = source->GetBinError(i_x, i_y);

      //for S(i-1) bin
      if (i_x > 1 ) {
        //GT1 = greater than 1
        vGT1  = source->GetBinContent(i_x - 1 , i_y);
        evGT1 = source->GetBinError(i_x - 1, i_y);
      }
      //for S(i+1)
      if (i_x < n_x)
      {
        //LTnx = lesser than nx
        vLTnx  = source->GetBinContent(i_x + 1, i_y);
        evLTnx = source->GetBinError(i_x + 1, i_y);
      }

//--------------- Setter--------------

      //T(i) =  2S(i) + S(i+1)/3  for first bin

      if (i_x == 1) {

        v = (2.0 * v + vLTnx) / 3.0;
        ev = TMath::Sqrt((2.0 / 3.0) * (2.0 / 3.0) * ev * ev + (1.0 / 3.0) * (1.0 / 3.0) * evLTnx * evLTnx );

        target->SetBinContent(i_x,  i_y , v);
        target->SetBinError(i_x,    i_y , ev);

      }
      //T(i) = S(i-1) + 2S(i)/3   for last bin
      else if (i_x == n_x ) {
        v = (vGT1 + 2.0 * v ) / 3.0;
        ev = TMath::Sqrt( ( (1.0 / 3.0) * (1.0 / 3.0) * evGT1 * evGT1  + 2.0 / 3.0) * (2.0 / 3.0) * ev * ev);

        target->SetBinContent(i_x,  i_y , v);
        target->SetBinError(i_x,    i_y , ev);

      }
      else
      {
        v = (vGT1 + 2.0 * v + vLTnx ) / 4.0;
        ev = TMath::Sqrt( ( (1.0 / 4.0) * (1.0 / 4.0) * evGT1 * evGT1  + 2.0 / 4.0) * (2.0 / 4.0) * ev * ev + (1.0 / 4.0) * (1.0 / 4.0) * evLTnx * evLTnx);

        target->SetBinContent(i_x,  i_y , v);
        target->SetBinError(i_x,    i_y , ev);
      }

     
//cout << "i_y :" << i_y << "i_y + nbins :" << i_y + nbins << "\t";
    }//i_y loop over

  }//i_x loop over
int rebinX = 1; //for CD
      int rebinY = 2; //for CD
      target->Rebin2D(rebinX, rebinY);
      target->Scale(1. / float(rebinX) / float(rebinY));

}

//-------------calculate B2 -------------------
//B2 = <Nch>*R2cd
void calculate_B2( const TH2D *hp, const TH2D *hm, const TH2D *hin,  TH2D *hout)
{
  const Int_t knXBinNp = hp->GetNbinsX();
  const Int_t knYBinNp = hp->GetNbinsY();

  //For N+ , N- & Nch
  Double_t Np = 0.0, Nn = 0.0, Nch = 0.0;

  for (Int_t ieta = 1; ieta <= knXBinNp; ++ieta)
  {
    for (Int_t iphi = 1; iphi <= knYBinNp; ++iphi)
    {
    	if (ieta % 2 == 0 || iphi % 2 == 0 )
    	{
    gStyle->SetHistLineWidth(0);
    	}
      Np += hp->GetBinContent(ieta, iphi);
      Nn += hm->GetBinContent(ieta, iphi);
      Nch = Np + Nn;
    }
  }
  Nch = Nch / (Double_t(knXBinNp * knYBinNp));
  //  cout << "<Nch>:" << Nch << endl;
  hout->Add(hin, 1.0);
  hout->Scale(Nch);
}
