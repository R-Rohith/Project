#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<stdbool.h>
#include<iostream>
using namespace std;
void Cmbn(const TString start,const TString end)
{
	TH1I h1i_n1_multPM;
	TH1D h1d_n1_ptM,h1d_n1_ptP;
	TH2D 	h2d_n1_etaPhiP,h2d_n1_etaPhiM,h2d_pt_etaPhiP,h2d_pt_etaPhiM,
		h2d_n2_eta1Phi1Eta2Phi2PP,h2d_n2_eta1Phi1Eta2Phi2PM,h2d_n2_eta1Phi1Eta2Phi2MM,
		h2d_ptpt_eta1Phi1Eta2Phi2PM,h2d_ptpt_eta1Phi1Eta2Phi2MM,h2d_ptpt_eta1Phi1Eta2Phi2PP,
		h2d_ptn_eta1Phi1Eta2Phi2PM,h2d_ptn_eta1Phi1Eta2Phi2MM,h2d_ptn_eta1Phi1Eta2Phi2PP,
		h2d_npt_eta1Phi1Eta2Phi2PM,h2d_npt_eta1Phi1Eta2Phi2MM,h2d_npt_eta1Phi1Eta2Phi2PP;
	
	TH1I *mult; TH1D *ptM,*ptP;
	TH2D *n1p,*n1m,*ptp,*ptm,*n2pp,*n2pm,*n2mm,*nptpp,*nptpm,*nptmm,*ptnpp,*ptnpm,*ptnmm,*ptptpp,*ptptpm,*ptptmm;
	int i=atoi(start),j=atoi(end);
	char vers[10],file[50];
//	vers[0]='0'+(i/100);
	vers[0]='1';
	vers[1]='\0';
	for(int ver=i;ver<=j;ver++)
	{
//		vers=citoa(ver,vers,10);
		strcpy(file,"AnalysisResults");
		strcat(strcat(file,vers),".root");
		cout<<file;
//		if(ver==2)break;
		TFile *ifile=new TFile(file,"read");
		mult=(TH1I*)ifile->Get("r2/h1i_n1_multPM");
		ptM=(TH1D*)ifile->Get("r2/h1d_n1_ptM");
		ptP=(TH1D*)ifile->Get("r2/h1d_n1_ptP");
		n1p=(TH2D*)ifile->Get("r2/h2d_n1_etaPhiP");
		n1m=(TH2D*)ifile->Get("r2/h2d_n1_etaPhiM");
		ptp=(TH2D*)ifile->Get("r2/h2d_pt_etaPhiP");
		ptm=(TH2D*)ifile->Get("r2/h2d_pt_etaPhiM");
		n2pp=(TH2D*)ifile->Get("r2/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm=(TH2D*)ifile->Get("r2/h2d_n2_eta1Phi1Eta2Phi2PM");
		n2mm=(TH2D*)ifile->Get("r2/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("r2/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm=(TH2D*)ifile->Get("r2/h2d_ptpt_eta1Phi1Eta2Phi2PM");
		nptmm=(TH2D*)ifile->Get("r2/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("r2/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm=(TH2D*)ifile->Get("r2/h2d_ptn_eta1Phi1Eta2Phi2PM");
		ptnmm=(TH2D*)ifile->Get("r2/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("r2/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm=(TH2D*)ifile->Get("r2/h2d_npt_eta1Phi1Eta2Phi2PM");
		ptptmm=(TH2D*)ifile->Get("r2/h2d_npt_eta1Phi1Eta2Phi2MM");
		if(ver==i)
		{
			h1i_n1_multPM=*mult;
			h1d_n1_ptM=*ptM;
			h1d_n1_ptP=*ptP;
			h2d_n1_etaPhiP=*n1p;
			h2d_n1_etaPhiM=*n1m;
			h2d_pt_etaPhiP=*ptp;
			h2d_pt_etaPhiM=*ptm;
			h2d_n2_eta1Phi1Eta2Phi2PP=*n2pp;
			h2d_n2_eta1Phi1Eta2Phi2PM=*n2pm;
			h2d_n2_eta1Phi1Eta2Phi2MM=*n2mm;
			h2d_ptpt_eta1Phi1Eta2Phi2PM=*ptptpm;
			h2d_ptpt_eta1Phi1Eta2Phi2MM=*ptptmm;
			h2d_ptpt_eta1Phi1Eta2Phi2PP=*ptptpp;
			h2d_ptn_eta1Phi1Eta2Phi2PM=*ptnpm;
			h2d_ptn_eta1Phi1Eta2Phi2MM=*ptnmm;
			h2d_ptn_eta1Phi1Eta2Phi2PP=*ptnpp;
			h2d_npt_eta1Phi1Eta2Phi2PM=*nptpm;
			h2d_npt_eta1Phi1Eta2Phi2MM=*nptmm;
			h2d_npt_eta1Phi1Eta2Phi2PP=*nptpp;
			continue;
		}
		h1i_n1_multPM=h1i_n1_multPM+*mult;
		h1d_n1_ptM=h1d_n1_ptM+*ptM;
		h1d_n1_ptP=h1d_n1_ptP+*ptP;
		h2d_n1_etaPhiP=h2d_n1_etaPhiP+*n1p;
		h2d_n1_etaPhiM=h2d_n1_etaPhiM+*n1m;
		h2d_pt_etaPhiP=h2d_pt_etaPhiP+*ptp;
		h2d_pt_etaPhiM=h2d_pt_etaPhiM+*ptm;
		h2d_n2_eta1Phi1Eta2Phi2PP=h2d_n2_eta1Phi1Eta2Phi2PP+*n2pp;
		h2d_n2_eta1Phi1Eta2Phi2PM=h2d_n2_eta1Phi1Eta2Phi2PM+*n2pm;
		h2d_n2_eta1Phi1Eta2Phi2MM=h2d_n2_eta1Phi1Eta2Phi2MM+*n2mm;
		h2d_ptpt_eta1Phi1Eta2Phi2PM=h2d_ptpt_eta1Phi1Eta2Phi2PM+*ptptpm;
		h2d_ptpt_eta1Phi1Eta2Phi2MM=h2d_ptpt_eta1Phi1Eta2Phi2MM+*ptptmm;
		h2d_ptpt_eta1Phi1Eta2Phi2PP=h2d_ptpt_eta1Phi1Eta2Phi2PP+*ptptpp;
		h2d_ptn_eta1Phi1Eta2Phi2PM=h2d_ptn_eta1Phi1Eta2Phi2PM+*ptnpm;
		h2d_ptn_eta1Phi1Eta2Phi2MM=h2d_ptn_eta1Phi1Eta2Phi2MM+*ptnmm;
		h2d_ptn_eta1Phi1Eta2Phi2PP=h2d_ptn_eta1Phi1Eta2Phi2PP+*ptnpp;
		h2d_npt_eta1Phi1Eta2Phi2PM=h2d_npt_eta1Phi1Eta2Phi2PM+*nptpm;
		h2d_npt_eta1Phi1Eta2Phi2MM=h2d_npt_eta1Phi1Eta2Phi2MM+*nptmm;
		h2d_npt_eta1Phi1Eta2Phi2PP=h2d_npt_eta1Phi1Eta2Phi2PP+*nptpp;
		vers[0]+=1;
	}
	TFile *ofile=new TFile("tst.root","recreate");
	ofile->cd();
		h1i_n1_multPM.Write();
		h1d_n1_ptM.Write();
		h1d_n1_ptP.Write();
		h2d_n1_etaPhiP.Write();
		h2d_n1_etaPhiM.Write();
		h2d_pt_etaPhiP.Write();
		h2d_pt_etaPhiM.Write();
		h2d_n2_eta1Phi1Eta2Phi2PP.Write();
		h2d_n2_eta1Phi1Eta2Phi2PM.Write();
		h2d_n2_eta1Phi1Eta2Phi2MM.Write();
		h2d_ptpt_eta1Phi1Eta2Phi2PM.Write();
		h2d_ptpt_eta1Phi1Eta2Phi2MM.Write();
		h2d_ptpt_eta1Phi1Eta2Phi2PP.Write();
		h2d_ptn_eta1Phi1Eta2Phi2PM.Write();
		h2d_ptn_eta1Phi1Eta2Phi2MM.Write();
		h2d_ptn_eta1Phi1Eta2Phi2PP.Write();
		h2d_npt_eta1Phi1Eta2Phi2PM.Write();
		h2d_npt_eta1Phi1Eta2Phi2MM.Write();
		h2d_npt_eta1Phi1Eta2Phi2PP.Write();
}

