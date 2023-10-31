#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<stdbool.h>
#include<iostream>
using namespace std;
void Cmbn(int i,int j)//(const char* start,const TString end)
{
	void itoa(int,char*);

	TH1I h1i_n1_multPM;
	TH1D h1d_n1_ptM,h1d_n1_ptP;
	TH2D 	h2d_n1_etaPhiP,h2d_n1_etaPhiM,h2d_pt_etaPhiP,h2d_pt_etaPhiM,
		h2d_n2_eta1Phi1Eta2Phi2PP,h2d_n2_eta1Phi1Eta2Phi2PM,h2d_n2_eta1Phi1Eta2Phi2MM,
		h2d_ptpt_eta1Phi1Eta2Phi2PM,h2d_ptpt_eta1Phi1Eta2Phi2MM,h2d_ptpt_eta1Phi1Eta2Phi2PP,
		h2d_ptn_eta1Phi1Eta2Phi2PM,h2d_ptn_eta1Phi1Eta2Phi2MM,h2d_ptn_eta1Phi1Eta2Phi2PP,
		h2d_npt_eta1Phi1Eta2Phi2PM,h2d_npt_eta1Phi1Eta2Phi2MM,h2d_npt_eta1Phi1Eta2Phi2PP;
	
	TH1I *mult; TH1D *ptM,*ptP;
	TH2D *n1p,*n1m,*ptp,*ptm,*n2pp,*n2pm,*n2mm,*nptpp,*nptpm,*nptmm,*ptnpp,*ptnpm,*ptnmm,*ptptpp,*ptptpm,*ptptmm;
//	int i=atoi(start),j=atoi(end);
	char vers[5],file[50];
//	vers[0]=48+i/10;
//	vers[1]=48+i%10;
//	vers[2]='\0';
	for(int ver=i;ver<=j;ver++)
	{
		vers[0]='\0';
		itoa(ver,vers);
//		vers=citoa(ver,vers,10);
		strcpy(file,"");
		strcat(strcat(file,vers),"/AnalysisResults.root");
//		cout<<vers;
		cout<<"\nProcessing: "<<file<<endl;
//		if(ver==2)break;
		TFile *ifile=new TFile(file,"read");
		if (!ifile || ifile->IsZombie()) {
		   std::cerr << "Error opening file" << endl;
		   exit(-1);
		}
		mult=(TH1I*)ifile->Get("r2p24id/h1i_n1_multPM");
		ptM=(TH1D*)ifile->Get("r2p24id/h1d_n1_ptM");
		ptP=(TH1D*)ifile->Get("r2p24id/h1d_n1_ptP");
		n1p=(TH2D*)ifile->Get("r2p24id/h2d_n1_etaPhiP");
		n1m=(TH2D*)ifile->Get("r2p24id/h2d_n1_etaPhiM");
		ptp=(TH2D*)ifile->Get("r2p24id/h2d_pt_etaPhiP");
		ptm=(TH2D*)ifile->Get("r2p24id/h2d_pt_etaPhiM");
		n2pp=(TH2D*)ifile->Get("r2p24id/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm=(TH2D*)ifile->Get("r2p24id/h2d_n2_eta1Phi1Eta2Phi2PM");
		n2mm=(TH2D*)ifile->Get("r2p24id/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("r2p24id/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm=(TH2D*)ifile->Get("r2p24id/h2d_ptpt_eta1Phi1Eta2Phi2PM");
		nptmm=(TH2D*)ifile->Get("r2p24id/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("r2p24id/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm=(TH2D*)ifile->Get("r2p24id/h2d_ptn_eta1Phi1Eta2Phi2PM");
		ptnmm=(TH2D*)ifile->Get("r2p24id/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("r2p24id/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm=(TH2D*)ifile->Get("r2p24id/h2d_npt_eta1Phi1Eta2Phi2PM");
		ptptmm=(TH2D*)ifile->Get("r2p24id/h2d_npt_eta1Phi1Eta2Phi2MM");
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
		}
		else
		{
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
		}
//		if(vers[0]!=48+ver/10)
//		{
//			vers[0]+=1;
//			vers[1]='0';
//		}
//		else
//			vers[1]+=1;
		delete ifile;

	}
	TFile *ofile=new TFile("AnalysisResults.root","recreate");
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
		delete ofile;
}
void itoa(int num,char* str)
{
	char ch[5]="5";
	if(num==0) return;
	ch[0]=48+num%10;
//	cout<<ch<<endl;
	strcat(ch,str);
	strcpy(str,ch);
	itoa(num/10,str);
}

