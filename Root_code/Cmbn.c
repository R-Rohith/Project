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
	TH1D h1d_n1_ptM,h1d_n1_ptP,h1d_n1_phi,h1d_n1_eta,h1d_n1_pt;
	TH2D 	h2d_n1_etaPhiP,h2d_n1_etaPhiM,h2d_pt_etaPhiP,h2d_pt_etaPhiM,
		h2d_n2_eta1Phi1Eta2Phi2PP,h2d_n2_eta1Phi1Eta2Phi2PM,h2d_n2_eta1Phi1Eta2Phi2MM,
		h2d_ptpt_eta1Phi1Eta2Phi2PM,h2d_ptpt_eta1Phi1Eta2Phi2MM,h2d_ptpt_eta1Phi1Eta2Phi2PP,
		h2d_ptn_eta1Phi1Eta2Phi2PM,h2d_ptn_eta1Phi1Eta2Phi2MM,h2d_ptn_eta1Phi1Eta2Phi2PP,
		h2d_npt_eta1Phi1Eta2Phi2PM,h2d_npt_eta1Phi1Eta2Phi2MM,h2d_npt_eta1Phi1Eta2Phi2PP;
	
	TH1I *mult; TH1D *ptM,*ptP,*phi,*eta,*pt;
	TH2D *n1p,*n1m,*ptp,*ptm,*n2pp,*n2pm,*n2mm,*nptpp,*nptpm,*nptmm,*ptnpp,*ptnpm,*ptnmm,*ptptpp,*ptptpm,*ptptmm;
//	int i=atoi(start),j=atoi(end);
	char vers[5],file[50];
//	vers[0]=48+i/10;
//	vers[1]=48+i%10;
//	vers[2]='\0';
	for(int num =0;num<4;num++){
	for(int ver=i;ver<=j;ver++)
	{
		vers[0]='\0';
		itoa(ver,vers);
//		vers=citoa(ver,vers,10);
		strcpy(file,"../");
		strcat(strcat(file,vers),"/AnalysisResults.root");
//		cout<<vers;
		cout<<"\nprocessing: "<<file<<endl;
//		if(ver==2)break;
		TFile *ifile=new TFile(file,"read");
		if (!ifile || ifile->IsZombie()) {
		   std::cerr << "Error opening file" << endl;
		   exit(-1);
		}
		if(num==0){cout<<"charged\n";
		phi=(TH1D*)ifile->Get("r2p24ch/phi");
		eta=(TH1D*)ifile->Get("r2p24ch/eta");
		pt=(TH1D*)ifile->Get("r2p24ch/pt1");
		mult=(TH1I*)ifile->Get("r2p24ch/h1i_n1_multPM");
		ptM=(TH1D*)ifile->Get("r2p24ch/h1d_n1_ptM");
		ptP=(TH1D*)ifile->Get("r2p24ch/h1d_n1_ptP");
		n1p=(TH2D*)ifile->Get("r2p24ch/h2d_n1_etaPhiP");
		n1m=(TH2D*)ifile->Get("r2p24ch/h2d_n1_etaPhiM");
		ptp=(TH2D*)ifile->Get("r2p24ch/h2d_pt_etaPhiP");
		ptm=(TH2D*)ifile->Get("r2p24ch/h2d_pt_etaPhiM");
		n2pp=(TH2D*)ifile->Get("r2p24ch/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm=(TH2D*)ifile->Get("r2p24ch/h2d_n2_eta1Phi1Eta2Phi2PM");
		n2mm=(TH2D*)ifile->Get("r2p24ch/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("r2p24ch/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm=(TH2D*)ifile->Get("r2p24ch/h2d_ptpt_eta1Phi1Eta2Phi2PM");
		nptmm=(TH2D*)ifile->Get("r2p24ch/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("r2p24ch/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm=(TH2D*)ifile->Get("r2p24ch/h2d_ptn_eta1Phi1Eta2Phi2PM");
		ptnmm=(TH2D*)ifile->Get("r2p24ch/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("r2p24ch/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm=(TH2D*)ifile->Get("r2p24ch/h2d_npt_eta1Phi1Eta2Phi2PM");
		ptptmm=(TH2D*)ifile->Get("r2p24ch/h2d_npt_eta1Phi1Eta2Phi2MM");
		}
		if(num==1){cout<<"pion\n";
		phi=(TH1D*)ifile->Get("r2p24pi/phi");
		eta=(TH1D*)ifile->Get("r2p24pi/eta");
		pt=(TH1D*)ifile->Get("r2p24pi/pt1");
		mult=(TH1I*)ifile->Get("r2p24pi/h1i_n1_multPM");
		ptM=(TH1D*)ifile->Get("r2p24pi/h1d_n1_ptM");
		ptP=(TH1D*)ifile->Get("r2p24pi/h1d_n1_ptP");
		n1p=(TH2D*)ifile->Get("r2p24pi/h2d_n1_etaPhiP");
		n1m=(TH2D*)ifile->Get("r2p24pi/h2d_n1_etaPhiM");
		ptp=(TH2D*)ifile->Get("r2p24pi/h2d_pt_etaPhiP");
		ptm=(TH2D*)ifile->Get("r2p24pi/h2d_pt_etaPhiM");
		n2pp=(TH2D*)ifile->Get("r2p24pi/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm=(TH2D*)ifile->Get("r2p24pi/h2d_n2_eta1Phi1Eta2Phi2PM");
		n2mm=(TH2D*)ifile->Get("r2p24pi/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("r2p24pi/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm=(TH2D*)ifile->Get("r2p24pi/h2d_ptpt_eta1Phi1Eta2Phi2PM");
		nptmm=(TH2D*)ifile->Get("r2p24pi/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("r2p24pi/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm=(TH2D*)ifile->Get("r2p24pi/h2d_ptn_eta1Phi1Eta2Phi2PM");
		ptnmm=(TH2D*)ifile->Get("r2p24pi/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("r2p24pi/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm=(TH2D*)ifile->Get("r2p24pi/h2d_npt_eta1Phi1Eta2Phi2PM");
		ptptmm=(TH2D*)ifile->Get("r2p24pi/h2d_npt_eta1Phi1Eta2Phi2MM");
		}
		if(num==2){cout<<"kaon\n";
		phi=(TH1D*)ifile->Get("r2p24ka/phi");
		eta=(TH1D*)ifile->Get("r2p24ka/eta");
		pt=(TH1D*)ifile->Get("r2p24ka/pt1");
		mult=(TH1I*)ifile->Get("r2p24ka/h1i_n1_multPM");
		ptM=(TH1D*)ifile->Get("r2p24ka/h1d_n1_ptM");
		ptP=(TH1D*)ifile->Get("r2p24ka/h1d_n1_ptP");
		n1p=(TH2D*)ifile->Get("r2p24ka/h2d_n1_etaPhiP");
		n1m=(TH2D*)ifile->Get("r2p24ka/h2d_n1_etaPhiM");
		ptp=(TH2D*)ifile->Get("r2p24ka/h2d_pt_etaPhiP");
		ptm=(TH2D*)ifile->Get("r2p24ka/h2d_pt_etaPhiM");
		n2pp=(TH2D*)ifile->Get("r2p24ka/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm=(TH2D*)ifile->Get("r2p24ka/h2d_n2_eta1Phi1Eta2Phi2PM");
		n2mm=(TH2D*)ifile->Get("r2p24ka/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("r2p24ka/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm=(TH2D*)ifile->Get("r2p24ka/h2d_ptpt_eta1Phi1Eta2Phi2PM");
		nptmm=(TH2D*)ifile->Get("r2p24ka/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("r2p24ka/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm=(TH2D*)ifile->Get("r2p24ka/h2d_ptn_eta1Phi1Eta2Phi2PM");
		ptnmm=(TH2D*)ifile->Get("r2p24ka/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("r2p24ka/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm=(TH2D*)ifile->Get("r2p24ka/h2d_npt_eta1Phi1Eta2Phi2PM");
		ptptmm=(TH2D*)ifile->Get("r2p24ka/h2d_npt_eta1Phi1Eta2Phi2MM");
		}
		if(num==3){cout<<"proton\n";
		phi=(TH1D*)ifile->Get("r2p24pr/phi");
		eta=(TH1D*)ifile->Get("r2p24pr/eta");
		pt=(TH1D*)ifile->Get("r2p24pr/pt1");
		mult=(TH1I*)ifile->Get("r2p24pr/h1i_n1_multPM");
		ptM=(TH1D*)ifile->Get("r2p24pr/h1d_n1_ptM");
		ptP=(TH1D*)ifile->Get("r2p24pr/h1d_n1_ptP");
		n1p=(TH2D*)ifile->Get("r2p24pr/h2d_n1_etaPhiP");
		n1m=(TH2D*)ifile->Get("r2p24pr/h2d_n1_etaPhiM");
		ptp=(TH2D*)ifile->Get("r2p24pr/h2d_pt_etaPhiP");
		ptm=(TH2D*)ifile->Get("r2p24pr/h2d_pt_etaPhiM");
		n2pp=(TH2D*)ifile->Get("r2p24pr/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm=(TH2D*)ifile->Get("r2p24pr/h2d_n2_eta1Phi1Eta2Phi2PM");
		n2mm=(TH2D*)ifile->Get("r2p24pr/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("r2p24pr/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm=(TH2D*)ifile->Get("r2p24pr/h2d_ptpt_eta1Phi1Eta2Phi2PM");
		nptmm=(TH2D*)ifile->Get("r2p24pr/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("r2p24pr/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm=(TH2D*)ifile->Get("r2p24pr/h2d_ptn_eta1Phi1Eta2Phi2PM");
		ptnmm=(TH2D*)ifile->Get("r2p24pr/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("r2p24pr/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm=(TH2D*)ifile->Get("r2p24pr/h2d_npt_eta1Phi1Eta2Phi2PM");
		ptptmm=(TH2D*)ifile->Get("r2p24pr/h2d_npt_eta1Phi1Eta2Phi2MM");
		}
		if(ver==i)
		{
			h1d_n1_phi=*phi;
			h1d_n1_eta=*eta;
			h1d_n1_pt=*pt;
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
			h1d_n1_phi=h1d_n1_phi+*phi;
			h1d_n1_eta=h1d_n1_eta+*eta;
			h1d_n1_pt=h1d_n1_pt+*pt;
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
	TFile *ofile;
	if(num==0)
	ofile=new TFile("AnalysisResults_Ch.root","recreate");
	if(num==1)
	ofile=new TFile("AnalysisResults_Pi.root","recreate");
	if(num==2)
	ofile=new TFile("AnalysisResults_Ka.root","recreate");
	if(num==3)
	ofile=new TFile("AnalysisResults_Pr.root","recreate");
	ofile->cd();
		h1d_n1_phi.Write();
		h1d_n1_eta.Write();
		h1d_n1_pt.Write();
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

/*

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
	TH1D h1d_n1_ptM1,h1d_n1_ptP1,h1d_n1_ptM2,h1d_n1_ptP2,h1d_n1_phi,h1d_n1_eta,h1d_n1_pt;
	TH2D 	h2d_n1_etaPhiP1,h2d_n1_etaPhiM1,h2d_pt_etaPhiP1,h2d_pt_etaPhiM1,h2d_n1_etaPhiP2,h2d_n1_etaPhiM2,h2d_pt_etaPhiP2,h2d_pt_etaPhiM2,
		h2d_n2_eta1Phi1Eta2Phi2PP,h2d_n2_eta1Phi1Eta2Phi2PM12,h2d_n2_eta1Phi1Eta2Phi2PM21,h2d_n2_eta1Phi1Eta2Phi2MM,
		h2d_ptpt_eta1Phi1Eta2Phi2PM12,h2d_ptpt_eta1Phi1Eta2Phi2PM21,h2d_ptpt_eta1Phi1Eta2Phi2MM,h2d_ptpt_eta1Phi1Eta2Phi2PP,
		h2d_ptn_eta1Phi1Eta2Phi2PM12,h2d_ptn_eta1Phi1Eta2Phi2PM21,h2d_ptn_eta1Phi1Eta2Phi2MM,h2d_ptn_eta1Phi1Eta2Phi2PP,
		h2d_npt_eta1Phi1Eta2Phi2PM12,h2d_npt_eta1Phi1Eta2Phi2PM21,h2d_npt_eta1Phi1Eta2Phi2MM,h2d_npt_eta1Phi1Eta2Phi2PP;
	
	TH1I *mult; TH1D *ptM1,*ptP1,*ptM2,*ptP2,*phi,*eta,*pt;
	TH2D *n1p1,*n1m1,*ptp1,*ptm1,*n1p2,*n1m2,*ptp2,*ptm2,*n2pp,*n2pm12,*n2pm21,*n2mm,*nptpp,*nptpm12,*nptpm21,*nptmm,*ptnpp,*ptnpm12,*ptnpm21,*ptnmm,*ptptpp,*ptptpm12,*ptptpm21,*ptptmm;
//	int i=atoi(start),j=atoi(end);
	char vers[5],file[50];
//	vers[0]=48+i/10;
//	vers[1]=48+i%10;
//	vers[2]='\0';
	for(int num=0;num<3;num++){
	for(int ver=i;ver<=j;ver++)
	{
		vers[0]='\0';
		itoa(ver,vers);
//		vers=citoa(ver,vers,10);
		strcpy(file,"../");
		strcat(strcat(file,vers),"/AnalysisResults.root");
//		cout<<vers;
		cout<<"\nProcessing: "<<file<<endl;
//		if(ver==2)break;
		TFile *ifile=new TFile(file,"read");
		if (!ifile || ifile->IsZombie()) {
		   std::cerr << "Error opening file" << endl;
		   exit(-1);
		}
		if(num==0){cout<<"Pion-Kaon\n";
		phi=(TH1D*)ifile->Get("crossr2p24pik/phi");
		eta=(TH1D*)ifile->Get("crossr2p24pik/eta");
		pt=(TH1D*)ifile->Get("crossr2p24pik/pt1");
		mult=(TH1I*)ifile->Get("crossr2p24pik/h1i_n1_multPM");
		ptM1=(TH1D*)ifile->Get("crossr2p24pik/h1d_n1_ptM1");
		ptP1=(TH1D*)ifile->Get("crossr2p24pik/h1d_n1_ptP1");
		ptM2=(TH1D*)ifile->Get("crossr2p24pik/h1d_n1_ptM2");
		ptP2=(TH1D*)ifile->Get("crossr2p24pik/h1d_n1_ptP2");
		n1p1=(TH2D*)ifile->Get("crossr2p24pik/h2d_n1_etaPhiP1");
		n1m1=(TH2D*)ifile->Get("crossr2p24pik/h2d_n1_etaPhiM1");
		ptp1=(TH2D*)ifile->Get("crossr2p24pik/h2d_pt_etaPhiP1");
		ptm1=(TH2D*)ifile->Get("crossr2p24pik/h2d_pt_etaPhiM1");
		n1p2=(TH2D*)ifile->Get("crossr2p24pik/h2d_n1_etaPhiP2");
		n1m2=(TH2D*)ifile->Get("crossr2p24pik/h2d_n1_etaPhiM2");
		ptp2=(TH2D*)ifile->Get("crossr2p24pik/h2d_pt_etaPhiP2");
		ptm2=(TH2D*)ifile->Get("crossr2p24pik/h2d_pt_etaPhiM2");
		n2pp=(TH2D*)ifile->Get("crossr2p24pik/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm12=(TH2D*)ifile->Get("crossr2p24pik/h2d_n2_eta1Phi1Eta2Phi2PM12");
		n2pm21=(TH2D*)ifile->Get("crossr2p24pik/h2d_n2_eta1Phi1Eta2Phi2PM21");
		n2mm=(TH2D*)ifile->Get("crossr2p24pik/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("crossr2p24pik/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm12=(TH2D*)ifile->Get("crossr2p24pik/h2d_ptpt_eta1Phi1Eta2Phi2PM12");
		nptpm21=(TH2D*)ifile->Get("crossr2p24pik/h2d_ptpt_eta1Phi1Eta2Phi2PM21");
		nptmm=(TH2D*)ifile->Get("crossr2p24pik/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("crossr2p24pik/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm12=(TH2D*)ifile->Get("crossr2p24pik/h2d_ptn_eta1Phi1Eta2Phi2PM12");
		ptnpm21=(TH2D*)ifile->Get("crossr2p24pik/h2d_ptn_eta1Phi1Eta2Phi2PM21");
		ptnmm=(TH2D*)ifile->Get("crossr2p24pik/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("crossr2p24pik/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm12=(TH2D*)ifile->Get("crossr2p24pik/h2d_npt_eta1Phi1Eta2Phi2PM12");
		ptptpm21=(TH2D*)ifile->Get("crossr2p24pik/h2d_npt_eta1Phi1Eta2Phi2PM21");
		ptptmm=(TH2D*)ifile->Get("crossr2p24pik/h2d_npt_eta1Phi1Eta2Phi2MM");
		}
		if(num==1){cout<<"Pion-Proton\n";
		phi=(TH1D*)ifile->Get("crossr2p24pip/phi");
		eta=(TH1D*)ifile->Get("crossr2p24pip/eta");
		pt=(TH1D*)ifile->Get("crossr2p24pip/pt1");
		mult=(TH1I*)ifile->Get("crossr2p24pip/h1i_n1_multPM");
		ptM1=(TH1D*)ifile->Get("crossr2p24pip/h1d_n1_ptM1");
		ptP1=(TH1D*)ifile->Get("crossr2p24pip/h1d_n1_ptP1");
		ptM2=(TH1D*)ifile->Get("crossr2p24pip/h1d_n1_ptM2");
		ptP2=(TH1D*)ifile->Get("crossr2p24pip/h1d_n1_ptP2");
		n1p1=(TH2D*)ifile->Get("crossr2p24pip/h2d_n1_etaPhiP1");
		n1m1=(TH2D*)ifile->Get("crossr2p24pip/h2d_n1_etaPhiM1");
		ptp1=(TH2D*)ifile->Get("crossr2p24pip/h2d_pt_etaPhiP1");
		ptm1=(TH2D*)ifile->Get("crossr2p24pip/h2d_pt_etaPhiM1");
		n1p2=(TH2D*)ifile->Get("crossr2p24pip/h2d_n1_etaPhiP2");
		n1m2=(TH2D*)ifile->Get("crossr2p24pip/h2d_n1_etaPhiM2");
		ptp2=(TH2D*)ifile->Get("crossr2p24pip/h2d_pt_etaPhiP2");
		ptm2=(TH2D*)ifile->Get("crossr2p24pip/h2d_pt_etaPhiM2");
		n2pp=(TH2D*)ifile->Get("crossr2p24pip/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm12=(TH2D*)ifile->Get("crossr2p24pip/h2d_n2_eta1Phi1Eta2Phi2PM12");
		n2pm21=(TH2D*)ifile->Get("crossr2p24pip/h2d_n2_eta1Phi1Eta2Phi2PM21");
		n2mm=(TH2D*)ifile->Get("crossr2p24pip/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("crossr2p24pip/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm12=(TH2D*)ifile->Get("crossr2p24pip/h2d_ptpt_eta1Phi1Eta2Phi2PM12");
		nptpm21=(TH2D*)ifile->Get("crossr2p24pip/h2d_ptpt_eta1Phi1Eta2Phi2PM21");
		nptmm=(TH2D*)ifile->Get("crossr2p24pip/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("crossr2p24pip/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm12=(TH2D*)ifile->Get("crossr2p24pip/h2d_ptn_eta1Phi1Eta2Phi2PM12");
		ptnpm21=(TH2D*)ifile->Get("crossr2p24pip/h2d_ptn_eta1Phi1Eta2Phi2PM21");
		ptnmm=(TH2D*)ifile->Get("crossr2p24pip/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("crossr2p24pip/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm12=(TH2D*)ifile->Get("crossr2p24pip/h2d_npt_eta1Phi1Eta2Phi2PM12");
		ptptpm21=(TH2D*)ifile->Get("crossr2p24pip/h2d_npt_eta1Phi1Eta2Phi2PM21");
		ptptmm=(TH2D*)ifile->Get("crossr2p24pip/h2d_npt_eta1Phi1Eta2Phi2MM");
		}
		if(num==2){cout<<"Proton-Kaon\n";
		phi=(TH1D*)ifile->Get("crossr2p24pk/phi");
		eta=(TH1D*)ifile->Get("crossr2p24pk/eta");
		pt=(TH1D*)ifile->Get("crossr2p24pk/pt1");
		mult=(TH1I*)ifile->Get("crossr2p24pk/h1i_n1_multPM");
		ptM1=(TH1D*)ifile->Get("crossr2p24pk/h1d_n1_ptM1");
		ptP1=(TH1D*)ifile->Get("crossr2p24pk/h1d_n1_ptP1");
		ptM2=(TH1D*)ifile->Get("crossr2p24pk/h1d_n1_ptM2");
		ptP2=(TH1D*)ifile->Get("crossr2p24pk/h1d_n1_ptP2");
		n1p1=(TH2D*)ifile->Get("crossr2p24pk/h2d_n1_etaPhiP1");
		n1m1=(TH2D*)ifile->Get("crossr2p24pk/h2d_n1_etaPhiM1");
		ptp1=(TH2D*)ifile->Get("crossr2p24pk/h2d_pt_etaPhiP1");
		ptm1=(TH2D*)ifile->Get("crossr2p24pk/h2d_pt_etaPhiM1");
		n1p2=(TH2D*)ifile->Get("crossr2p24pk/h2d_n1_etaPhiP2");
		n1m2=(TH2D*)ifile->Get("crossr2p24pk/h2d_n1_etaPhiM2");
		ptp2=(TH2D*)ifile->Get("crossr2p24pk/h2d_pt_etaPhiP2");
		ptm2=(TH2D*)ifile->Get("crossr2p24pk/h2d_pt_etaPhiM2");
		n2pp=(TH2D*)ifile->Get("crossr2p24pk/h2d_n2_eta1Phi1Eta2Phi2PP");
		n2pm12=(TH2D*)ifile->Get("crossr2p24pk/h2d_n2_eta1Phi1Eta2Phi2PM12");
		n2pm21=(TH2D*)ifile->Get("crossr2p24pk/h2d_n2_eta1Phi1Eta2Phi2PM21");
		n2mm=(TH2D*)ifile->Get("crossr2p24pk/h2d_n2_eta1Phi1Eta2Phi2MM");
		nptpp=(TH2D*)ifile->Get("crossr2p24pk/h2d_ptpt_eta1Phi1Eta2Phi2PP");
		nptpm12=(TH2D*)ifile->Get("crossr2p24pk/h2d_ptpt_eta1Phi1Eta2Phi2PM12");
		nptpm21=(TH2D*)ifile->Get("crossr2p24pk/h2d_ptpt_eta1Phi1Eta2Phi2PM21");
		nptmm=(TH2D*)ifile->Get("crossr2p24pk/h2d_ptpt_eta1Phi1Eta2Phi2MM");
		ptnpp=(TH2D*)ifile->Get("crossr2p24pk/h2d_ptn_eta1Phi1Eta2Phi2PP");
		ptnpm12=(TH2D*)ifile->Get("crossr2p24pk/h2d_ptn_eta1Phi1Eta2Phi2PM12");
		ptnpm21=(TH2D*)ifile->Get("crossr2p24pk/h2d_ptn_eta1Phi1Eta2Phi2PM21");
		ptnmm=(TH2D*)ifile->Get("crossr2p24pk/h2d_ptn_eta1Phi1Eta2Phi2MM");
		ptptpp=(TH2D*)ifile->Get("crossr2p24pk/h2d_npt_eta1Phi1Eta2Phi2PP");
		ptptpm12=(TH2D*)ifile->Get("crossr2p24pk/h2d_npt_eta1Phi1Eta2Phi2PM12");
		ptptpm21=(TH2D*)ifile->Get("crossr2p24pk/h2d_npt_eta1Phi1Eta2Phi2PM21");
		ptptmm=(TH2D*)ifile->Get("crossr2p24pk/h2d_npt_eta1Phi1Eta2Phi2MM");
		}
		if(ver==i)
		{
			h1d_n1_phi=*phi;
			h1d_n1_eta=*eta;
			h1d_n1_pt=*pt;
			h1i_n1_multPM=*mult;
			h1d_n1_ptM1=*ptM1;
			h1d_n1_ptP1=*ptP1;
			h2d_n1_etaPhiP1=*n1p1;
			h2d_n1_etaPhiM1=*n1m1;
			h2d_pt_etaPhiP1=*ptp1;
			h2d_pt_etaPhiM1=*ptm1;
			h1d_n1_ptM2=*ptM2;
			h1d_n1_ptP2=*ptP2;
			h2d_n1_etaPhiP2=*n1p2;
			h2d_n1_etaPhiM2=*n1m2;
			h2d_pt_etaPhiP2=*ptp2;
			h2d_pt_etaPhiM2=*ptm2;
			h2d_n2_eta1Phi1Eta2Phi2PP=*n2pp;
			h2d_n2_eta1Phi1Eta2Phi2PM12=*n2pm12;
			h2d_n2_eta1Phi1Eta2Phi2PM21=*n2pm21;
			h2d_n2_eta1Phi1Eta2Phi2MM=*n2mm;
			h2d_ptpt_eta1Phi1Eta2Phi2PM12=*ptptpm12;
			h2d_ptpt_eta1Phi1Eta2Phi2PM21=*ptptpm21;
			h2d_ptpt_eta1Phi1Eta2Phi2MM=*ptptmm;
			h2d_ptpt_eta1Phi1Eta2Phi2PP=*ptptpp;
			h2d_ptn_eta1Phi1Eta2Phi2PM12=*ptnpm12;
			h2d_ptn_eta1Phi1Eta2Phi2PM21=*ptnpm21;
			h2d_ptn_eta1Phi1Eta2Phi2MM=*ptnmm;
			h2d_ptn_eta1Phi1Eta2Phi2PP=*ptnpp;
			h2d_npt_eta1Phi1Eta2Phi2PM12=*nptpm12;
			h2d_npt_eta1Phi1Eta2Phi2PM21=*nptpm21;
			h2d_npt_eta1Phi1Eta2Phi2MM=*nptmm;
			h2d_npt_eta1Phi1Eta2Phi2PP=*nptpp;
		}
		else
		{
			h1d_n1_phi=h1d_n1_phi+*phi;
			h1d_n1_eta=h1d_n1_eta+*eta;
			h1d_n1_pt=h1d_n1_pt+*pt;
			h1i_n1_multPM=h1i_n1_multPM+*mult;
			h1d_n1_ptM1=h1d_n1_ptM1+*ptM1;
			h1d_n1_ptP1=h1d_n1_ptP1+*ptP1;
			h2d_n1_etaPhiP1=h2d_n1_etaPhiP1+*n1p1;
			h2d_n1_etaPhiM1=h2d_n1_etaPhiM1+*n1m1;
			h2d_pt_etaPhiP1=h2d_pt_etaPhiP1+*ptp1;
			h2d_pt_etaPhiM1=h2d_pt_etaPhiM1+*ptm1;
			h1d_n1_ptM2=h1d_n1_ptM2+*ptM2;
			h1d_n1_ptP2=h1d_n1_ptP2+*ptP2;
			h2d_n1_etaPhiP2=h2d_n1_etaPhiP2+*n1p2;
			h2d_n1_etaPhiM2=h2d_n1_etaPhiM2+*n1m2;
			h2d_pt_etaPhiP2=h2d_pt_etaPhiP2+*ptp2;
			h2d_pt_etaPhiM2=h2d_pt_etaPhiM2+*ptm2;
			h2d_n2_eta1Phi1Eta2Phi2PP=h2d_n2_eta1Phi1Eta2Phi2PP+*n2pp;
			h2d_n2_eta1Phi1Eta2Phi2PM12=h2d_n2_eta1Phi1Eta2Phi2PM12+*n2pm12;
			h2d_n2_eta1Phi1Eta2Phi2PM21=h2d_n2_eta1Phi1Eta2Phi2PM21+*n2pm21;
			h2d_n2_eta1Phi1Eta2Phi2MM=h2d_n2_eta1Phi1Eta2Phi2MM+*n2mm;
			h2d_ptpt_eta1Phi1Eta2Phi2PM12=h2d_ptpt_eta1Phi1Eta2Phi2PM12+*ptptpm12;
			h2d_ptpt_eta1Phi1Eta2Phi2PM21=h2d_ptpt_eta1Phi1Eta2Phi2PM21+*ptptpm21;
			h2d_ptpt_eta1Phi1Eta2Phi2MM=h2d_ptpt_eta1Phi1Eta2Phi2MM+*ptptmm;
			h2d_ptpt_eta1Phi1Eta2Phi2PP=h2d_ptpt_eta1Phi1Eta2Phi2PP+*ptptpp;
			h2d_ptn_eta1Phi1Eta2Phi2PM12=h2d_ptn_eta1Phi1Eta2Phi2PM12+*ptnpm12;
			h2d_ptn_eta1Phi1Eta2Phi2PM21=h2d_ptn_eta1Phi1Eta2Phi2PM21+*ptnpm21;
			h2d_ptn_eta1Phi1Eta2Phi2MM=h2d_ptn_eta1Phi1Eta2Phi2MM+*ptnmm;
			h2d_ptn_eta1Phi1Eta2Phi2PP=h2d_ptn_eta1Phi1Eta2Phi2PP+*ptnpp;
			h2d_npt_eta1Phi1Eta2Phi2PM12=h2d_npt_eta1Phi1Eta2Phi2PM12+*nptpm12;
			h2d_npt_eta1Phi1Eta2Phi2PM21=h2d_npt_eta1Phi1Eta2Phi2PM21+*nptpm21;
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
	TFile *ofile;
	if(num==0)
	ofile=new TFile("AnalysisResults_pik.root","recreate");
	if(num==1)
	ofile=new TFile("AnalysisResults_pip.root","recreate");
	if(num==2)
	ofile=new TFile("AnalysisResults_pk.root","recreate");
	ofile->cd();
		h1d_n1_phi.Write();
		h1d_n1_eta.Write();
		h1d_n1_pt.Write();
		h1i_n1_multPM.Write();
		h1d_n1_ptM1.Write();
		h1d_n1_ptP1.Write();
		h2d_n1_etaPhiP1.Write();
		h2d_n1_etaPhiM1.Write();
		h2d_pt_etaPhiP1.Write();
		h2d_pt_etaPhiM1.Write();
		h1d_n1_ptM2.Write();
		h1d_n1_ptP2.Write();
		h2d_n1_etaPhiP2.Write();
		h2d_n1_etaPhiM2.Write();
		h2d_pt_etaPhiP2.Write();
		h2d_pt_etaPhiM2.Write();
		h2d_n2_eta1Phi1Eta2Phi2PP.Write();
		h2d_n2_eta1Phi1Eta2Phi2PM12.Write();
		h2d_n2_eta1Phi1Eta2Phi2PM21.Write();
		h2d_n2_eta1Phi1Eta2Phi2MM.Write();
		h2d_ptpt_eta1Phi1Eta2Phi2PM12.Write();
		h2d_ptpt_eta1Phi1Eta2Phi2PM21.Write();
		h2d_ptpt_eta1Phi1Eta2Phi2MM.Write();
		h2d_ptpt_eta1Phi1Eta2Phi2PP.Write();
		h2d_ptn_eta1Phi1Eta2Phi2PM12.Write();
		h2d_ptn_eta1Phi1Eta2Phi2PM21.Write();
		h2d_ptn_eta1Phi1Eta2Phi2MM.Write();
		h2d_ptn_eta1Phi1Eta2Phi2PP.Write();
		h2d_npt_eta1Phi1Eta2Phi2PM12.Write();
		h2d_npt_eta1Phi1Eta2Phi2PM21.Write();
		h2d_npt_eta1Phi1Eta2Phi2MM.Write();
		h2d_npt_eta1Phi1Eta2Phi2PP.Write();
		delete ofile;
}}
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


*/
