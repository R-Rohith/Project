
#include<iostream>
#include<string.h>

#include<TFile.h>
#include<TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
void PID_efficiency(const char *arg)
{
TFile *infile;
infile=new TFile("AnalysisResults_MC.root","read");
TH1F *purity1=new TH1F("Purity1","PID Purity;p_{T}",100,0,2),*efficiency1=new TH1F("Efficiency1","PID Efficiency;p_{T}",100,0,2),*purity2,*efficiency2;
TH1I *curr,*pure,*reco;
float num,den,part;
char str[50];
if((strcmp(arg,"pi")==0)||(strcmp(arg,"ka")==0))
    part=0.6;
else if(strcmp(arg,"pr")==0)
    part=1.1;
else return;
auto c1=new TCanvas("c1","c1",900,600);
sprintf(str,"fill-flags-table/pt%s",arg);
curr=(TH1I*)infile->Get(str);
sprintf(str,"fill-flags-table/pureidpt%s",arg);
pure=(TH1I*)infile->Get(str);
sprintf(str,"fill-flags-table/recopt%s",arg);
reco=(TH1I*)infile->Get(str);
for(int i=1;i<101;i++){
    den=curr->GetBinContent(i);
    num=pure->GetBinContent(i);
//    cout<<num<<" -- "<<den<<endl;
    if(den==0)
        num=0;
    else
        num=num/den;
//    cout<<num<<" "<<i<<endl;
    purity1->SetBinContent(i,num);
}
purity1->Draw();
sprintf(str,"MCPurity_%s.pdf",arg);
c1->SaveAs(str);
purity2=(TH1F*)purity1->Clone("Purity2");
purity1->GetXaxis()->SetRangeUser(0.2,part);
purity2->GetXaxis()->SetRangeUser(part,2.0);
purity1->Draw();
sprintf(str,"MCPurity1_%s.pdf",arg);
c1->SaveAs(str);
purity2->Draw();
sprintf(str,"MCPurity2_%s.pdf",arg);
c1->SaveAs(str);
for(int i=1;i<101;i++){
    den=reco->GetBinContent(i);
    num=curr->GetBinContent(i);
//    cout<<num<<" -- "<<den<<endl;
    if(den==0)
        num=0;
    else
        num=num/den;
//    cout<<num<<" "<<i<<endl;
    efficiency1->SetBinContent(i,num);
}
efficiency1->Draw();
sprintf(str,"MCEfficiency_%s.pdf",arg);
c1->SaveAs(str);
efficiency2=(TH1F*)efficiency1->Clone("Efficiency2");
efficiency1->GetXaxis()->SetRangeUser(0.2,part);
efficiency2->GetXaxis()->SetRangeUser(part,2.0);
efficiency1->Draw();
sprintf(str,"MCEfficiency1_%s.pdf",arg);
c1->SaveAs(str);
efficiency2->Draw();
sprintf(str,"MCEfficiency2_%s.pdf",arg);
c1->SaveAs(str);
}
