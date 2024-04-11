#include<iostream>
#include<TFile.h>
#include<TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include<TF1.h>
#include<TMath.h>
#include<unistd.h>
void BB()
{
    TFile *infile=new TFile("BB&beta.root","read");
    TH2D *beta,*BB;
    auto canv=new TCanvas("canv","canv",900,600);
    beta=(TH2D*)infile->Get("r2/Beta");
    beta->Draw();
    canv->SetLogz();
auto betafn1=new TF1("beta","(x/TMath::Sqrt(0.000000261+x*x))",0.2,2.0);
betafn1->SetLineColor(6);
betafn1->SetTitle("electron");
auto betafn2=new TF1("beta","(x/TMath::Sqrt(0.01948+x*x))",0.2,2.0);
betafn2->SetLineColor(2);
betafn2->SetTitle("pion");
auto betafn3=new TF1("beta","(x/TMath::Sqrt(0.2+x*x))",0.2,2.0);
betafn3->SetLineColor(3);
betafn3->SetTitle("kaon");
auto betafn4=new TF1("beta","(x/TMath::Sqrt(0.8+x*x))",0.2,2.0);
betafn4->SetLineColor(4);
betafn4->SetTitle("proton");
betafn1->Draw("Same");
betafn2->Draw("Same");
betafn3->Draw("Same");
betafn4->Draw("Same");
gPad->BuildLegend();
canv->Update();
//sleep(3);
return;
BB=(TH2D*)infile->Get("r2/BB");
BB->Draw();
double c1=3,c2=1000000;
auto bbfn1=new TF1("beta","[0]*((1+(0.000000261/(x*x)))*(TMath::Log(([1]*x*x)/0.000000261))-1)",0.2,2.0);
bbfn1->SetParameters(c1,c2);
bbfn1->SetLineColor(6);
bbfn1->SetTitle("electron");
auto bbfn2=new TF1("bb","[0]*((1+(0.01948/(x*x)))*(TMath::Log(([1]*x*x)/0.01948))-1)",0.2,2.0);
bbfn2->SetParameters(c1,c2);
bbfn2->SetLineColor(2);
bbfn2->SetTitle("pion");
auto bbfn3=new TF1("bb","[0]*((1+(0.244/(x*x)))*(TMath::Log(([1]*x*x)/0.244))-1)",0.2,2.0);
bbfn3->SetParameters(c1,c2);
bbfn3->SetLineColor(3);
bbfn3->SetTitle("kaon");
auto bbfn4=new TF1("bb","[0]*((1+(0.880354/(x*x)))*(TMath::Log(([1]*x*x)/0.880354))-1)",0.2,2.0);
bbfn4->SetParameters(c1,c2);
bbfn4->SetLineColor(4);
bbfn4->SetTitle("proton");
bbfn1->Draw("Same");
bbfn2->Draw("Same");
bbfn3->Draw("Same");
bbfn4->Draw("Same");
gPad->BuildLegend();
canv->Update();
}
