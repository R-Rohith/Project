#include<iostream>

#include<TFile.h>
#include<TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

void ExtractPlots()
{
    char input_filename[]="R2P2B2_prc_.root";
    TFile *infile=new TFile(input_filename,"read");
    TH2D *Thist2d;
    TH1D *projx,*projy;
    const char histname[8][30]={"h2d_r2_DetaDphiUS","h2d_r2_DetaDphiLS","h2d_r2CI_DetaDphi","h2d_r2CD_DetaDphi","h2d_p2DptDpt_DetaDphiUS","h2d_p2DptDpt_DetaDphiLS","h2d_p2DptDptCI_DetaDphi","h2d_p2DptDptCD_DetaDphi"};
    auto c1=new TCanvas("c1","c1",900,600);//use to specify dimensions
    char stringbffr[50];
    for(int i=0;i<8;i++)//revert to 0-8
    {
        Thist2d=(TH2D*)infile->Get(histname[i]);
        projx=Thist2d->ProjectionX();
        projy=Thist2d->ProjectionY();
        Thist2d->Draw("SURF3Z");// SPEC a(1,44,180)");
//        Thist2d->Draw("SURF3 SAME");
        sprintf(stringbffr,"%s.pdf",histname[i]);
        c1->SaveAs(stringbffr);
        projx->Draw("E");
        sprintf(stringbffr,"%s_projX.pdf",histname[i]);
        c1->SaveAs(stringbffr);
        projy->Draw("E");
        sprintf(stringbffr,"%s_projY.pdf",histname[i]);
        c1->SaveAs(stringbffr);
    }
    delete infile;
}//TPad::SaveAs
