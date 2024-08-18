// force update
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

#include <vector>
#include <map>
#include <algorithm>


void SCEDraw1 (TCanvas* canv, const char* name, TH1F* h1, const char* outfile, bool logy);
void SCEDraw1_2D (TCanvas* canv, const char* name, TH2F* h1, const char* outfile);


int nshowers=10;
double h_s=0.9;
double h_c=0.7;



void DualReadoutToy() {

  TRandom rrr;
  TH1F *fff = new TH1F("fff","shower em fraction",100,0.,1.);
  TH1F *sss = new TH1F("sss","shower scintillation",100,0.,1.);
  TH1F *ccc = new TH1F("ccc","shower cherenkov",100,0.,1.);
  TH2F *sscc = new TH2F("sscc","cheren versus scint", 100,0.,1.,100,0.,1.);

  for (int i=0;i<nshowers;i++) {
    // pick a value for fraction EM in the shower
    double a=-1.;
    while((a<0)||(a>1) ) {
      a = rrr.Gaus(0.7,0.3);
    }
    fff->Fill(a);
    double SSS=(1+h_s*(1-a));
    sss->Fill(SSS);
    double CCC=(1+h_s*(1-a));
    ccc->Fill(CCC);
    sscc->Fill(SSS,CCC);
  }


  TCanvas* c1;
  SCEDraw1(c1,"c1",fff,"junk1.png",0);
  TCanvas* c2;
  SCEDraw1(c2,"c2",sss,"junk2.png",0);
  TCanvas* c3;
  SCEDraw1(c3,"c3",ccc,"junk3.png",0);
  TCanvas* c4;
  SCEDraw1_2D(c4,"c4",sscc,"junk4.png");

}



void SCEDraw1 (TCanvas* canv,  const char* name,TH1F* h1, const char* outfile, bool logy) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);
  if(logy) canv->SetLogy();

  h1->SetLineColor(kGreen);
  h1->SetLineWidth(1);
  h1->SetStats(111111);
  h1->Draw("HIST");


  canv->Print(outfile,".png");
  canv->Update();

  return;
}



void SCEDraw1_2D (TCanvas* canv,  const char* name,TH2F* h1, const char* outfile) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);

  h1->SetLineColor(kGreen);
  h1->SetLineWidth(kGreen);
  h1->SetStats(111111);
  h1->Draw("");

  canv->Print(outfile,".png");
  canv->Update();

  return;
}
