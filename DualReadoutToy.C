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
void dotoy(bool doplot,double h_s,double h_c,double nscint,double ncer,double fmean,double frms, TH1F* fff, TH1F* sss, TH1F* ccc, TH2F* sscc, TH1F* ddd, TH1F* cov,double &sssm,double &cccm,double &sigmaS, double &sigmaC, double &sigmaD, double &acov);

int nshowers=10000;


TRandom rrr;


void DualReadoutToy() {
  TH1F *fff = new TH1F("fff","shower em fraction",300,0.,2.0);
  TH1F *sss = new TH1F("sss","shower scintillation",300,0.,2.0);
  TH1F *ccc = new TH1F("ccc","shower cherenkov",300,0.,2.0);
  TH2F *sscc = new TH2F("sscc","cheren versus scint", 1000,0.,2.0,1000,0.,2.0);
  TH1F *ddd = new TH1F("ddd","dual readout",300,0.,2.0);
  TH1F *cov = new TH1F("cov","covariance",300,-2.,2.0);
  double h_s=0.9;
  double h_c=0.3;
  double nscint=1000;
  double ncer=1000;
  double fmean=0.5;
  double frms=0.2;
  bool doplot=1;
  double sssm,cccm,sigmaS,sigmaC,sigmaD,acov;

  dotoy(doplot,h_s,h_c,nscint,ncer,fmean,frms,fff,sss,ccc,sscc,ddd,cov,sssm,cccm,sigmaS,sigmaC,sigmaD,acov);
  std::cout<<"mean scint is "<<sssm<<" while predicted is "<<(fmean+(1-fmean)*h_s)<<std::endl;
  std::cout<<"mean cer is "<<cccm<<" while predicted is "<<(fmean+(1-fmean)*h_c)<<std::endl;
  std::cout<<"scint res is "<<sigmaS<<std::endl;
  std::cout<<"cere res is "<<sigmaC<<std::endl;
  double HHH=(1-h_s)*(1-h_c);
  double precov=HHH*frms*frms;
  std::cout<<"cov is "<<acov<<" while predicted is "<<precov<<std::endl;
  double term1= (1-h_c)*(1-h_c)*sigmaS*sigmaS;
  double term2=(1-h_s)*(1-h_s)*sigmaC*sigmaC;
  double sum12= term1+term2;
  double term3_true= 2*(1-h_s)*(1-h_c)*acov;
  double term3_formula= 2*(1-h_s)*(1-h_c)*precov;
  std::cout<<"1 2 sum 3_true are "<<term1<<" "<<term2<<" "<<sum12<<" "<<term3_true<<std::endl;

  double dualpred = (1/(h_s-h_c))*sqrt(term1+term2-term3_true);
  double dualpred2 = (1/(h_s-h_c))*sqrt(term1+term2-term3_formula);
  std::cout<<"dual res is "<<sigmaD<<" predicted using true cov is  "<<dualpred<<" predicted using formula cov is "<<dualpred2<<std::endl;



  TH2F *covcheck = new TH2F("covcheck","covcheck", 1000,0.,0.01,1000,0.,0.01);
  doplot=0;
  for(int j=1;j<100;j++) {
    double frestry=(1/100.)*j;
    dotoy(doplot,h_s,h_c,nscint,ncer,fmean,frestry,fff,sss,ccc,sscc,ddd,cov,sssm,cccm,sigmaS,sigmaC,sigmaD,acov);
  double HHH=(1-h_s)*(1-h_c);
  double precov=HHH*frestry*frestry; 
  covcheck->Fill(acov,precov);
  std::cout<<"cov is "<<acov<<" while predicted is "<<precov<<std::endl;
  }


  TCanvas* c1;
  SCEDraw1(c1,"c1",fff,"junk1.png",0);
  TCanvas* c2;
  SCEDraw1(c2,"c2",sss,"junk2.png",0);
  TCanvas* c3;
  SCEDraw1(c3,"c3",ccc,"junk3.png",0);
  TCanvas* c4;
  SCEDraw1_2D(c4,"c4",sscc,"junk4.png");
  TCanvas* c5;
  SCEDraw1(c5,"c5",ddd,"junk5.png",0);
  TCanvas* c6;
  SCEDraw1(c6,"c6",cov,"junk6.png",0);
  TCanvas* c7;
  SCEDraw1_2D(c7,"c7",covcheck,"junk7.png");

}


void dotoy(bool doplot,double h_s,double h_c,double nscint,double ncer,double fmean,double frms, TH1F* fff, TH1F* sss, TH1F* ccc, TH2F* sscc, TH1F* ddd, TH1F* cov,double &sssm, double& cccm, double &sigmaS, double &sigmaC, double &sigmaD, double &acov) {


  double HHH=(1-h_s)*(1-h_c);
  acov=0.;
  double pmeans=(fmean+(1-fmean)*h_s);
  double pmeanc=(fmean+(1-fmean)*h_c);
  for (int i=0;i<nshowers;i++) {
    // pick a value for fraction EM in the shower
    double FFF=-1.;
    while((FFF<0)||(FFF>1) ) {
      FFF = rrr.Gaus(fmean,frms);
    }
    if(doplot) fff->Fill(FFF);
    double SSS=(FFF+h_s*(1-FFF));
    SSS=SSS*(1+rrr.Gaus(0.,1/sqrt(nscint)));
    if(doplot) sss->Fill(SSS);
    double CCC=(FFF+h_c*(1-FFF));    
    CCC=CCC*(1+rrr.Gaus(0.,1/sqrt(ncer)));
    if(doplot) ccc->Fill(CCC);
    if(doplot) sscc->Fill(SSS,CCC);

    double DDD=((1-h_c)*SSS - (1-h_s)*CCC)/(h_s-h_c);
    if(doplot) ddd->Fill(DDD);

    if(doplot) cov->Fill((SSS-pmeans)*(CCC-(fmean-pmeanc)));
    acov+=(SSS-pmeans)*(CCC-pmeanc);
  }
  acov=acov/(nshowers-1);

  sigmaS=sss->GetRMS();
  sigmaC=ccc->GetRMS();
  sigmaD=ddd->GetRMS();
  sssm=sss->GetMean();
  cccm=ccc->GetMean();


  //double covmean = cov->GetMean();
  //int iii=cov->GetMaximumBin();
  //double covmean=cov->GetBinCenter(iii);
 
 

  return;
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
