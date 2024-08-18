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


int nshowers=10000;
double h_s=0.9;
double h_c=0.6;
double nscint=10;
double ncer=10;
double fmean=0.5;
double frms=0.03;


void DualReadoutToy() {

  TRandom rrr;
  TH1F *fff = new TH1F("fff","shower em fraction",300,0.,2.0);
  TH1F *sss = new TH1F("sss","shower scintillation",300,0.,2.0);
  TH1F *ccc = new TH1F("ccc","shower cherenkov",300,0.,2.0);
  TH2F *sscc = new TH2F("sscc","cheren versus scint", 1000,0.,2.0,1000,0.,2.0);
  TH1F *ddd = new TH1F("ddd","dual readout",300,0.,2.0);
  TH1F *cov = new TH1F("cov","covariance",300,0.,2.0);


  double acov=0.;
  for (int i=0;i<nshowers;i++) {
    // pick a value for fraction EM in the shower
    double FFF=-1.;
 
    while((FFF<0)||(FFF>1) ) {
      FFF = rrr.Gaus(fmean,frms);
    }
    fff->Fill(FFF);
    double SSS=(FFF+h_s*(1-FFF));
    SSS=SSS*(1+rrr.Gaus(0.,1/sqrt(nscint)));
    sss->Fill(SSS);
    double CCC=(FFF+h_c*(1-FFF));    CCC=CCC*(1+rrr.Gaus(0.,1/sqrt(ncer)));
    ccc->Fill(CCC);
    sscc->Fill(SSS,CCC);
    //std::cout<<"FFF SSS CCC are "<<FFF<<" "<<SSS<<" "<<CCC<<std::endl;

    double DDD=((1-h_c)*SSS - (1-h_s)*CCC)/(h_s-h_c);
    ddd->Fill(DDD);
    //std::cout<<"DDD is "<<DDD<<std::endl;

    cov->Fill((SSS-(fmean-(1-fmean)*h_s))*(CCC-(fmean-(1-fmean)*h_c)));
    acov+=(SSS-(fmean-(1-fmean)*h_s))*(CCC-(fmean-(1-fmean)*h_c));
  }
  acov=acov/(nshowers-1);
  std::cout<<"acov is "<<acov<<" mean of cov is "<<cov->GetMean()<<" rms of cov is "<<cov->GetRMS()<<std::endl;


  double sigmaS=sss->GetRMS();
  double sigmaC=ccc->GetRMS();
  double sigmaD=ddd->GetRMS();
  std::cout<<"scint res is "<<sigmaS<<std::endl;
  std::cout<<"cere res is "<<sigmaC<<std::endl;
  std::cout<<"dual res is "<<sigmaD<<std::endl;


  std::cout<<" mean scint is "<<sss->GetMean()<<" while predicted is "<<(fmean+(1-fmean)*h_s)<<std::endl;
  std::cout<<" mean cer is "<<ccc->GetMean()<<" while predicted is "<<(fmean+(1-fmean)*h_c)<<std::endl;
  double covmean = cov->GetMean();
  int iii=cov->GetMaximumBin();
  covmean=cov->GetBinCenter(iii);
  double precov=h_c*(1-h_s)+h_s*(1-h_c)+(1/sqrt(2))*(1-h_s)*(1-h_c);
  std::cout<<" mean cov is "<<covmean<<" while predicted is "<<precov<<std::endl;
  std::cout<<"ratio is "<<covmean/precov<<std::endl;


  double term1= (1-h_c)*(1-h_c)*sigmaS*sigmaS;
  double term2=(1-h_s)*(1-h_s)*sigmaC*sigmaC;
  double sum12= term1+term2;
  double term3= 2*(1-h_s)*(1-h_c)*acov;
  std::cout<<"1 2 sum 3 are "<<term1<<" "<<term2<<" "<<sum12<<" "<<term3<<std::endl;
 
  double dualpred = (1/(h_s-h_c))*sqrt(term1+term2-term3);


  std::cout<<"predicted dual resolution "<<dualpred<<std::endl;

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
