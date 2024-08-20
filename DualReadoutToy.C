

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
void SCEDraw1_2D_2 (TCanvas* canv, const char* name, TH2F* h1, const char* outfile);
void dotoy(bool doplot,double h_s,double h_c,double nscint,double ncer,double fmean,double frms, double &sssm,double &cccm,double &sigmaS, double &sigmaC, double &sigmaD, double &acov);

int nshowers=10000;


TRandom rrr;


void DualReadoutToy() {

  double h_s=0.9;
  double h_c=0.5;
  double nscint=10000;
  double ncer=10000;
  double fmean=0.6;
  double frms=0.2;
  double sssm,cccm,sigmaS,sigmaC,sigmaD,acov;

  dotoy(1,h_s,h_c,nscint,ncer,fmean,frms,sssm,cccm,sigmaS,sigmaC,sigmaD,acov);
  std::cout<<"mean scint is "<<sssm<<" while predicted is "<<(fmean+(1-fmean)*h_s)<<std::endl;
  std::cout<<"mean cer is "<<cccm<<" while predicted is "<<(fmean+(1-fmean)*h_c)<<std::endl;
  std::cout<<"scint res is "<<sigmaS<<std::endl;
  std::cout<<"cere res is "<<sigmaC<<std::endl;
  double precov=(1-h_s)*(1-h_c)*frms*frms;
  //double precov=(1-h_s)*(1-h_c)*(frms*frms*fmean*fmean/(frms*frms+fmean*fmean));
  std::cout<<"cov is "<<acov<<" while predicted is "<<precov<<std::endl;
  double term1= (1-h_c)*(1-h_c)*sigmaS*sigmaS;
  double term2=(1-h_s)*(1-h_s)*sigmaC*sigmaC;
  double sum12= term1+term2;
  double term3_true= 2*(1-h_s)*(1-h_c)*acov;
  double term3_formula= 2*(1-h_s)*(1-h_c)*precov;
  std::cout<<"1 2 sum 3_true 3_for are "<<term1<<" "<<term2<<" "<<sum12<<" "<<term3_true<<" "<<term3_formula<<std::endl;

  double dualpred = (1/(h_s-h_c))*sqrt(term1+term2-term3_true);
  double dualpred2 = (1/(h_s-h_c))*sqrt(term1+term2-term3_formula);
  std::cout<<"dual res is "<<sigmaD<<" predicted using true cov is  "<<dualpred<<" predicted using formula cov is "<<dualpred2<<std::endl;


  
  TH2F *covcheck = new TH2F("covcheck","covcheck", 1000,0.,0.02,1000,0.,0.02);
  TH2F *covcheckf1 = new TH2F("true-pred vs f","covcheckf1", 1000,0.,0.6,1000,-0.01,0.05);
  TH1F *dualcheck = new TH1F("dualcheck true","dualcheck", 1000,-0.05,0.05);
  TH1F *dualcheckf = new TH1F("dualcheck formula","dualcheckf", 1000,-0.05,0.05);
  TH1F *dualcheckab = new TH1F("dualcheck formulaab","dualcheckab", 1000,-0.05,0.05);
  TH2F *scintdual = new TH2F("scint res vs dual res","scintdual", 1000,0.,0.1,1000,0.,0.1);
  int jmax=0;

  int npts=200;
  double range=min(fmean,1-fmean)/2.;
  for(int j=1;j<npts;j++) {
    double frestry=j*(range/npts);
    std::cout<<"frestry is "<<frestry<<std::endl;
    //    std::cout<<std::endl;
    dotoy(0,h_s,h_c,nscint,ncer,fmean,frestry,sssm,cccm,sigmaS,sigmaC,sigmaD,acov);

    //precov=(1-h_s)*(1-h_c)*(frestry*frestry*fmean*fmean/(frestry*frestry+fmean*fmean));
    covcheck->Fill(acov,precov);
    precov=(1-h_s)*(1-h_c)*frestry*frestry; 
    //    std::cout<<"cov is "<<acov<<" while predicted is "<<precov<<std::endl;

    term1= (1-h_c)*(1-h_c)*sigmaS*sigmaS;
    term2=(1-h_s)*(1-h_s)*sigmaC*sigmaC;
    term3_true= 2*(1-h_s)*(1-h_c)*acov;
    term3_formula= 2*(1-h_s)*(1-h_c)*precov;

    if(term3_formula>term1+term2) {
      std::cout<<"invalid prediction frestry="<<frestry<<std::endl;
      jmax=j-1;
    }
    if(jmax!=0) term3_formula= 2*(1-h_s)*(1-h_c)*(1-h_s)*(1-h_c)*(1./40.)*(1./40.)*jmax*jmax;
    double dualpreda = (1/(h_s-h_c))*sqrt(term1+term2-term3_true);
    double dualpredb = (1/(h_s-h_c))*sqrt(term1+term2-term3_formula);
    dualcheck->Fill(sigmaD-dualpreda);
    dualcheckf->Fill(sigmaD-dualpredb);
    dualcheckab->Fill(dualpreda-dualpredb);
    //    std::cout<<"sigma D and pre "<<sigmaD<<" "<<dualpreda<<std::endl;
    scintdual->Fill(sigmaS,sigmaD);
    covcheckf1->Fill(frestry,sigmaD-dualpredb);
  }
  std::cout<<"ult cov is "<<acov<<std::endl;

  TCanvas* c7;
  SCEDraw1_2D(c7,"c7",covcheck,"junk7.png");
  TCanvas* c7a;
  SCEDraw1_2D(c7a,"c7a",covcheckf1,"junk7a.png");
  TCanvas* c8;
  SCEDraw1(c8,"c8",dualcheck,"junk8.png",0);
  TCanvas* c9;
  SCEDraw1(c9,"c9",dualcheckf,"junk9.png",0);
  TCanvas* c91;
  SCEDraw1(c91,"c91",dualcheckab,"junk9a.png",0);
  TCanvas* c10;
  SCEDraw1_2D_2(c10,"c10",scintdual,"junk10.png");


  TH2F *scintdual2 = new TH2F("scint2 res vs dual res","scintdual2", 1000,0.,0.1,1000,0.,0.1);
  TH2F *dualcheckha = new TH2F("dual true v p ha","dualcheckha", 1000,0.,0.1,1000,0.,0.1);
  for(int j=1;j<500;j++) {
    double nnn=100.*j;

    dotoy(0,h_s,h_c,nnn,nnn,fmean,frms,sssm,cccm,sigmaS,sigmaC,sigmaD,acov);
    scintdual2->Fill(sigmaS,sigmaD);

    precov=(1-h_s)*(1-h_c)*(frms*frms*fmean*fmean/(frms*frms+fmean*fmean));
    term1= (1-h_c)*(1-h_c)*sigmaS*sigmaS;
    term2=(1-h_s)*(1-h_s)*sigmaC*sigmaC;
    term3_true= 2*(1-h_s)*(1-h_c)*acov;
    term3_formula= 2*(1-h_s)*(1-h_c)*precov;
    double dualpreda = (1/(h_s-h_c))*sqrt(term1+term2-term3_true);
    double dualpredb = (1/(h_s-h_c))*sqrt(term1+term2-term3_formula);
    dualcheckha->Fill(dualpreda,dualpredb);

  }
  std::cout<<"last dualres is "<<sigmaD<<std::endl;
  TCanvas* c11;
  SCEDraw1_2D_2(c11,"c11",scintdual2,"junk11.png");
  TCanvas* c11a;
  SCEDraw1_2D_2(c11a,"c11a",dualcheckha,"junk11a.png");



  TH2F *scintdual3 = new TH2F("cer mean vs dual res","scintdual3", 1000,0.,1.0,1000,0.,0.1);
  for(int j=1;j<10;j++) {
    double fmeanaa=0.2+0.05*j;

    dotoy(0,h_s,h_c,nscint,ncer,fmeanaa,frms,sssm,cccm,sigmaS,sigmaC,sigmaD,acov);
    scintdual3->Fill(cccm,sigmaD);

  }
  TCanvas* c12;
  SCEDraw1_2D(c12,"c12",scintdual3,"junk12.png");

  int ijunk=0;
  std::cout<<"input an integer"<<std::endl;
  std::cin>>ijunk;

  
}


void dotoy(bool doplot, double h_s,double h_c,double nscint,double ncer,double fmean,double frms, double &sssm, double& cccm, double &sigmaS, double &sigmaC, double &sigmaD, double &acov) {

  TH1F *fff = new TH1F("fff","shower em fraction",300,0.,2.0);
  TH1F *sss = new TH1F("sss","shower scintillation",300,0.,2.0);
  TH1F *ccc = new TH1F("ccc","shower cherenkov",300,0.,2.0);
  TH2F *sscc = new TH2F("sscc","cheren versus scint", 1000,0.,2.0,1000,0.,2.0);
  TH1F *ddd = new TH1F("ddd","dual readout",900,0.,2.0);
  TH1F *cov = new TH1F("cov","covariance",3000,-2.,2.0);
  fff->Reset();
  sss->Reset();
  ccc->Reset();
  sscc->Reset();
  ddd->Reset();
  cov->Reset();

  //  std::cout<<"frms is "<<frms<<std::endl;


  acov=0.;
  double pmeans=(fmean+(1-fmean)*h_s);
  double pmeanc=(fmean+(1-fmean)*h_c);
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
    double CCC=(FFF+h_c*(1-FFF));    
    CCC=CCC*(1+rrr.Gaus(0.,1/sqrt(ncer)));
    ccc->Fill(CCC);
    sscc->Fill(SSS,CCC);

    double DDD=((1-h_c)*SSS - (1-h_s)*CCC)/(h_s-h_c);
    ddd->Fill(DDD);

    cov->Fill((SSS-pmeans)*(CCC-(fmean-pmeanc)));
    acov+=(SSS-pmeans)*(CCC-pmeanc);
  }
  acov=acov/(nshowers-1);

  sigmaS=sss->GetRMS();
  sigmaC=ccc->GetRMS();
  sigmaD=ddd->GetRMS();
  sssm=sss->GetMean();
  cccm=ccc->GetMean();
  //  std::cout<<"sigma S C D "<<sigmaS<<" "<<sigmaC<<" "<<sigmaD<<std::endl;

  //double covmean = cov->GetMean();
  //int iii=cov->GetMaximumBin();
  //double covmean=cov->GetBinCenter(iii);
 
  if(doplot) {
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

  delete fff;
  delete sss;
  delete ccc;
  delete sscc;
  delete ddd;
  delete cov;

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
  //canv->Update();

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
  //canv->Update();

  return;
}

void SCEDraw1_2D_2 (TCanvas* canv,  const char* name,TH2F* h1, const char* outfile) {

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


  TLine line = TLine(0.,0.,1.,1.);
  line.SetLineColor(kYellow);
  line.SetLineWidth(2);
  line.Draw();


  canv->cd(0);
  canv->Modified();
  canv->Update();
canv->Print(outfile,".png");
  return;
}
