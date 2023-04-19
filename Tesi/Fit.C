#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <THnSparse.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TRandom.h"

Double_t ffitf(Double_t *x, Double_t *par){
  Double_t fitval = par[3] + par[2] * x[0] + par[1] * x[0]*x[0] + par[0] * x[0] * x[0] * x[0];
//p3= d p2=c p1=b p0=a ax3+bx2+cx+d
  if (2.265 < x[0] && x[0] < 2.305) fitval = 0;
  return fitval;

}
Double_t ffun2(Double_t *x, Double_t *par){
  Double_t val = par[3] + par[2] * x[0] + par[1] * x[0]*x[0] + par[0] * x[0] * x[0] * x[0];
//p3= d p2=c p1=b p0=a ax3+bx2+cx+d

  return val;

}
void Fit(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111);
  gStyle->SetOptStat(10);

  Int_t rebin = 12;
  Int_t nbins = 1000/rebin;

  //fit limits
  Float_t min = 2.145;
  Float_t max = 2.425;
  
  TFile *f = new TFile("TMVAApp_BDT_SigmacPt_20220504_0_1.root");
  TH2F *h2 =(TH2F*)f->Get("MVA_BDT_vs_InvMass");
  h2->GetXaxis()->SetRangeUser(0,1);
  TH1F *hInvMass = (TH1F*) h2-> ProjectionY("hInvMass");
  hInvMass->Rebin(rebin);

  TF1 *f1 = new TF1("ffitf",ffitf,min,max,4);
  //func->Draw();
  f1->SetParameters(2.45632e+05,-1.71793e+06,3.99246e+06,-3.07494e+06); 
  f1->SetParNames("a","b","c","d");//
  hInvMass->Fit(f1,"r");
  //hInvMass->Draw();
  Double_t p[4];
  f1->GetParameters(p);
  TF1 * f2 = new TF1("pol3", ffun2 ,2.05,2.55, 4 );
  f2->SetParameters(p);
  //f2->Draw();

//fit range [2.08675,2.486]
  hInvMass->GetXaxis()->SetRangeUser(2.05,2.5);
  TH1F* hNoise = new TH1F ("hNoise", "Fondo", nbins, 2.05, 2.55);
Float_t a = 0.;
Int_t b = 0;
  for(int i = 0; i < nbins; ++i){
    a = hInvMass->GetBinCenter(i);
    b = static_cast<Int_t> (f2->Eval(a));
    hNoise->SetBinContent(i, b);
  }
//hNoise->Draw();

  TH1F* hSignal = new TH1F ("hSignal","Segnale Massa Invariante", nbins, 2.05, 2.55);

for(int i = 0; i<nbins; ++i){
  a = hInvMass->GetBinCenter(i);
  if (a>min && a<max){ // hInvMass->GetBinContent(i)!=0
  b = hInvMass->GetBinContent(i)-hNoise->GetBinContent(i);
  hSignal->SetBinContent(i, b);
  hSignal->SetBinError(i,hInvMass->GetBinError(i));
  }
}
  //hSignal->Add(hNoise, hInvMass, -1, 1);
  hSignal->Sumw2();
  //hSignal->Fit("gaus");
  TF1* fitS = new TF1("fitS", "gaus", min,max);
  fitS->SetLineColor(kRed);
  fitS->SetParameter(1,2.287);
  fitS->FixParameter(2,0.0076);
  fitS->SetParLimits(1,2.2,2.35);
  

 //TCanvas *c1 = new TCanvas("c1");
  gPad->SetGrid(1,1);

  hSignal->Fit("fitS","","ep",min,max);
  
  hSignal->GetXaxis()->SetTitle("m_{inv}(pK^{0}_{S})[GeV/#it{c}^{2}]");
  hSignal->GetYaxis()->SetTitle("Entries/6.0 MeV/#it{c}^{2}");
  hSignal->GetYaxis()->SetLabelFont(42);
  hSignal->GetXaxis()->SetTitleFont(42);
  hSignal->GetYaxis()->SetTitleFont(42);
  hSignal->GetYaxis()->SetTitleOffset(1.25);
  hSignal->GetXaxis()->SetRangeUser(min,max);
  fitS->SetLineWidth(3);
  hSignal->Draw();
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerColor(kBlack);
  hSignal->SetLineColor(kBlack);
  fitS->Draw("same");
  
   
  //*/
  Float_t binwidth = hSignal->GetBinWidth(1);
  Float_t signal = fitS->Integral(min,max)/binwidth; //N di lambda c
  Float_t errSignal = fitS->IntegralError(min,max)/binwidth;

  std::cout<<"Signal = " <<signal <<"+-"<<errSignal<<std::endl;
  


//*/std::setprecision(4)<<
}
