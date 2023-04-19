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

Double_t fitf(Double_t *x, Double_t *par){
  Double_t fitval = par[0] + par[1] * x[0];
// p1=b p0=a a+bx
  if (2.265 < x[0] && x[0] < 2.305) fitval = 0;
  return fitval;
}

Double_t fun2(Double_t *x, Double_t *par){
  Double_t fitval = par[0] + par[1] * x[0];
// p1=b p0=a a+bx
  return fitval;
}
Double_t mygaus(Double_t* x, Double_t* par)
{
  Double_t total = par[0] * (TMath::Gaus(x[0], par[1], par[2], 0));
  return total;
}


void FitRot(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111);
  gStyle->SetOptStat(10);

  Int_t rebin = 12;
  Int_t nbins = 1000/rebin;

  Float_t min = 2.145;
  Float_t max = 2.425;

  TFile *f = new TFile("TMVAApp_BDT_SigmacPt_20220329_0_1_RotationalBackground.root");
  TH2F *hr2 =(TH2F*)f->Get("MVA_BDT_vs_InvMass");
  hr2->GetXaxis()->SetRangeUser(0,1);
  TH1F *hRotBG = (TH1F*) hr2-> ProjectionY("hInvMass");
  hRotBG->Rebin(rebin);
  
  TFile* fi = new TFile("TMVAApp_BDT_SigmacPt_20220504_0_1.root");
  TH2F * h2 = (TH2F*)fi->Get("MVA_BDT_vs_InvMass");
  h2->GetXaxis()->SetRangeUser(0,1);
  TH1F *hInvMass = (TH1F*) h2-> ProjectionY("hInvMass");
  hInvMass->Rebin(rebin);




  TH1F* hRatio = (TH1F*)hInvMass->Clone("hRatio");
  hRatio->SetTitle("Rapporto tra dati reali e fondo rotazionale");
  hRatio->Divide(hRotBG);

  hRatio->GetXaxis()->SetTitle("m_{inv}(pK^{0}_{S})[GeV/#it{c}^{2}]");
/*  hRatio->Draw();
  */ 

Float_t a = 0.;
Int_t b = 0;


  TF1* f1 = new TF1("f1", fitf, min, max, 2);
  f1->SetParNames("a","b");

  hRatio->Fit(f1,"r");
   Double_t p[2];
  f1->GetParameters(p);

  TF1* f2 = new TF1("flin", fun2, 2.05, 2.5, 2);
  f2->SetParameters(p);

  TH1F *hBackground = (TH1F*)hInvMass->Clone("hBackground");
  hBackground->SetTitle("Fondo riscalato");
  TH1F *hSignal = (TH1F*)hInvMass->Clone("hSignal");
  hSignal->SetTitle("Segnale Massa invariante");

  //TH1F* hBackground = new TH1F("hBackground","Fondo riscalato", nbins, 2.05, 2.5);
  for(int i = 0; i< nbins; ++i){
   a = hInvMass->GetBinCenter(i);
   b = static_cast<Int_t> (f2->Eval(a)*hRotBG->GetBinContent(i));
   hBackground-> SetBinContent(i,b);

  }

  //TH1F* hSignal = new TH1F("hSignal", "Segnale di Massa invariante", nbins, 2.05, 2.5);
  for(int i = 0; i< nbins; ++i){
    a = hInvMass->GetBinCenter(i);
    if(a>min && a<max){
      b = hInvMass->GetBinContent(i) - hBackground->GetBinContent(i);
      hSignal->SetBinContent(i, b);
      hSignal->SetBinError(i, hInvMass->GetBinError(i));
    }
    else{ hSignal->SetBinContent(i, 0);
          hSignal->SetBinError(i,0);
    }
  }
  
  TF1* fitS = new TF1("fitS", mygaus, min, max, 3);
  fitS->SetName("fitS");
  fitS->SetParName(0, "const");
  fitS->SetParameter(0, 5.0);
  fitS->SetParName(1, "Mass");
  fitS->SetParameter(1, 2.2865);
  fitS->SetParLimits(1, 2.2, 2.3);
  fitS->SetParName(2, "Sigma");
  fitS->FixParameter(2, 0.0076);
  fitS->SetLineColor(kRed);

  hSignal->Fit("fitS","","ep",min,max);
  
  TCanvas *c2 = new TCanvas("c2");
  gPad->SetGrid(1,1);

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



  fitS->SetLineWidth(3);
  hSignal->Draw();
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerColor(kBlack);
  hSignal->SetLineColor(kBlack);
  fitS->Draw("same");

  Float_t binwidth = hSignal->GetBinWidth(1);
  Float_t signal = fitS->Integral(min,max)/binwidth; //N di lambda c
  Float_t errSignal = fitS->IntegralError(min,max)/binwidth;

  std::cout<<"Signal = " <<signal <<"+-"<<errSignal<<std::endl;
//*/<<std::setprecision(4)
}