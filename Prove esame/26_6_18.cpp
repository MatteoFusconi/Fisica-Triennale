//QUESITO 1

void q1(){
gRandom->SetSeed();
TH1F* h1 = new TH1F("h1", "h1", 1000, 0., 10.);
TH1F* h2 = new TH1F("h2", "h2", 1000, 0., 10.);
TH1F* h3 = new TH1F("h3", "h3", 1000, 0., 10.);
TF1* eff = new TF1("eff", ".1* x * exp( -x)", 0., 10.);

Float_t x, y;
for(int i = 0; i<1E7; ++i){
  x = gRandom->Gaus(5, 1);
  h1->Fill(x);
  y = gRandom->Rndm();
  if(y < eff->Eval(x)) 
    h2->Fill(x);
}
h3->Sumw2();
h3->Divide(h2, h1, 1, 1, "B");
h3->Draw("E");

}

//QUESITO 2
Double_t func(Double_t* x, Double_t* par) {
Double_t xx = x[0];
Double_t val = par[0]*TMath::Gaus(xx, par[1], par[2]) + par[3]*TMath::Exp(par[4]);
return val;
}
void q2(){
gRandom->SetSeed();
gStyle->SetOptFit(111);
TH1F* h1 = new TH1F("h1", "h1", 500, 0., 5.);
TH1F* h2 = new TH1F("h2", "h2", 500, 0., 5.);
TH1F* h3 = new TH1F("h3", "h3", 500, 0., 5.);

Float_t x;
for(int i = 0; i<1E5; ++i){
  x = gRandom->Gaus(2, .5);
  h1->Fill(x);
  if(i<1E4){
  x = gRandom->Exp(1);
  h2->Fill(x);
}
}
h3->Sumw2();
h3->Add(h1, h2, 1, 1);
TF1* f1 = new TF1("f1", func, 0., 5., 5);
f1->SetParameters(1,1,1,1,1);
h3->Fit(f1, "Q");
h3->Draw();

}

void q3(){
gRandom->SetSeed();
TH1F* h1 = new TH1F("h1", "h1", 100, 0., 10.);
TF1* f1 = new TF1("f1", "sin(x) + x^2", 0., 10.);
h1->FillRandom("f1", 1E5);
h1->Draw();
}






