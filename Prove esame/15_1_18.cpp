//QUESITO 1
void q1(){
TH1F* h1 = new TH1F ("h1", "h1", 1000, 0., 5.);
TH1F* h2 = new TH1F ("h2", "h2", 1000, 0., 5.);
TH1F* hEff = new TH1F ("hEff", "hEff", 1000, 0., 5.);

Float_t x, y;
TF1* eff = new TF1("eff", "x/5", 0., 5.);

for (int i = 0; i < 1E7; ++i){
  x = gRandom->Exp(1);
  y = gRandom->Rndm();
  h1->Fill(x);
  if(y < eff->Eval(x)) h2->Fill(x);
}
hEff->Sumw2();
hEff->Divide(h2, h1, 1, 1,"B");

hEff->Draw("E");
}

//QUESITO 2
Double_t func(Double_t* x, Double_t* par){
  Double_t xx= x[0];
  Double_t val = par[0]*(TMath::Gaus(xx, par[1],par[2])) + par[3];
  return val;
}

void q2(){
gStyle->SetOptFit(111);
gRandom->SetSeed();
TH1F* h1 = new TH1F("h1", "h1", 500, 0., 5.);
TH1F* h2 = new TH1F("h2", "h2", 500, 0., 5.);
TH1F* hSum = new TH1F("hSum", "hSum", 500, 0., 5.);
Float_t x, y;
for(Int_t i = 0; i < 1E6; ++i){
  x = gRandom->Gaus(2.5, .25);
  h1->Fill(x);
  if(i<1E4){
    y = gRandom->Uniform(0, 5);
    h2->Fill(y);
  }
}
hSum->Sumw2();
hSum->Add(h1, h2, 1, 1);
TF1* f1 = new TF1("f1", func, 0., 5., 4);
//TF1* f1 = new TF1("f1", "[0]*(Gaus(x, [1], [2])) + [3]", 0., 5., 4);
f1->SetParameters(1,1,1,1);
hSum->Fit(f1, "Q");
for (int i = 0; i<4; ++i){
  cout << f1->GetParameter(i) << "+/-" << f1->GetParError(i) << endl;
  }
cout<< "Chisquare = " <<   f1->GetChisquare()/f1->GetNDF()<<endl;
hSum->Draw();
}

//QUESITO 3

void q3(){
gRandom->SetSeed();

TH1F* h1 = new TH1F("h1", "h1", 100, 0., 10.);
TF1* f1 = new TF1("f1", "sqrt(x)+ x**2", 0., 10.);
h1->FillRandom("f1", 1E7);
h1->Draw();
}
