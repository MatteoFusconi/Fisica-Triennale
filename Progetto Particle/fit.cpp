
void fit(){

/////////////////////////////////
////Cosmetics////////////////////
/////////////////////////////////
gROOT->SetStyle("Plain");
gStyle->SetOptFit(1111);
gStyle->SetOptStat(2210);

/////////////////////////////////
////Get Histograms from file/////
/////////////////////////////////
 TFile * file = new TFile("lab2.root", "READ");
 
 TH1F * hfTypes = (TH1F*)file->Get("HistTypes");
 TH1F * hfPhi = (TH1F*)file->Get("HistPhi");
 TH1F * hfTheta = (TH1F*)file->Get("HistTheta");
 TH1F * hfP = (TH1F*)file->Get("HistP"); 
 
 TH1F * h1 = (TH1F*)file->Get("HistIM2"); 
 TH1F * h2 = (TH1F*)file->Get("HistIM1");
 TH1F * h3 = (TH1F*)file->Get("HistIM3");
 TH1F * h4 = (TH1F*)file->Get("HistIM4");
 TH1F * h5 = (TH1F*)file->Get("HistIMDecay");
 
 /////////////////////////////////
 ////Print Particles Abundances///
 /////////////////////////////////
 const char* names[7] {"Pi+", "Pi-", "K+", "K-", "p+", "p-", "K*"};
 for (int i = 0; i < 7; ++i)
 cout << names[i] << " generated: "<< hfTypes->GetBinContent(i + 1) << " +/- " << hfTypes->GetBinError(i + 1) << endl; 
 
 /////////////////////////////////
 ////Histogram Subtractions///////
 ///////////////////////////////// 
 
 TH1F * hK3_4 = new TH1F("hK3_4", "K* 3-4", 160, 0, 4);
 hK3_4->Sumw2();
 hK3_4->Add(h3, h4, 1, -1);
 TH1F * hK1_2 = new TH1F("hK1_2", "K* 1-2", 160, 0, 4);
 hK1_2->Sumw2();
 hK1_2->Add(h1, h2, 1, -1);
 
 hK1_2->SetEntries(hK1_2->Integral());
 hK3_4->SetEntries(hK3_4->Integral());
 
 /////////////////////////////////
 ////Fitting//////////////////////
 ///////////////////////////////// 
 hfPhi->Fit("pol0", "Q");
 hfTheta->Fit("pol0", "Q");
 hfP->Fit("expo", "Q");
 
 hK1_2-> Fit("gaus","RSQ", "", 0.6, 1.2);
 hK3_4-> Fit("gaus","RSQ", "", 0.6, 1.2);
 h5-> Fit("gaus", "Q");
 

 /////////////////////////////////
 ////Storage//////////////////////
 ///////////////////////////////// 
 TFile * file2 = new TFile("fit.root", "RECREATE");
 
 hfPhi->Write();
 hfTheta->Write();
 hfP->Write();
 h5->Write();
 hK3_4->Write();
 hK1_2->Write();
 
 file2->Close();
 
 
 /////////////////////////////////
 ////Display//////////////////////
 /////////////////////////////////  
 
 //First Canvas
 TCanvas* can1 = new TCanvas("can1", "", 200, 100, 1000, 600);
  can1->Divide(2,2);
  
  can1->cd(1);
  hfTypes->SetTitle("Particle Abundances");
  hfTypes->GetXaxis()->SetTitle("ParticleType");
  hfTypes->GetYaxis()->SetTitle("Abundance");
  gPad->SetGrid();
  gPad->SetFillColor(42); //17
  gPad->SetFrameFillColor(20); //19
  hfTypes->SetFillColor(36);
  hfTypes->DrawCopy();
  
  can1->cd(2);
  hfP->SetTitle("Impulse Distribution");
  hfP->GetXaxis()->SetTitle("Impulse");
  hfP->GetYaxis()->SetTitle("Occurrences");
  hfP->SetFillColor(46);
  gPad->SetGrid();
  gPad->SetFillColor(42);
  gPad->SetFrameFillColor(20);
  hfP->DrawCopy();
  
  can1->cd(3);
  hfPhi->SetTitle("Polar Angle Distribution");
  hfPhi->GetXaxis()->SetTitle("Angle(rad)");
  hfPhi->GetYaxis()->SetTitle("Occurrences");
  gPad->SetGrid();
  gPad->SetFillColor(42);
  gPad->SetFrameFillColor(20);
  hfPhi->SetFillColor(34);
  hfPhi->DrawCopy();  
  
  can1->cd(4);
  hfTheta->SetTitle("Azimutal Angle Distribution");
  hfTheta->GetXaxis()->SetTitle("Angle(rad)");
  hfTheta->GetYaxis()->SetTitle("Occurrences");
  gPad->SetGrid();
  gPad->SetFillColor(42);
  gPad->SetFrameFillColor(20);
  hfTheta->SetFillColor(34);
  hfTheta->DrawCopy();    
  
  //Second Canvas
  TCanvas* can2 = new TCanvas("can2", "", 200, 100, 1000, 600);
  can2->Divide(2,2);
  
  can2->cd(1);
  hK1_2->GetXaxis()->SetRangeUser(0, 2);
  hK1_2->GetXaxis()->SetTitle("Invariant Mass");
  hK1_2->GetYaxis()->SetTitle("Occurrences");
  hK1_2->SetTitle("Subtraction 1");
  gPad->SetGrid();
  gPad->SetFillColor(42);
  gPad->SetFrameFillColor(19);
  hK1_2->SetFillColor(30);
  hK1_2->DrawCopy("hist");
  hK1_2->GetFunction("gaus")->DrawCopy("same");
  
  can2->cd(2);
  hK3_4->GetXaxis()->SetRangeUser(0, 2);
  hK3_4->GetXaxis()->SetTitle("Invariant Mass");
  hK3_4->GetYaxis()->SetTitle("Occurrences");
  hK3_4->SetTitle("Subtraction 2");
  gPad->SetGrid();
  gPad->SetFillColor(42);
  gPad->SetFrameFillColor(19);
  hK3_4->SetFillColor(30);
  hK3_4->DrawCopy("hist");
  hK3_4->GetFunction("gaus")->DrawCopy("same");

  
  can2->cd(3);
  h5->GetXaxis()->SetRangeUser(0, 2);
  h5->GetXaxis()->SetTitle("Invariant Mass");
  h5->GetYaxis()->SetTitle("Occurrences");
  h5->SetTitle("Theoretical K*");
  gPad->SetGrid();
  gPad->SetFillColor(42);
  gPad->SetFrameFillColor(19);
  h5->SetFillColor(30);
  h5->DrawCopy();
  
  can2->cd(4);
  gPad->SetFillColor(42);

  can1->Print("FirstCanvas.gif");
  can1->Print("FirstCanvas.C");
  can1->Print("FirstCanvas.root");  
  
  can2->Print("SecondCanvas.gif");
  can2->Print("SecondCanvas.C");
  can2->Print("SecondCanvas.root");
 
}
