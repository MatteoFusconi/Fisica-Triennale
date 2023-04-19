#include "Particle.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"


void gen(){
gRandom->SetSeed();
TH1F* histTypes = new TH1F("HistTypes","Particles Types Generated", 7, 0, 7);
TH1F* histPhi = new TH1F("HistPhi", "Distribution  Phi", 100, 0, 2*TMath::Pi());
TH1F* histTheta = new TH1F("HistTheta", "Distribution  Theta", 50, 0, TMath::Pi());
TH1F* histP = new TH1F("HistP", "Impulse", 500, 0, 5);
TH1F* histPt = new TH1F("HistPt", "Transverse Impulse", 500, 0, 5);
TH1F* histEnergy = new TH1F("HistEnergy", "Energy", 500, 0, 5);
TH1F* histIMall = new TH1F("HistIMall", "Invariant Mass between every particle", 160, 0, 4);
TH1F* histIM1 = new TH1F("HistIM1", "Invariant Mass same charge", 160, 0, 4);
TH1F* histIM2 = new TH1F("HistIM2", "Invariant Mass opposite charge", 160, 0, 4);
TH1F* histIM3 = new TH1F("HistIM3", "Invariant Mass K Pi opposite charge ", 160, 0, 4);
TH1F* histIM4 = new TH1F("HistIM4", "Invariant Mass K Pi same charge", 160, 0, 4);
TH1F* histIMDecay = new TH1F("HistIMDecay", "Invariant Mass Decay", 160, 0, 4);


Particle::AddParticleType("Pi+", 0.13957, 1);		//index = 0  Pi+
Particle::AddParticleType("Pi-", 0.13957, -1);			//1  Pi-
Particle::AddParticleType("K+", 0.49367, 1);			//2  K+
Particle::AddParticleType("K-", 0.49367, -1);			//3  K-
Particle::AddParticleType("p+", 0.93827, 1);			//4  p+
Particle::AddParticleType("p-", 0.93827, -1);			//5  p-
Particle::AddParticleType("K*", 0.89166, 0, 0.050);		//6  K*


int nEvents = 100000;
int nPartForEvent = 120; 

Particle particle[nPartForEvent];

for (int ev = 0; ev < nEvents; ++ev){
  int count = 0;
  for (int i = 0; i < 100; ++i){
    double phi, theta, P, Px, Py, Pz;
    phi = 2 * TMath::Pi() * gRandom->Uniform();
    theta = TMath::Pi() * gRandom->Uniform();
    P = gRandom->Exp(1); // 1GeV
    Px = P * TMath::Sin(theta) * TMath::Cos(phi);
    Py = P * TMath::Sin(theta) * TMath::Sin(phi);
    Pz = P * TMath::Cos(theta);
    particle[i].SetP(Px, Py, Pz);
    
    double gen = gRandom->Uniform();
    if (gen < .4)		particle[i].SetIndex(0);
    else if (gen < .8)		particle[i].SetIndex(1);
    else if (gen < .85)	particle[i].SetIndex(2);
    else if (gen < .9)		particle[i].SetIndex(3);
    else if (gen < .945)	particle[i].SetIndex(4);
    else if (gen < .99)	particle[i].SetIndex(5);
    else{			particle[i].SetIndex(6);
      if(gen > .995){
	particle[100+count].SetIndex(0);
	particle[101+count].SetIndex(3);
        }
      else{
        particle[100+count].SetIndex(1);
	particle[101+count].SetIndex(2);
        }
      particle[i].Decay2body(particle[100+count], particle[101+count]);
      histIMDecay->Fill(particle[100 + count].InvMass(particle[101 + count]));
        count+=2;
    }
    
    histTypes->Fill(particle[i].GetIndex());
    histPhi->Fill(phi);
    histTheta->Fill(theta);
    histP->Fill(P);
    histPt->Fill(sqrt(pow(Px,2)+pow(Py,2)));
    histEnergy->Fill(particle[i].Energy());
    
  }
  for (int i = 0; i < 100 + count - 1; ++i){
    for(int j = i + 1; j < 100 + count; ++j){
      histIMall->Fill(particle[i].InvMass(particle[j]));
      
      if (particle[i].GetCharge() * particle[j].GetCharge() == 1)
        histIM1->Fill(particle[i].InvMass(particle[j]));
      else if (particle[i].GetCharge() * particle[j].GetCharge() == -1)
        histIM2->Fill(particle[i].InvMass(particle[j]));
      
      if((particle[i].GetIndex() == 0 &&  particle[j].GetIndex() == 3) || 
      	 (particle[i].GetIndex() == 3 &&  particle[j].GetIndex() == 0) ||
      	 (particle[i].GetIndex() == 1 &&  particle[j].GetIndex() == 2) ||        
      	 (particle[i].GetIndex() == 2 &&  particle[j].GetIndex() == 1))
      	histIM3->Fill(particle[i].InvMass(particle[j]));
      	
      else if((particle[i].GetIndex() == 0 &&  particle[j].GetIndex() == 2) || 
      	 (particle[i].GetIndex() == 2 &&  particle[j].GetIndex() == 0) ||
      	 (particle[i].GetIndex() == 1 &&  particle[j].GetIndex() == 3) ||        
      	 (particle[i].GetIndex() == 3 &&  particle[j].GetIndex() == 1))
      	histIM4->Fill(particle[i].InvMass(particle[j]));	     	
    }
  }
  
}


TFile* lab2 = new TFile("lab2.root", "RECREATE");
 histTypes->Write();
 histPhi->Write();
 histTheta->Write();
 histP->Write();
 histPt->Write();
 histEnergy->Write();
 histIMall->Write();
 histIM1->Write();
 histIM2->Write();
 histIM3->Write();
 histIM4->Write();
 histIMDecay->Write();
lab2->Close();
 
}
