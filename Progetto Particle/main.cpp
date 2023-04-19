#include "Particle.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"


int main(){
gRandom->SetSeed();
TH1D* histTypes = new TH1D("HistTypes","Particles Types Generated", 7, 0, 7);
TH2D* histAngles = new TH2D("HistAngles", "Distribution  Angle", 100, 0, 2*TMath::Pi(), 50, 0, TMath::Pi());
TH1D* histP = new TH1D("HistP", "Impulse", 500, 0, 5);
TH1D* histPt = new TH1D("HistPt", "Transverse Impulse", 500, 0, 5);
TH1D* histEnergy = new TH1D("HistEnergy", "Energy", 500, 0, 5);
TH1D* histIMall = new TH1D("HistIMall", "Invariant Mass between every particle", 500, 0, 4);
TH1D* histIMsc = new TH1D("HistIMsc", "Invariant Mass between every same charged particle", 500, 0, 4);
TH1D* histIMoc= new TH1D("HistIMoc", "Invariant Mass between every opposite charged particle", 500, 0, 4);
TH1D* histIMKPoc = new TH1D("HistIMKPoc", "Invariant Mass between K+ Pi- ", 500, 0, 4);
TH1D* histIMKPsc = new TH1D("HistIMKPsc", "Invariant Mass between K+ Pi+", 500, 0, 4);
TH1D* histIMDecay = new TH1D("HistIMDecay", "Invariant Mass between Products of Decay", 500, 0, 4);


Particle::AddParticleType("Pi+", 0.13957, 1);		//index = 0  Pi+
Particle::AddParticleType("Pi-", 0.13957, -1);			//1  Pi-
Particle::AddParticleType("K+", 0.49367, 1);			//2  K+
Particle::AddParticleType("K-", 0.49367, -1);			//3  K-
Particle::AddParticleType("p+", 0.93827, 1);			//4  p+
Particle::AddParticleType("p-", 0.93827, -1);			//5  p-
Particle::AddParticleType("K*", 0.89166, 0, 0.050);		//6  K*


int nEvents = 100000;
int nPartForEvent = 110;

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
        count+=2;
    }
    
    histTypes->Fill(particle[i].GetIndex());
    histAngles->Fill(phi, theta);
    histP->Fill(P);
    histPt->Fill(sqrt(pow(Px,2)+pow(Py,2)));
    histEnergy->Fill(particle[i].Energy());
    
  }
  for (int i = 0; i < 100 + count - 1; ++i){
    for(int j = i + 1; j < 100 + count - 1; ++j){
      histIMall->Fill(particle[i].InvMass(particle[j]));
      
      if (particle[i].GetCharge() * particle[j].GetCharge() == 1)
        histIMsc->Fill(particle[i].InvMass(particle[j]));
      else if(particle[i].GetCharge() * particle[j].GetCharge() == -1)
        histIMoc->Fill(particle[i].InvMass(particle[j]));
      
      if((particle[i].GetIndex() == 0 &&  particle[j].GetIndex() == 3) || 
      	 (particle[i].GetIndex() == 3 &&  particle[j].GetIndex() == 0) ||
      	 (particle[i].GetIndex() == 1 &&  particle[j].GetIndex() == 2) ||        
      	 (particle[i].GetIndex() == 2 &&  particle[j].GetIndex() == 1))
      	histIMKPoc->Fill(particle[i].InvMass(particle[j]));
      	
      else if((particle[i].GetIndex() == 0 &&  particle[j].GetIndex() == 2) || 
      	 (particle[i].GetIndex() == 2 &&  particle[j].GetIndex() == 0) ||
      	 (particle[i].GetIndex() == 1 &&  particle[j].GetIndex() == 3) ||        
      	 (particle[i].GetIndex() == 3 &&  particle[j].GetIndex() == 1))
      	histIMKPsc->Fill(particle[i].InvMass(particle[j]));	
      	
     
     
    }  	     
  }
  
  for (int m  = 0; m < count; ++m){
    for(int n = m + 1; n < count; ++n){
      histIMDecay->Fill(particle[100 + m].InvMass(particle[100 + n]));
  }
}
}


TFile* lab = new TFile("lab.root", "RECREATE");
 histTypes->Write();
 histAngles->Write();
 histP->Write();
 histPt->Write();
 histEnergy->Write();
 histIMall->Write();
 histIMsc->Write();
 histIMoc->Write();
 histIMKPoc->Write();
 histIMKPsc->Write();
 histIMDecay->Write();
lab->Close();
 
}


