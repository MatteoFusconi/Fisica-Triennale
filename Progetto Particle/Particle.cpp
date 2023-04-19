#include "Particle.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // !M_PI



ParticleType* Particle::fParticleType[fMaxNumParticleType];
int Particle::fNParticleType = 0;

Particle::Particle(const char* name, double Px = 0., double Py = 0., double Pz = 0.): fPx(Px), fPy(Py), fPz(Pz) {
	fIndex = FindParticle(name);
	if (fIndex == -1)
		std::cout << "Error! There is no such Particle named " << name << std::endl;
}
Particle::Particle(): fIndex (-1), fPx(0), fPy(0), fPz(0){}

int Particle::FindParticle(const char* parname) {
	for (int i = 0; i < fNParticleType; ++i) {
		if (fParticleType[i]->GetName() == parname)
			return i;
	}
	return -1;
}

const char* Particle::GetName() const {
	return fParticleType[fIndex]->GetName();
}
double Particle::GetMass() const {
	return fParticleType[fIndex]->GetMass();
}
int Particle::GetCharge() const {
	return fParticleType[fIndex]->GetCharge();
}

int Particle::GetIndex() const {
	return fIndex;
}

double Particle::GetPx() const {
	return fPx;
}
double Particle::GetPy() const {
	return fPy;
}
double Particle::GetPz() const {
	return fPz;
}

void Particle::AddParticleType(const char* aname, double amass, int ach, double awi) {
	if (FindParticle(aname) == -1) {
		if (awi == 0)
			fParticleType[fNParticleType] = new ParticleType(aname, amass, ach);
		else
			fParticleType[fNParticleType] = new ResonanceType(aname, amass, ach, awi);
		++fNParticleType;
	}
	else
		std::cout << "Particle named " << aname << " already exists" << std::endl;
}

void Particle::SetIndex(int i) {
	fIndex = i;
}
void Particle::SetIndex(const char* name) {
	fIndex = FindParticle(name);
}
void Particle::SetP(double px, double py, double pz) {
	fPx = px;
	fPy = py;
	fPz = pz;
}

void Particle::PrintParticleTypes() {
	for (int i = 0; i < fNParticleType; ++i) {	
		std::cout << "Particle " << i << ":\n";
		fParticleType[i]->Print();
		std::cout << "\n\n";
	}
}

void Particle::Print() const {
	std::cout << "Particle " << fIndex << ":\n";
	fParticleType[fIndex]->Print();
	std::cout << "Momentum = (" << fPx << ", " << fPy << ", " << fPz << ")\n";
}

double Particle::Energy() const {
	return sqrt(pow(fParticleType[fIndex]->GetMass(), 2) +
		pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2));
}

double Particle::InvMass(Particle& p) const {
	return sqrt(pow(Energy() + p.Energy(), 2) - pow(fPx + p.GetPx(), 2) -
		pow(fPy + p.GetPy(), 2) - pow(fPz + p.GetPz(), 2));
}

int Particle::Decay2body(Particle& dau1, Particle& dau2) const {
	if (GetMass() == 0.0) {
		printf("Decayment cannot be preformed if mass is zero\n");
		return 1;
	}

	double massMot = GetMass();
	double massDau1 = dau1.GetMass();
	double massDau2 = dau2.GetMass();

	if (fIndex > -1) { // add width effect

	  // gaussian random numbers

		float x1, x2, w, y1, y2;

		double invnum = 1. / RAND_MAX;
		do {
			x1 = 2.0 * rand() * invnum - 1.0;
			x2 = 2.0 * rand() * invnum - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w)) / w);
		y1 = x1 * w;
		y2 = x2 * w;

		massMot += fParticleType[fIndex]->GetWidth() * y1;

	}

	if (massMot < massDau1 + massDau2) {
		printf("Decayment cannot be preformed because mass is too low in this channel\n");
		return 2;
	}

	double pout = sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

	double norm = 2 * M_PI / RAND_MAX;

	double phi = rand() * norm;
	double theta = rand() * norm * 0.5 - M_PI / 2.;
	dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));
	dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));

	double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

	double bx = fPx / energy;
	double by = fPy / energy;
	double bz = fPz / energy;

	dau1.Boost(bx, by, bz);
	dau2.Boost(bx, by, bz);

	return 0;
}

void Particle::Boost(double bx, double by, double bz)
{

	double energy = Energy();

	//Boost this Lorentz vector
	double b2 = bx * bx + by * by + bz * bz;
	double gamma = 1.0 / sqrt(1.0 - b2);
	double bp = bx * fPx + by * fPy + bz * fPz;
	double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

	fPx += gamma2 * bp * bx + gamma * bx * energy;
	fPy += gamma2 * bp * by + gamma * by * energy;
	fPz += gamma2 * bp * bz + gamma * bz * energy;
}