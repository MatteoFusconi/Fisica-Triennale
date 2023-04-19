#include "ParticleTypes.h"
#include <iostream>

ParticleType::ParticleType(const char* name, double mass, int ch): fName(name), fMass(mass), fCharge(ch) {}

const char* ParticleType::GetName() const {
	return fName;
}

double ParticleType::GetMass() const {
	return fMass;
}

int ParticleType::GetCharge() const {
	return fCharge;
}

double ParticleType::GetWidth() const {
	return 0;
}

void ParticleType::Print() const {
	std::cout << "Name: \t" << fName << "\n\tMass = " 
		<< fMass << "\n\tCharge = " << fCharge << std::endl;
}

//Resonance Type

ResonanceType::ResonanceType(const char* name, double mass, int ch, double wi) : ParticleType(name, mass, ch), fWidth(wi){}

double ResonanceType::GetWidth() const {
	return fWidth;
}

void ResonanceType::Print() const {
	ParticleType::Print();
	std::cout << "\tWidth: " << fWidth << std::endl;
}