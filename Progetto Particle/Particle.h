#ifndef PARTICLE_H
#define PARTICLE_H

#include "ParticleTypes.h"

class Particle {
public:
	Particle(const char* name, double Px, double Py, double Pz);
	Particle();

	const char* GetName() const;
	double GetMass() const;
	int GetCharge() const;
	int GetIndex() const;
	double GetPx() const;
	double GetPy() const;
	double GetPz() const;

	static void AddParticleType(const char* aname, double amass, int ach, double awi = 0);
	
	void SetIndex(int i);
	void SetIndex(const char* name);
	void SetP(double px, double py, double pz);
	
	static void PrintParticleTypes();
	void Print() const;
	double Energy() const;
	double InvMass(Particle& p) const;

	int Decay2body(Particle& dau1, Particle& dau2) const;

private:
	static const int fMaxNumParticleType = 10;
	static ParticleType* fParticleType[fMaxNumParticleType];
	static int fNParticleType;

	int fIndex;
	double fPx, fPy, fPz;

	static int FindParticle(const char* parname);

	void Boost(double bx, double by, double bz);
};






#endif