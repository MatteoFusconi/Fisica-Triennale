#ifndef PARTICLETYPES_H
#define	PARTICLETYPES_H

class ParticleType {
public:
	ParticleType(const char* name, double mass, int ch);

	const char* GetName() const;
	double GetMass() const;
	int GetCharge() const;
	virtual double GetWidth() const;

	void Print() const;

private:
	const char* fName;
	const double fMass;
	const int fCharge;
};

class ResonanceType : public ParticleType {
public:
	ResonanceType(const char* name, double mass, int ch, double wi);

	double GetWidth() const;

	void Print() const;

private:
	const double fWidth;
};


#endif