// asteroid.h
// Root of the classes' hierarchy of asteroid shape models.
// Purely abstract.

// Author: michael.mueller@dlr.de
// Date:   Nov. 18 2003

// body file: asteroid.cpp

// REMARK:
// Note that spinAxis is SpinState rather than reference to one ->
// it is save to assign it to a temporary object! (it's copied, anyway!)



#if !defined(AFX_asteroid_H__58D2ED80_92B4_47D9_99C2_804BB4447285__INCLUDED_)
#define AFX_asteroid_H__58D2ED80_92B4_47D9_99C2_804BB4447285__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<stdexcept>
#include<new>
#include<map>
#include "SpinState.h"
#include "constants.h"


class asteroid  
{
public:
	class fitFileMJy;	// defined after asteroid is defined
	class fitFileSI;    // fluxes in W/m^2/mu, as opposed to mJy (fitFileMJy) 
	class fluxRecord;   // defined at end of asteroid

	class ThermalModel;		// abstract base classes for thermal + optical modelling
	class ScatteringModel; 
protected:
	// protected ctor -- base-class!!!
	asteroid(double pV, 
		  const SpinState& spinAxis, // length is angular velocity in h^{-1}
		  double H, double G=0.15, double emissivity=0.9)
		  throw(std::invalid_argument);
public:
	// converts flux (in multiples of solar flux!!!) into magnitudes
	static inline double flux2mag(double flux // in units of solar flux
		) throw (std::invalid_argument);
	// converts a given value of thermal inertia (J/m^2/K/s^{1/2}) into the thermal parameter (Lebofsky/Spencer)
	double calculateThermalParameter (double ThermalInertiaSI, double rAU) const;
	// output in J/m^2/K/sqrt(s)
	double calculateThermalInertia (double ThermalParameter, double rAU) const;

protected:
// data members:	
	// effectice diameter (km)
	double diameter;
	inline void setDiameter()
		{diameter=1329.*pow(10., -H/5.)/sqrt(pV);};
	// geometric albedo
	double pV;
	// bolometric Bond albedo
	double bond;
	void setBond()
		{bond = AsteroidConstants::q(G) * pV;};
	double emissivity;

	// length is omega / h^{-1}
	SpinState spinAxis;
	
	double H, G;
	double period;

public:  
// methods to change / read out the data members
	inline double getDiameter() const throw()
		{return diameter;};
	virtual inline void setPV (double newAlbedo)
		throw(std::invalid_argument);
	double getPV() const throw()
		{return pV;};
	double getBond() const throw()
		{return bond;};
	inline void setEmissivity(double emissivity)
		throw(std::invalid_argument);
	double getEmissivity() const throw()
		{return emissivity;};
	virtual inline void setH(double H);
	double getH() const throw()
		{return H;};
	inline void setG(double G);
	double getG() const throw()
		{return G;};	
protected:
	// spins vector vec along z-axis 'forward in time', i.e. clockwise (t/h)
	inline void spin (vector3& vec, double tH) const throw();
	// spins the vector vec, returns a new one (keeping vec const)
	inline vector3 spin (const vector3& vec, double tH) const throw();

	// transforms unit-vector given in ecliptic coordinates into bodyfixed system
	// (potential length is ignored -> output is of unit length!!!!)
	// Parm JD0=0 -> Leave at zero time!
	vector3 bodyFix (double ecl_latitude, double ecl_longitude, double JD=0) const throw();
	// transforms vector given in ecliptic coordinates into bodyfixed system
	// (potential length is ignored -> output is of unit length!!!!)
	vector3 bodyFix (const eclipticVector& vec, double JD=0) const throw();
public:
	virtual ~asteroid() throw() {};
private:
	// prohibitive
	asteroid();
	asteroid(const asteroid& );
	asteroid& operator= (const asteroid&);
}; // class asteroid



class asteroid::ThermalModel 
{ // abstract basis class for thermal modelling
public:
	virtual ~ThermalModel()throw(){};
protected:
	ThermalModel (double emissivity, double bond);
	ThermalModel (const asteroid& shape);
	
	const double emissivity, bond;
	
	// returns subsolar temperature depending on heliocentric distance (in AU!)
	inline double TSS (double rAU // heliocentric distance in AU
		, double cosSun           // cos(Sun's height over horizon)
		, double eta=1            // beaming parameter
		) const
	{return sqrt(sqrt( (1-bond)*cosSun*AsteroidConstants::solar_over_stefan_boltzmann/(emissivity*eta)) /rAU);};
	// looks up
	inline double TSS (double rAU) const;
	
	mutable double accuracyGoal;
	mutable std::map<double, double> TSS_lookup;
	static const double DefaultAccuracyGoal;
public:
	virtual void setAccuracy(double accuracyGoal) const;
	double getAccuracy() const
	{return accuracyGoal;};	
	double getEmissivity() const throw()
	{return emissivity;};
	double getBond() const throw()
	{return bond;};
	virtual void setBond(double bond) throw(std::invalid_argument)
	{
		if (this->bond == bond)
			return;
		if (bond<0 || bond>1)
			throw std::invalid_argument("asteroid::ThermalModel::setBond: Incompatible Bond albedo entered (<0 or >1)!");
		const_cast<double&>(this->bond) = bond;
		TSS_lookup.clear();
	};
	virtual void setEmissivity(double emissivity) throw (std::invalid_argument)
	{
		if (this->emissivity == emissivity)
			return;
		if (emissivity<0 || emissivity>1)
			throw std::invalid_argument("asteroid::ThermalModel::setEmissivity: Incompatibel emissivity entered (<0 or >1)!");
		const_cast<double&>(this->emissivity) = emissivity;
		TSS_lookup.clear();
	};
	double getTSS(double rAU) const
		{return TSS(rAU);};
private:
	inline void check();
	ThermalModel();
	ThermalModel (const ThermalModel&);
	ThermalModel& operator=(const ThermalModel&);
}; // class asteroid::ThermalModel



class ScatteringModel
{ // abstract basis class for the description of optical scattering models
	public:
		virtual ~ScatteringModel()=0;
	protected:
		ScatteringModel(){};
	private:
		ScatteringModel (const ScatteringModel&);
		ScatteringModel& operator= (const ScatteringModel&);
}; // class asteroid::ScatteringModel



class asteroid::fluxRecord
// Don't use this as an end-user!
// Use classes fitFileMJy or fitFileSI!
{
public:
	// reads in a list of flux records from inName
	// returns #fluxes, records are stored in that array (must be delete[]d at the end!)
	// file syntax: n = number of data lines,
	//	then: n lines with the format:
	//			1.	ephemeris as given by horizons using: 
	//				JD, heliocentric ecl. coordinates, r, delta, obs.centric ecl. coordinates,
	//			2. lambdaMu  flux/(W/m^2/mu) sigma/(W/m^2/mu)
	// IMPORTANT: The CSV-option must be turned ON in horizons, the last comma must be copied into the input file,
	//		after the ephemeris follows the actual data WITHOUT commas
	static const unsigned int readSI (fluxRecord**& records, const char inName[]);
	
	double JD, lambdaMu, fluxSI, sigmaSI;
	eclipticVector Sun2Aster, Earth2Aster;
	
	fluxRecord(double JD, double lambdaMu, double fluxSI, double sigmaSI,
		double heliocentricLongitude, double heliocentricLatitude, double rAU,
		double geocentricLongitude, double geocentricLatitude, double deltaAU)
		:   JD(JD), lambdaMu(lambdaMu), fluxSI(fluxSI), sigmaSI(sigmaSI),
		Sun2Aster(heliocentricLongitude, heliocentricLatitude, rAU),
		Earth2Aster(geocentricLongitude, geocentricLatitude, deltaAU)
	{};
	
	fluxRecord(double JD, double lambdaMu, double fluxSI, double sigmaSI,
		const eclipticVector& Sun2Aster, 
		const eclipticVector& Earth2Aster)
		: JD(JD), lambdaMu(lambdaMu), fluxSI(fluxSI), sigmaSI(sigmaSI),
		Sun2Aster(Sun2Aster), Earth2Aster(Earth2Aster)
	{};
	
	virtual ~fluxRecord(){};
}; // class asteroid::fluxRecord



class asteroid::fitFileSI
{
public:
	fitFileSI (const char inName[])
		: nRecords(fluxRecord::readSI(records, inName))
	{};
	virtual ~fitFileSI();
	
	const fluxRecord& operator[] (unsigned int i) const
	{
		if (i>=nRecords) throw std::invalid_argument("");
		return *records[i];
	};
	const unsigned int nRecords;
protected:
	fluxRecord** records;
}; // class asteroid::fitFileSI



class asteroid::fitFileMJy : public asteroid::fitFileSI
{ // reads in file as if fluxes were in W/m^2/mu - converts them in ctor.
public:
	fitFileMJy(const char inName[]);
	virtual ~fitFileMJy(){};
}; //class asteroid::fitFileMJy






// !!!!!!!!!!!!!!!!!!!!!!!!!
// Inline implementations:
// !!!!!!!!!!!!!!!!!!!!!!!!!







inline void asteroid::setPV(double newAlbedo)
	throw(std::invalid_argument)
{
	if (newAlbedo<=0)
		throw std::invalid_argument("Albedo must be nonnegative!");
	if (newAlbedo>3)
		throw std::invalid_argument
		  ("A geometric albedo higher than 3 hardly makes sense. \nIf you disagree, please consult the developer!");
	pV=newAlbedo;
	bond=AsteroidConstants::q(G)*pV;
	setDiameter();
}; //asteroid::setAlbedo


inline void asteroid::setH(double H)
{
	this->H=H;
	setDiameter();
};
	

inline void asteroid::setG(double G)
{
	this->G = G;
	setBond();
};


inline vector3 asteroid::spin(const vector3 &vec, double tH) const throw()
{
	vector3 dummy=vec;
	spin(dummy, tH);
	return dummy;	
} // vector3 asteroid::spin(const vector3 &vec, double t)



inline void asteroid::spin(vector3 &vec, double tH) const throw()
{
	vec.rzcw(spinAxis.getLength()*tH); // rotate clockwise!!!!!!! 
				// (Asteroid rotates ccw -> in asterocentric system fixed vector3s rotated cw)
} // void asteroid::spin(vector3 &vec, double t)


inline double asteroid::flux2mag(double flux) throw(std::invalid_argument)
{
	if (flux==0.) throw std::invalid_argument("flux=0 in flux2mag encountered");
	return -26.74 - 2.5*log10(flux);
};


inline void asteroid::setEmissivity(double emissivity) throw(std::invalid_argument)
{
	if (emissivity<=0 || emissivity > 1)
		throw std::invalid_argument("Error in asteroid::setEmissivity: 0<emissivity<=1 expected!");
	this->emissivity = emissivity;
};


inline double asteroid::ThermalModel::TSS (double rAU) const
{
	std::map<double,double>::iterator it = TSS_lookup.find(rAU);
	if (it != TSS_lookup.end())
		return it->second; // rAU already in lookup-list
	else
	{
		double result = sqrt(sqrt( (1-bond)*AsteroidConstants::solar_over_stefan_boltzmann/(emissivity)) /rAU);
		TSS_lookup.insert(std::pair<double,double>(rAU, result));
		return result;
	};
};




#endif // !defined(AFX_asteroid_H__58D2ED80_92B4_47D9_99C2_804BB4447285__INCLUDED_)
