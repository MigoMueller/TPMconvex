// shape.h
// Root of the classes' hierarchy of asteroid shape models.
// Purely abstract.

// Author: michael.mueller@dlr.de
// Date:   Nov. 18 2003

// body file: shape.cpp

// REMARK:
// Note that spinAxis is SpinState rather than reference to one ->
// it is save to assign it to a temporary object! (it's copied, anyway!)


// BAUSTELLE: The optional parameter phi0 is not handled consistently
// in fact: not at all...


#if !defined(AFX_SHAPE_H__58D2ED80_92B4_47D9_99C2_804BB4447285__INCLUDED_)
#define AFX_SHAPE_H__58D2ED80_92B4_47D9_99C2_804BB4447285__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<stdexcept>
#include<new>
#include "SpinState.h"
#include "../constants.h"


class shape  
{
public:
	class fitFile;		// defined at end of shape
	class fluxRecord;   // defined at end of shape
protected:
	// effectice diameter (km)
	double diameter;
	double pV;
	double bond;
	// length is omega / h^{-1}
	const SpinState spinAxis;
	const double H;
	const double G;
	const double period;
	double emissivity;
public:
	shape(double pV, 
		  const SpinState& spinAxis, // length is angular velocity in h^{-1}
		  double H, double G=0.15)
		  throw(std::invalid_argument);

	// sets the internal variable pV
	virtual inline void setPV (double newAlbedo)
		throw(std::invalid_argument);
	double getPV() const
		{return pV;};
	double getBond() const
		{return bond;};
	inline void setEmissivity(double emissivity)
		throw(std::invalid_argument);
	double getEmissivity() const throw()
		{return emissivity;};

	// returns subsolar temperature depending on heliocentric distance (in AU!)
	inline double TSS (double rAU // heliocentric distance in AU
		, double cosSun=1         // cos(Sun's height over horizon)
		, double eta=1            // beaming parameter
		) const
		{return sqrt(sqrt( (1-getBond())*cosSun*AsteroidConstants::solar_over_stefan_boltzmann/getEmissivity()/eta) /rAU);};

	void plotLambertianLC(const char* const filename,
								  const eclipticVector& Sun2Aster,
								  const eclipticVector& Earth2Aster,
								  const unsigned int plotPoints = 40,
								  unsigned int periods=1)
								  const 
								  throw (std::invalid_argument)
		{plotMikkoLC(filename, 0, 1, Sun2Aster, Earth2Aster, plotPoints, periods);};
	void plotLommelSeeligerLC(    const char* const filename,
								  const eclipticVector& Sun2Aster,
								  const eclipticVector& Earth2Aster,
								  const unsigned int plotPoints = 40,
								  unsigned int periods=1)
								  const
								  throw (std::invalid_argument)
		{plotMikkoLC(filename, 1, 0, Sun2Aster, Earth2Aster, plotPoints, periods);};
	virtual void plotGeometricLC (const char* const filename,
								  const eclipticVector& Sun2Aster,
								  const eclipticVector& Earth2Aster,
								  const unsigned int plotPoints = 40,
								  unsigned int periods=1)
								  const
								  throw (std::invalid_argument)
								=0;
	virtual void plotMikkoLC     (const char* const filename,
								  double LS, double Lambert,
								  const eclipticVector& Sun2Aster,
								  const eclipticVector& Earth2Aster,
								  const unsigned int plotPoints = 40,
								  unsigned int periods=1)
								  const
								  throw (std::invalid_argument)
								=0;

	// outputs a file to be plotted by the freeware-raytracer POVray
	// The result will show the temperature distribution on the surface according to NEATM
	// Light source is the Sun, the observer is on Earth.
	virtual void plotNEATM_POV (const char* const filename, 
								const eclipticVector& Sun2Aster, 
								const eclipticVector& Earth2Aster, 
								double eta, double rotationalPhase,
								int nColors) 
								const
								throw (std::invalid_argument, std::bad_alloc)
							=0;

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

	static inline double flux2mag(double flux) throw (std::invalid_argument);

	inline void setDiameter()
		{diameter=1329.*pow(10., -H/5.)/sqrt(pV);};
public:
	virtual ~shape() throw() {};
private:
	// prohibitive
	shape();
	shape(const shape& );
	shape& operator= (const shape&);

public:
	class fluxRecord
	{
	public:
		// reads in a list of flux records from inName
		// returns #fluxes, records are stored in that array (must be delete[]d at the end!)
		// file syntax: n = number of data lines,
		//	then: n lines with the format:
		//			1.	ephemeris as given by horizons using: 
		//				JD, heliocentric ecl. coordinates, r, delta, obs.centric ecl. coordinates,
		//			2. lambdaMu  flux/mJy sigma/mJy
		// IMPORTANT: The CSV-option must be turned ON in horizons, the last comma must be copied into the input file,
		//		after the ephemeris follows the actual data WITHOUT commas
		static const unsigned int read(fluxRecord**& records, const char inName[]);
		
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
	}; // class shape::fluxRecord
		
		
	class fitFile
	{
	public:
		fitFile(const char inName[]) 
			:
		nRecords(fluxRecord::read(records, inName))
		{};
		virtual ~fitFile();
		
		const fluxRecord& operator[] (unsigned int i) const
		{
			if (i>=nRecords)
				throw std::invalid_argument("");
			return *records[i];
		};
		
		const unsigned int nRecords;
	protected:
		fluxRecord** records;
	}; //class shape::fitFile

	
}; // class shape


// !!!!!!!!!!!!!!!!!!!!!!!!!
// Inline implementations:
// !!!!!!!!!!!!!!!!!!!!!!!!!


inline void shape::setPV(double newAlbedo)
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
}; //shape::setAlbedo



inline vector3 shape::spin(const vector3 &vec, double tH) const throw()
{
	vector3 dummy=vec;
	spin(dummy, tH);
	return dummy;	
} // vector3 shape::spin(const vector3 &vec, double t)



inline void shape::spin(vector3 &vec, double tH) const throw()
{
	vec.rzcw(spinAxis.getLength()*tH); // rotate clockwise!!!!!!! 
				// (Asteroid rotates ccw -> in asterocentric system fixed vector3s rotated cw)
} // void shape::spin(vector3 &vec, double t)


inline double shape::flux2mag(double flux) throw(std::invalid_argument)
{
	if (flux==0.) throw std::invalid_argument("flux=0 in flux2mag encountered");
	return -26.74 - 2.5*log10(flux);
};


inline void shape::setEmissivity(double emissivity) throw(std::invalid_argument)
{
	if (emissivity<=0 || emissivity > 1)
		throw std::invalid_argument("Error in shape::setEmissivity: 0<emissivity<=1 expected!");
	this->emissivity = emissivity;
};


#endif // !defined(AFX_SHAPE_H__58D2ED80_92B4_47D9_99C2_804BB4447285__INCLUDED_)
