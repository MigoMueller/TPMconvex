// End member of the class hierarchy for asteroid shape models.
// Describes a convex shape model made up from vertices and triangular facets.
// Input can be: class ObjFile or class MikkoFile (Kaasalainen's syntax)

// Author: michael.mueller@dlr.de
// Date:   Nov. 19 2003

// body file: ConvexTriangulated.cpp

#if !defined(AFX_CONVEXTRIANGULATED_H__AF1A20B5_E3B1_4E4E_8595_99E45C138FC9__INCLUDED_)
#define AFX_CONVEXTRIANGULATED_H__AF1A20B5_E3B1_4E4E_8595_99E45C138FC9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include<sstream>
#include "TriangulatedShape.h"
#include "ObjFile.h"
#include "MikkoFile.h"


class ConvexTriangulated : public TriangulatedShape  
{
public:
	ConvexTriangulated(const ObjFile& obj,	const SpinState& spinAxis,
				double pv, double H, double G=0.15)
		throw(std::bad_alloc, std::invalid_argument);

	ConvexTriangulated(MikkoFile& mikko, const SpinState& spinAxis,
				double pv, double H, double G=0.15)
		throw(std::bad_alloc, std::invalid_argument);

	virtual void plotGeometricLC (const char* const filename,
								  const eclipticVector& Sun2Aster,
								  const eclipticVector& Earth2Aster,
								  const unsigned int plotPoints = 40,
								  unsigned int periods=1)
								  const
								  throw (std::invalid_argument);
	virtual void plotMikkoLC     (const char* const filename,
								  double LS, double Lambert,
								  const eclipticVector& Sun2Aster,
								  const eclipticVector& Earth2Aster,
								  const unsigned int plotPoints = 40,
								  unsigned int periods=1)
								  const
								  throw (std::invalid_argument);
	virtual void plotThermalLag  (const char* const filename,
								  const eclipticVector& Sun2Aster,
								  const eclipticVector& Earth2Aster,
								  double lambdaMu, // wavelength in microns
								  double openingAngle, // opening angle of hemisphere's segments (deg)
								  double craterDensity, // 0<=craterDensity<=1
								  const unsigned int plotPoints=40,
								  unsigned int periods=1)
								  const
								  throw (std::invalid_argument);

	virtual void plotNEATM_POV (const char* const filename, 
								const eclipticVector& Sun2Aster, 
								const eclipticVector& Earth2Aster, 
								double eta, double rotationalPhase,
								int nColors) 
								const
								throw (std::invalid_argument, std::bad_alloc);
	virtual ~ConvexTriangulated() throw();

private:
	inline void fillArrays(const ObjFile& obj) throw();
public:
	// Input: measured data (fitFile), model fluxes for one albedo as produced by outputFluxes 
	//			(list of opening angles MUST include 0 as first entry!)
	// Ouput: File of the specified name, containing, for each non-zero opening angle, 
	//		  the best-fit crater density + ChiSquared
	static void calculateChiSquared(const fitFile& data, const char modelFluxes[], const char outName[]);
//	inline static void calculateChiSquared(const fitFile& data, 
//		const std::stringstream& modelFluxes, const std::stringstream& outName)
//		{calculateChiSquared (data, modelFluxes.str().c_str(), outName.str().c_str());};
	// outputs list of fluxes using LagFlux3 -> albedo constant, but craters' opening angle varies
	void outputFluxes (const char outName[], const fitFile& ephemerides, 
		unsigned int nOpeningAngles, const double openingAngles[], double accuracy=0.01) const;
//	inline void outputFluxes (const std::stringstream& outName, const fitFile& ephemerides,
//		unsigned int nOpeningAngles, const double openingAngles[], double accuracy=0.01) const
//		{outputFluxes (outName.str().c_str(), ephemerides, nOpeningAngles, openingAngles, accuracy);};
	// returns thermal flux assuming 100% crater density
	// uses LagerrosLookedup
	double LagFluxLookedup(const eclipticVector &Sun2Aster, const eclipticVector &Earth2Aster, 
					       double JD, double lambdaMu, double openingAngle) const;
	double LagFluxLookedup(const vector3 &A2Sun, const vector3 &A2Earth, 
					       double lambdaMu, double rAU, double deltaKM, double openingAngle) const;

	
	// ???????????????????????????????????


	void plotNEATMfitFileFRM(int n, double lambdaMu[], const char outname[], 
		const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster, 
		double JD) const;

	
	// ???????????????????????????????????


	void FRMspectrum (int n, const double* const lambdaMu, double* const flux, 
		const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster, 
		double JD) const;

	// returns thermal flux for infinite thermal parameter (cf. FRM.dvi in 'My Documents\FRM')
	double FRMflux(const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster, 
		double JD, double lambdaMu) const;
	double FRMflux (const vector3& A2Sun, const vector3& A2Earth, double lambdaMu, double rAU, double deltaKM) const;

	// returns the thermal flux (W/m^2/mu) at lambdaMu WITHOUT craters, thermal inertia or anything
	double LambertFlux (const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster,
		double JD, double lambdaMu) const;
	// thermal flux (smooth) with unit vector3 A2Sun, A2Earth in corotating coordinate system
	double LambertFlux (const vector3& A2Sun, const vector3& A2Earth, 
		double lambdaMu, double rAU, double deltaKM) const;

	// returns the thermal flux (W/m^2/mu) at lambdaMu WITHOUT craters, thermal inertia or anything
	double NEATMFlux (const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster,
		double JD, double lambdaMu, double eta) const;
	// thermal flux (smooth) with unit vector3 A2Sun, A2Earth in corotating coordinate system
	double NEATMFlux (const vector3& A2Sun, const vector3& A2Earth, 
		double lambdaMu, double rAU, double deltaKM, double eta) const;

	void makeHomFile (const char* const outName, double startJD, 
					  const eclipticVector& Sun2Aster, 
					  const eclipticVector& Earth2Aster, 
					  const unsigned int plotPoints) const;
protected:
	double NEATMFacet (int facet, const vector3& A2Sun, const vector3& A2Earth, 
		double lambdaMu, double rAU, double eta) const;
	double planckXU (double X,  double u4) const;
	double LambertFacet(int i, const vector3& A2Sun, const vector3&  A2Earth, double lambdaMu, double rAU) const;
	double FRMFacet(int facet, const vector3& A2Sun, const vector3& A2Earth, double lambdaMu, double rAU) const;
	double GeometricFacet (int i, vector3& A2Sun, vector3& A2Earth) const;
	double LSFacet (int i, vector3& A2Sun, vector3& A2Earth) const;
	double LambertFacetOptical(int i, vector3& A2Sun, vector3& A2Earth) const;
};

#endif // !defined(AFX_CONVEXTRIANGULATED_H__AF1A20B5_E3B1_4E4E_8595_99E45C138FC9__INCLUDED_)
