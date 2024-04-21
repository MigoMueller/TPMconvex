// NEATMFile.h
// Reads in and stores the content of a NEATM-fit file (as suitable for NEATMfit, too),
// is also capable of fitting eta and pV to it, 
// using ConvexTriangulated combined with an approximation to a sphere obj (static const).

// Author: michael.mueller@dlr.de
// Date:   2004 May 11

// not supposed to be base class to anything --> non-virtual dtor!


#if !defined(AFX_NEATMFILE_H__58AEDE22_4F46_4D92_A76F_BCAD5F629AE3__INCLUDED_)
#define AFX_NEATMFILE_H__58AEDE22_4F46_4D92_A76F_BCAD5F629AE3__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<map>
#include<fstream>
#include<string>
#include "ObjFile.h"
#include "ConvexTriangulated.h"


class NEATMFile  
{
public:
	// ctor reading in the file
	NEATMFile(const char*const fileName);
	~NEATMFile(){};

	// fits eta and pV to the file
	struct etaPV; 
	etaPV fit(double accuracy=0.01) const;
	struct etaPV
	{
		double eta, pV;
		etaPV(double eta, double pV): eta(eta), pV(pV){};
		etaPV average(const etaPV& rhs) const 
			{return etaPV (0.5*(eta+rhs.eta), 0.5*(pV+rhs.pV));};
		etaPV(): eta(0), pV(0){};
	}; // struct etaPV

	// parameters
	const double H, G, alphaDeg, rAU, deltaAU;

	// Wavelengths, fluxes, and errors are stored in the multimap 'data',
	// several occurances of the same wavelength are accepted!
	// Some typedefs before:

	struct ValueSigma
	{
		double value, sigma;
		ValueSigma(double value, double sigma) : value(value), sigma(sigma){};
		ValueSigma() : value(0),sigma(0) {};
	};

	// one data line -- can be directly added to 'data'
	struct dataLine
	{
		double lambdaMu, fluxSI, sigmaSI;
		dataLine (double lambdaMu, double fluxSI, double sigmaSI)
			: lambdaMu(lambdaMu), fluxSI(fluxSI), sigmaSI(sigmaSI) {};
		operator const std::pair<const double, ValueSigma> () const 
			{return std::pair<const double, ValueSigma> (lambdaMu, ValueSigma(fluxSI, sigmaSI));};
	};	

	typedef std::multimap<double, ValueSigma> fluxMap;

	// Here we go!!!
	const fluxMap data;

private:
	void skipComments (std::ifstream& in) const;
	NEATMFile();
	NEATMFile(const NEATMFile&);
	NEATMFile operator=(const NEATMFile&);
protected:
	// returns chi^2 for fix values of pV and eta
	double chi2 (ConvexTriangulated& asteroid, double pV, double eta) const;

	// For fixed eta, updates pV until its relative change is less than accuracy
	// Returns a fair approximation to chi^2 for the best-fit pV 
	double fitPV (ConvexTriangulated& asteroid, double eta, double& pV, double accuracy=0.1) const;
	inline double fitPV (ConvexTriangulated& asteroid, etaPV& parms, double accuracy=0.1) const
		{return fitPV(asteroid, parms.eta, parms.pV, accuracy);};

	const eclipticVector Sun2Aster, Earth2Aster;

	static const std::string sphereFile;
	static const ObjFile::proxy obj;
	static const double etaMin, etaStep;
	static const double etaMAX;
	static const double pvDefault;
	static const double pvMAX;

	// #steps taken after first coarse sweep thru parameter space (for fitting)
	static const int nRefinements;
};

#endif // !defined(AFX_NEATMFILE_H__58AEDE22_4F46_4D92_A76F_BCAD5F629AE3__INCLUDED_)
