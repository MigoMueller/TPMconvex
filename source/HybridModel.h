// HybridModel
// Author: Michael.Mueller@dlr.de
// Date:   2005 June 14
//
// Abstract base class to describe thermal models with two sub-models
// (typically: one with and one without craters).
// Not meant to be instanciated, rather as a baseclass to, e.g., InertiaTimesCraterLagerros.
//
// Implementations should feature two thermal models (model0 and model1),
// which in the following will be called 'noCraters' and 'Craters'.
//
// fluxModFactors will mix the fluxes from both models, using 
// the variable (defined here) craterDensity:
// flux = (1-craterDensity)*flux(noCraters) + craterDensity*flux(Craters)
//
// subclasses must define fluxCraters and fluxNoCraters,
// with fluxPrelimininaries to possibly speed up the numerics.
// These are defined as 'const', so whatever fluxPreliminaries changes must be 'mutable'
//
// Additionally, HybridModel features the function bothFluxesModFactors,
// returning a std::pair of fluxes from both models.
// For the sake of simplicity, a type vectorPair is typedefed (=return type).

#if !defined (__Triangulated_Convex_H)
#define __Triangulated_Convex_H


//#pragma once
#include "TriangulatedConvex.h"

class HybridModel :
	public TriangulatedConvex::ThermalModelConvex
{
public:
	typedef std::pair<std::vector<double>, std::vector<double> > vectorPair;

	HybridModel(double emissivity, double bond, double craterDensity=1)
		: TriangulatedConvex::ThermalModelConvex(emissivity, bond),
		  craterDensity(craterDensity)
	{};
	HybridModel(const TriangulatedConvex& shape, double craterDensity=1)
		: TriangulatedConvex::ThermalModelConvex(shape),
		  craterDensity(craterDensity)
	{};
	virtual ~HybridModel() throw () {};

	virtual std::vector<double> fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const;
	
	// see explanation in header of file
	vectorPair bothFluxesModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const
	{
		fluxPreliminaries(A2Sun, A2Earth, dA, lambdaMu, rAU);
		std::vector<double> noCraters = fluxNoCraters(A2Sun, A2Earth, dA, lambdaMu, rAU);
		std::vector<double> Craters   = fluxCraters(A2Sun, A2Earth, dA, lambdaMu, rAU);
		return vectorPair(noCraters, Craters);
	};

protected:
	const double craterDensity;
	// to facilitate numerical speed up in subclass
	// maybe some calculations need to be done for both models? This is the place for these!
	// Note: defined as 'const', so whatever fluxPreliminaries changes must be 'mutable'.
	virtual void fluxPreliminaries (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const
		=0;
	virtual std::vector<double> fluxNoCraters (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const
		=0;
	virtual std::vector<double> fluxCraters (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const
		=0;
};

#endif
