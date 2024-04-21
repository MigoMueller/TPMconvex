// LagerrosCraterInertia.h: interface for the LagerrosCraterInertia class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LAGERROSCRATERINERTIA_H__C761197E_6901_4794_9CAF_7F0B086B14F1__INCLUDED_)
#define AFX_LAGERROSCRATERINERTIA_H__C761197E_6901_4794_9CAF_7F0B086B14F1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "HemisphericNoInertia.h"
#include "ThermalInertiaOnlyConvex.h"

{
public:
	virtual ~InertiaTimesCrater(){};

	InertiaTimesCrater (double emissivity, double uncorrected_bond, 
		double ThermalParameter, int nTime, int nZ, double zMax, double accuracyGoal,
		double gammaDeg)
	:
		ThermalModelConvex (emissivity, uncorrected_bond),
		craters (emissivity, uncorrected_bond, gammaDeg),
		inertia (emissivity, uncorrected_bond, ThermalParameter, nTime, nZ, zMax, accuracyGoal)
	{};

//	virtual double fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
//			double lambdaMu, double rAU) const;
	virtual std::vector<double> fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
			const std::vector<double>& lambdaMu, double rAU) const;

protected:
	HemisphericNoInertia      craters;
	ThermalInertiaOnlyConvex  inertia;
};

class LagerrosCraterInertia : 
	public TriangulatedConvex::ThermalModelConvex  
{
public:
	virtual ~LagerrosCraterInertia(){};

	LagerrosCraterInertia (double emissivity, double uncorrected_bond, 
		double ThermalParameter, int nTime, int nZ, double zMax, double accuracyGoal,
		double gammaDeg)
		:
	ThermalModelConvex (emissivity, uncorrected_bond),
		craters (emissivity, uncorrected_bond, gammaDeg),
		inertia (emissivity, uncorrected_bond, ThermalParameter, nTime, nZ, zMax, accuracyGoal)
	{};


};

#endif // !defined(AFX_LAGERROSCRATERINERTIA_H__C761197E_6901_4794_9CAF_7F0B086B14F1__INCLUDED_)
