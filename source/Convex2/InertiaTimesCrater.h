// Author: michael.mueller@dlr.de
// Date:   2004 Aug 30

// Simplistic treatment of thermal inertia and craters:
// Calculate both separately, then multiply corresponding correction factors.

// Should be OK as 1st order approximation.

// Night-side: only thermal inertia, no craters.


#if !defined(AFX_INERTIATIMESCRATER_H__7A5B3B82_6C27_4ACF_91F4_31226285EF55__INCLUDED_)
#define AFX_INERTIATIMESCRATER_H__7A5B3B82_6C27_4ACF_91F4_31226285EF55__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "HemisphericNoInertia.h"
#include "ThermalInertiaOnlyConvex.h"

class InertiaTimesCrater : 
	public TriangulatedConvex::ThermalModelConvex
{
public:
	virtual ~InertiaTimesCrater()throw(){};

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

#endif // !defined(AFX_INERTIATIMESCRATER_H__7A5B3B82_6C27_4ACF_91F4_31226285EF55__INCLUDED_)
