// Author: michael.mueller@dlr.de
// Date:   Nov. 15 2004

// Approximative treatment of thermal inertia in hemispheric craters
// Approximation used: cf. Lagerros A&A 332 1123-1132 (1998), eqn. 23:
// 1. Determine temperature on flat surface with and without thermal conduction
// 2. Calculate ratio w/wo conduction
// 3. Multiply temperature inside the crater with that ratio 
//
// Method used here: Instead of multiplying the temperature, multiply the wavelength 
// This is equivalent inside the integral over the Planck equation proper:
// 1/(exp( something/(lambda*T) )-1)
// Of course, the pre-factor lambda^{-5} must NOT be multiplied!

#if !defined(AFX_INERTIATIMESCRATERLAGERROS_H__D0F0A32F_E463_44EE_AC59_B39F537E36F7__INCLUDED_)
#define AFX_INERTIATIMESCRATERLAGERROS_H__D0F0A32F_E463_44EE_AC59_B39F537E36F7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "HybridModel.h"
#include "HemisphericNoInertia.h"
#include "ThermalInertiaOnlyConvex.h"

class InertiaTimesCraterLagerros : 
//	public TriangulatedConvex::ThermalModelConvex
	public HybridModel
{
public:
	virtual ~InertiaTimesCraterLagerros()throw(){};

	InertiaTimesCraterLagerros (double emissivity, double uncorrected_bond, 
		double ThermalParameter, int nTime, int nZ, double zMax, double accuracyGoal,
		double gammaDeg,
		double craterDensity = 1)
	:
		//ThermalModelConvex (emissivity, uncorrected_bond),
		HybridModel(emissivity, uncorrected_bond, craterDensity),
		craters (emissivity, uncorrected_bond, gammaDeg),
		inertia (emissivity, uncorrected_bond, ThermalParameter, nTime, nZ, zMax, accuracyGoal),
		dummyReturnVector(0)
	{};

	InertiaTimesCraterLagerros (const TriangulatedConvex& aster,
		double ThermalParameter, int nTime, int nZ, double zMax, double accuracyGoal,
		double gammaDeg,
		double craterDensity = 1)
	:
		HybridModel(aster, craterDensity),
		craters (aster.getEmissivity(), aster.getBond(), gammaDeg),
		inertia (aster.getEmissivity(), aster.getBond(), ThermalParameter, nTime, nZ, zMax, accuracyGoal),
		dummyReturnVector(0)
	{};
    void setAccuracies(double accuracyGoalInertia, double accuracyGoalCrater)
    {
        inertia.setAccuracy(accuracyGoalInertia);
        craters.setAccuracy(accuracyGoalCrater);
    };
protected:
//	virtual double fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
//			double lambdaMu, double rAU) const;
//	virtual std::vector<double> fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
//  const std::vector<double>& lambdaMu, double rAU) const;
	virtual void fluxPreliminaries (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const;
	virtual std::vector<double> fluxNoCraters (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const;
	virtual std::vector<double> fluxCraters (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const;
	mutable double areaCosEarth, cosSun, u_inertia;
	mutable std::vector<double>::size_type N;
	
	// dummy used when facet on nightside:
	// fluxNoCraters will save its result here if cosSun<0, fluxCraters will copy from here and delete afterwards
	mutable std::vector<double>* dummyReturnVector; // for nightside

    // not properly implemented; don't use!
    virtual vectorN ThermalLightCurveModFactors(const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		double lambdaMu, double rAU, unsigned int nPointsEach) const;
 //   {throw std::logic_error(
	//"Routine InertiaTimesCraterLagerros::ThermalLightCurveModFactors should not have been called!");};

protected:
	HemisphericNoInertia      craters;
	ThermalInertiaOnlyConvex  inertia;
};

#endif // !defined(AFX_INERTIATIMESCRATERLAGERROS_H__D0F0A32F_E463_44EE_AC59_B39F537E36F7__INCLUDED_)
