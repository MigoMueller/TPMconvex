

#if !defined(AFX_THERMALINERTIAONLYCONVEX_H__43E6342E_48AF_45CB_8457_900EC030B3BA__INCLUDED_)
#define AFX_THERMALINERTIAONLYCONVEX_H__43E6342E_48AF_45CB_8457_900EC030B3BA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "TriangulatedConvex.h"
#include "ThermalInertia.h"
#include <vector>


class ThermalInertiaOnlyConvex : 
	public ThermalInertia, 
	public TriangulatedConvex::ThermalModelConvex  
{
public:
	virtual ~ThermalInertiaOnlyConvex()throw(){};
	
	ThermalInertiaOnlyConvex (double emissivity, double bond, 
		double ThermalParameter, int nTime, int nZ, double zMax, double accuracyGoal)
		: ThermalInertia(ThermalParameter, nTime, nZ, zMax, accuracyGoal),
		TriangulatedConvex::ThermalModelConvex(emissivity, bond)
	{};

	ThermalInertiaOnlyConvex (const TriangulatedConvex& shape,
		double ThermalParameter, int nTime, int nZ, double zMax, double accuracyGoal)
		: ThermalInertia(ThermalParameter, nTime, nZ, zMax, accuracyGoal),
		TriangulatedConvex::ThermalModelConvex(shape)
	{};

//	virtual double fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
//		double lambdaMu, double rAU) const;
	virtual std::vector<double> fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const std::vector<double>& lambdaMu, double rAU) const;

	virtual vectorN ThermalLightCurveModFactors(const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		double lambdaMu, double rAU, unsigned int nPointsEach) const;

	// calculate surface temperature / TSS
	double u (const vector3& A2Sun, const vector3& dA) const;

	// calculate contribution to Yarkovsky recoil force modulo factors
	vector3 Yarkovsky_u4dA (const vector3& A2Sun, const vector3& dA) const;
protected:
	// auxiliary functions:
	inline double averageU () const;
	void solvePDE (const TimeStep& firstGuess) const;
	void cook(const vector3& A2Sun, const vector3& dA) const;
public:
	void calculate_inradiation (const vector3& A2Sun, const vector3& dA) const;
	vectorN uLightCurve(const vector3& A2Sun, const vector3& dA) const;

	virtual void setAccuracy (double accuracyGoal) const
		{ThermalModelConvex::setAccuracy(accuracyGoal);	ThermalInertia::setAccuracy(accuracyGoal);};
};

#endif // !defined(AFX_THERMALINERTIAONLYCONVEX_H__43E6342E_48AF_45CB_8457_900EC030B3BA__INCLUDED_)
