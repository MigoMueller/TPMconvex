// Abstract basis class
// Describes craters in hemispherical shapes as introduced by Buhl, continued (e.g.) by Lagerros.
// Uses mainly Lagerros' parametrizations + methods.

// One improvement over Lagerros is incorporated (unpublished, so far):
// Multiple reflection is treated to *all* orders instead of to first order only -
// this is achieved my re-defining the member variable multipleReflection0.
// For the maths, see Lagerros.dvi in "My documents".


#if !defined(AFX_HemisphericCrater_H__11E065EE_231D_4EC7_BF6F_FE06121F6E67__INCLUDED_)
#define AFX_HemisphericCrater_H__11E065EE_231D_4EC7_BF6F_FE06121F6E67__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<cmath>
#include "vector3.h"

class HemisphericCrater  
{	
protected:
	HemisphericCrater (double emissivity, double uncorrected_bond, double gammaDeg);
	const double gammaDeg; // crater's opening angle in degrees (0<gamma<180)
	const double gammaRad; // in rad
	const double cosGamma2, sinGamma2, cosGamma4, sinGamma4, cosGamma, sinGamma;
	const double S;
	const double multipleReflection0;    // without mu (set in ctor)
	mutable double multipleReflection;   // including mu (set when integrating)

	// returns whether at the location of unitNormal unitOutbound is visible
	// both vectors must be of unit norm and point outwards
	inline bool sees(const vector3& unitNormal, const vector3& unitOutbound) const;

	// BAUSTELLE: replace by lookup-file!!!
	virtual inline double planck (double u4) const
		{return 1. / (exp(X/pow(u4, 0.25)) - 1);};

	mutable vector3 A2Earth, A2Sun;
	// must be set in deriving class (as soon as lambda and TSS are known...)!
	// X = h*c/(k*TSS*lambda)
	mutable double X;

	static double Simpson (const double oldTrapez, const double newTrapez)
		{return (4*newTrapez - oldTrapez)/3.;};
public:
	virtual ~HemisphericCrater()throw(){};
private:
	HemisphericCrater();
	HemisphericCrater(const HemisphericCrater&);
	HemisphericCrater& operator= (const HemisphericCrater&);
};



inline bool HemisphericCrater::sees(const vector3& unitNormal, const vector3& unitOutbound) const
//	{return 2*unitNormal.dot(unitOutbound) * unitOutbound.z <= cosGamma2 - unitNormal.z;};
	{return 2*unitNormal.dot(unitOutbound) * unitOutbound.z >= unitNormal.z - cosGamma2;};

#endif // !defined(AFX_HemisphericCrater_H__11E065EE_231D_4EC7_BF6F_FE06121F6E67__INCLUDED_)
