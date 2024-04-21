// LagerrosConversions.h: interface for the LagerrosConversions class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LAGERROSCONVERSIONS_H__04BBCE46_7B42_416B_B4AE_AB82D4FC1FFB__INCLUDED_)
#define AFX_LAGERROSCONVERSIONS_H__04BBCE46_7B42_416B_B4AE_AB82D4FC1FFB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<stdexcept>

class LagerrosConversions  
{
public:
	// calculates the crater density needed to get a Hapke thetabar
	// with Buhl-type craters (segments of hemispheres) of specified opening angle
	static double DensityForThetabar (double openingAngleDeg, double thetabarDeg);
	// returns Hapke's thetabar for given values of craters' opening angle and crater density
	static double HapkeRoughness(double density, double openingAngleDeg);

	LagerrosConversions();
	virtual ~LagerrosConversions();
	

};

#endif // !defined(AFX_LAGERROSCONVERSIONS_H__04BBCE46_7B42_416B_B4AE_AB82D4FC1FFB__INCLUDED_)
