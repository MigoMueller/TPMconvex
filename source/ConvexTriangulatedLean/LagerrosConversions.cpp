// LagerrosConversions.cpp: implementation of the LagerrosConversions class.
//
//////////////////////////////////////////////////////////////////////

#include "LagerrosConversions.h"
#include<cmath>
#include<sstream>
#include "vector3.h" // for PI
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

LagerrosConversions::LagerrosConversions()
{}


LagerrosConversions::~LagerrosConversions()
{}

double LagerrosConversions::HapkeRoughness(double density, double openingAngleDeg)
{
	if (density<0 || density>1)
	{
		ostringstream dummy;
		dummy<<"Invalid density "<<density<<": must be between 0 and 1!"<<ends;
		throw invalid_argument(dummy.str());
	};
	if (openingAngleDeg<=0 || openingAngleDeg >= 180)
	{
		ostringstream dummy;
		dummy<<"Invalid opening angle "<<openingAngleDeg<<": must be >=0 and <= 180 deg!"<<ends;
		throw invalid_argument(dummy.str());
	};

	double gamma2=openingAngleDeg*vector3::DEGREES*0.5; // gamma/2 in rad
	double sinGamma2=sin(gamma2);
	double cosGamma2=cos(gamma2);

	double tanThetaBar = 2*density/vector3::PI * (sinGamma2 + log(cosGamma2/(1+sinGamma2))) / (cosGamma2-1);
	return atan(tanThetaBar)*vector3::RAD;
}; // LagerrosConversion::HapkeRoughness



double LagerrosConversions::DensityForThetabar(double openingAngleDeg, double thetabarDeg)
{
	if (openingAngleDeg <=0 || openingAngleDeg >= 180)
	{
		ostringstream dummy;
		dummy<<"Invalid opening angle "<<openingAngleDeg<<": must be >=0 and <= 180 deg!"<<ends;
		throw invalid_argument(dummy.str());
	};

	if (thetabarDeg <=0 || thetabarDeg >= 180)
	{
		ostringstream dummy;
		dummy<<"Invalid thetabar "<<thetabarDeg<<": must be >=0 and <= 180 deg!"<<ends;
		throw invalid_argument(dummy.str());
	};

	const double gamma2=openingAngleDeg*vector3::DEGREES*0.5; // gamma/2 in rad
	const double sinGamma2=sin(gamma2);
	const double cosGamma2=cos(gamma2);

	// tan(thetabar) assuming density=100%
	const double tanThetaBar1 = 2/vector3::PI * (sinGamma2 + log(cosGamma2/(1+sinGamma2))) / (cosGamma2-1);
	// tan(thetabar) as wished by user
	const double tanThetaBar  = tan(thetabarDeg*vector3::DEGREES);

	return tanThetaBar/tanThetaBar1;
}; 
