#include "HemisphericNoInertia.h"
#include<sstream>
using namespace std;


const int HemisphericNoInertia::maxTrapez = 10;


HemisphericNoInertia::HemisphericNoInertia (double emissivity, double uncorrected_bond, double gammaDeg)
: 
HemisphericCrater (emissivity, uncorrected_bond, gammaDeg),
TriangulatedConvex::ThermalModelConvexNoInertia (emissivity, uncorrected_bond),
TempTerm0 ( S * ( (1-S)*bond+emissivity ) / (1-bond*S)  )
{
	init();
}; // HemisphericNoInertia::ctor



HemisphericNoInertia::HemisphericNoInertia (const TriangulatedConvex& shape, double gammaDeg)
:
HemisphericCrater (shape.getEmissivity(), shape.getBond(), gammaDeg),
TriangulatedConvex::ThermalModelConvexNoInertia(shape),
TempTerm0 ( S * ( (1-S)*bond+emissivity ) / (1-bond*S)  )
{ 
	init();
}; // HemisphericNoInertia::ctor



inline void HemisphericNoInertia::init()
{}; // HemisphericNoInertia::init



inline double HemisphericNoInertia::planck (const vector3& unitNormal) const
{
	if (sees(unitNormal, A2Sun))
		return HemisphericCrater::planck( unitNormal.dot(A2Sun) + TempTerm);
	else
		return planck0;
}; // HemisphericNoInertia::planck (const vector3&)



double HemisphericNoInertia::fluxModFactorsPerArea (double cosSun, double cosEarth,
		                                            double cosAzi, double lambdaMu, double rAU) const
{
// check whether there's anything to do at all:
	if (cosSun < 0.0004) // don't calculate anything for Sun ridiculously low
		return 0;
	if (cosEarth < 0.0004) // same with Earth
		return 0;
	if (lambdaMu <= 0)
	{
		ostringstream dummy;
		dummy<<"Encountered a wavelength of "<<lambdaMu<<" microns -- must be nonnegative!!!"<<ends;
		throw logic_error (dummy.str());
	};
	if (rAU <= 0)
	{
		ostringstream dummy;
		dummy<<"Encountered a heliocentric distance of "<<rAU<<" AU -- must be nonnegative!!!"<<ends;
		throw logic_error(dummy.str());
	};
// check cosines for numerical noise, correct if necessary
	if (fabs(cosAzi)>1)
	{
		if (fabs (fabs(cosAzi)-1) < 1e-4) // |cos| = 1 +- epsilon --> numeric noise!
		{
			if (cosAzi > 0)
				cosAzi = 0.99999;
			else
				cosAzi = -0.99999;
		}
		else
		{
			ostringstream dummy;
			dummy<<"Encountered cosAzimuth of "<<cosAzi<<", must be element of [-1,1]!"<<ends;
			throw logic_error (dummy.str());
		};
	};
	if (cosEarth > 1)
	{
		if (cosEarth > 1.00004)
		{
			ostringstream dummy;
			dummy<<"Encountered cosEarth of "<<cosEarth<<", must be <= 1!"<<ends;
			throw logic_error (dummy.str());
		};
		cosEarth = 0.99999;
	};
	if (cosSun > 1)
	{
		if (cosSun > 1.00004)
		{
			ostringstream dummy;
			dummy<<"Encountered cosSun of "<<cosSun<<", must be <=1!"<<ends;
			throw logic_error (dummy.str());
		};
		cosSun = 0.99999;
	};

// done with checking; set constants:
	{
		double sinSun   = sqrt(1-cosSun*cosSun);
		double sinEarth = sqrt(1-cosEarth*cosEarth);
		double sinAzi   = sqrt(1-cosAzi*cosAzi);

		A2Sun   = vector3(sinSun, 0, cosSun);
		A2Earth = vector3(sinEarth * cosAzi, sinEarth*sinAzi, cosEarth);
	};
	X = AsteroidConstants::hck/(lambdaMu*TSS(rAU));
	TempTerm = TempTerm0 * cosSun;
	multipleReflection = multipleReflection0 * cosEarth;
	planck0 = HemisphericCrater::planck(TempTerm);

// done setting constants, start actual computations:
	// estimates for integral / (2Pi * 2S)
	// 2Pi from azimuth-integration, 2S from integral over cos(theta=polar angle)
	double oldTrapez=0, newTrapez=0;
	double oldSimpson=0, newSimpson=0;

	double limit = 1-2*S;     // integration over cosTheta from 1-2S .. 1
	if (limit < 0) limit = 0; // just in case
	
	// 1st estimate by hand:
	{ // contribution from bottom
		double m = cosEarth;
		if (m*m < S) //bottom invisible
			m = 0;
		if (cosSun*cosSun < S) // bottom shadowed
			oldTrapez = planck0 * (m + multipleReflection);
		else
			oldTrapez = HemisphericCrater::planck (cosSun + TempTerm) * (m + multipleReflection);
	};
	oldTrapez += phiIntegral (limit); // contribution from rim
	oldTrapez *= 0.5;

	double start = 1-S;			// middle point between 1-2S and 1
	double cosStep   = 2*S;		// initial step size
	unsigned long nPoints = 1;	// all these three parameters will be changed in trapezTheta
								// fit automatically for next iteration!

	newTrapez = trapezTheta (oldTrapez, start, cosStep, nPoints);
	oldSimpson = Simpson(oldTrapez, newTrapez);

	oldTrapez=newTrapez;
	newTrapez = trapezTheta(oldTrapez, start, cosStep, nPoints);
	newSimpson = Simpson(oldTrapez, newTrapez);

	for (int i=2; i<maxTrapez && fabs(newSimpson-oldSimpson)/newSimpson > accuracyGoal; i++)
	{
		oldTrapez = newTrapez;
		oldSimpson = newSimpson;
		newTrapez = trapezTheta(oldTrapez, start, cosStep,nPoints);
		newSimpson = Simpson(oldTrapez, newTrapez);
	};
	return newSimpson / (1-S); 
	// integral = newSimpson * 2Pi * 2S
	// area     = 4Pi * S * (1-S)
	// return integral/area
} // HemisphericNoInertia::fluxModFactorsPerArea



double HemisphericNoInertia::phiIntegral(double cosTheta) const
{
	const double sinTheta = sqrt(1-cosTheta*cosTheta);
	
	double oldTrapez = integrand (vector3(-sinTheta, 0, cosTheta));
	oldTrapez+= integrand (vector3( sinTheta, 0, cosTheta));
	oldTrapez+= integrand (vector3(0, -sinTheta, cosTheta));
	oldTrapez+= integrand (vector3(0,  sinTheta, cosTheta));
	oldTrapez *= 0.25; // trapez with four inner points

	double cosStep = 0; // phiStep = PI/2
	double sinStep = 1; 
	unsigned long nPoints = 4; // four new points to be added
	double newTrapez = trapezPhi (nPoints, oldTrapez, cosStep, sinStep, cosTheta, sinTheta);
	double oldSimpson = Simpson(oldTrapez, newTrapez);

	oldTrapez = newTrapez;
	newTrapez = trapezPhi (nPoints, oldTrapez, cosStep, sinStep, cosTheta, sinTheta);
	double newSimpson = Simpson(oldTrapez, newTrapez);

	if (newSimpson == 0)
		return 0; // suppose everything is invisible -- prevents crashing due to dividing by zero!
	
	for (int i=2; i<maxTrapez && fabs(newSimpson-oldSimpson)/newSimpson > accuracyGoal; i++)
	{
		oldTrapez=newTrapez;
		oldSimpson=newSimpson;
		newTrapez = trapezPhi (nPoints, oldTrapez, cosStep, sinStep, cosTheta, sinTheta);
		newSimpson = Simpson(oldTrapez, newTrapez);
	};
	return newSimpson;
} // HemisphericNoInertia::phiIntegral



inline double HemisphericNoInertia::integrand(const vector3 &unitNormal) const
{
	double m;
	if (sees(unitNormal, A2Earth))
		m = A2Earth.dot(unitNormal);
	else
		m=0;
	return planck(unitNormal) * (m + multipleReflection);
} //HemisphericNoInertia::integrand



// nPoints = #points so far = #points to be added
// phiStep = 2Pi / nPoints
// cosStep = cos(phiStep), sinStep = sin(phiStep)

// will be updated for next iteration!!!
double HemisphericNoInertia::trapezPhi(unsigned long &nPoints, 
									   double oldTrapez, double &cosStep, double &sinStep,
									   double cosTheta, double sinTheta) const
{
	double cosPhi0=sqrt( (1+cosStep)*0.5 ); // cos(phi0) = cos (phiStep/2)
	double sinPhi0=sqrt( (1-cosStep)*0.5 ); // sin(phi0) = sin (phiStep/2)

	double value=0, dummy;
	double cosPhi=cosPhi0, sinPhi=sinPhi0;
	for (unsigned int i=0; i<nPoints; i++)
	{
		value += integrand(vector3(-sinTheta*cosPhi, -sinTheta*sinPhi, cosTheta));
		dummy  = sinPhi*cosStep + cosPhi*sinStep;
		cosPhi = cosPhi*cosStep - sinPhi*sinStep;
		sinPhi = dummy;
	}
	value/=nPoints; 
	value += oldTrapez;

	cosStep=cosPhi0;		// step size is now half of old step size
	sinStep=sinPhi0;
	nPoints *= 2;

	return value*0.5;
} //HemisphericNoInertia::trapezPhi



double HemisphericNoInertia::trapezTheta(double oldTrapez, 
										 double &start, double &stepSize, unsigned long &nPoints) const
{
	double value = 0;
	double cosTheta=start;
	for (unsigned long i=0; i<nPoints; i++)
	{
		//sinTheta=sqrt(1-cosTheta*cosTheta);
		value += phiIntegral(cosTheta);
		cosTheta+=stepSize;
	}
	value/=nPoints;
	value+=oldTrapez; 
	value*=0.5;

	stepSize *= 0.5;				// new stepsize is half of the old stepsize
	start    -= 0.5*stepSize;		// new start point is half a NEW step size smaller
	nPoints  *= 2;					// twice as many points, now

	return value;
}; // HemisphericNoInertia::trapezTheta
