#include "InertiaTimesCraterLagerros.h"
#include<iostream>
using namespace std;


void InertiaTimesCraterLagerros::fluxPreliminaries (
	const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
	const std::vector<double>& lambdaMu, double rAU) const
{
	N = lambdaMu.size();
	areaCosEarth = dA.dot(A2Earth);
	if (areaCosEarth <= 0) // facet invisible
		return; // both fluxCraters and fluxNoCraters will return 0
	cosSun = dA.dot(A2Sun)/dA.modulus();
	if (cosSun <= 0) // facet on night side
		return; // don't need u_inertia in this case
	u_inertia = inertia.u(A2Sun, dA);
};



vector<double> InertiaTimesCraterLagerros::fluxNoCraters (
	const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
	const std::vector<double>& lambdaMu, double rAU) const
{
	if (areaCosEarth <= 0)
		return vector<double> (N,0);
	if (cosSun <= 0)
	{
		// night-side, so fluxCraters won't calculate anything
		// save result here, so fluxCraters can copy it
		delete dummyReturnVector;
		dummyReturnVector = new vector<double> (inertia.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU));
		return *dummyReturnVector;
	}
	else
	{
		double X=AsteroidConstants::hck/(u_inertia*TSS(rAU));
		vector<double> result (N,0);
		for (unsigned int i=0; i<N; i++)
			result[i] = areaCosEarth / ( exp(X/lambdaMu[i]) - 1);
		return result;
	};
}; // InertiaTimesCraterLagerros::fluxNoCraters



vector<double> InertiaTimesCraterLagerros::fluxCraters (
	const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
	const std::vector<double>& lambdaMu, double rAU) const
{
	if (areaCosEarth <= 0)
		return vector<double> (N,0);
	if (cosSun <= 0)
	{
		// night side, so craters are disregarded
		if (craterDensity < 0.99999) // make boundary somewhat fuzzy
			// fluxNoCraters has done the job, already
			return *dummyReturnVector;
		else
			return inertia.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);
	}
	else
	{
		double u_Lambert = sqrt(sqrt( A2Sun.dot(dA)/dA.modulus()  )); // u^4 = cosSun
		double temperature_correction = u_inertia/u_Lambert;
		// prevent divergences from terminator:
		if (temperature_correction > 1.4)
		{
			if (cosSun > 0.1)
				cerr<<"Encountered large temperature correction ("<<100*(temperature_correction-1)
				<<"%) at cosSun = "<<cosSun<<endl;
			temperature_correction = 1.4; 
		};
		vector<double> result(N,0);
		for (unsigned int i=0; i<N; i++)
			// this is the central piece of this routine:
			// call crater-function with wavelength*temperature_correction
			result[i] = craters.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu[i]*temperature_correction, rAU);
		return result;
	};
}; // InertiaTimesCraterLagerros::fluxCraters



/*
vector<double> InertiaTimesCraterLagerros::fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		const vector<double>& lambdaMu, double rAU) const
{
	const int N = lambdaMu.size();
	if (dA.dot(A2Earth) < 0) // invisible
		return vector<double> (N,0);
	if (dA.dot(A2Sun) <= 0) // night side -- forget about craters
		return inertia.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);
	else
	{
		double u_inertia = inertia.u(A2Sun, dA);
		double u_Lambert = sqrt(sqrt( A2Sun.dot(dA)/dA.modulus()  )); // u^4 = cosSun
		double temperature_correction = u_inertia/u_Lambert;
		// prevent divergences from terminator:
		if (temperature_correction > 1.3)
		{
			double cosSun = pow(u_Lambert, 4);  //= dA.dot(A2Sun)/dA.modulus()
			if (cosSun > 0.1)
				cerr<<"Encountered large temperature correction ("<<100*(temperature_correction-1)
				    <<"%) at cosSun = "<<cosSun<<endl;
			temperature_correction = 1.3; 
		};
		vector<double> result(N,0);
		for (int i=0; i<N; i++)
			// this is the central piece of this routine:
			// call crater-function with wavelength*temperature_correction
			result[i] = craters.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu[i]*temperature_correction, rAU);
		return result;
	};
}; // InertiaTimesCraterLagerros::fluxModFactors
*/ // now taken care of by HybridModel


// not properly implemented, don't use!
vectorN InertiaTimesCraterLagerros::ThermalLightCurveModFactors
	(const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
	double lambdaMu, double rAU, unsigned int nPointsEach) const
{
	const unsigned int N = inertia.getNTime(); // disregard nPointsEach

	vector<vector3> normals;
	normals.reserve(N);
	normals.push_back(dA);
	for (unsigned int i=1; i<N; i++)
	{
		normals.push_back(normals.back());
		normals.back().rzcw(-inertia.getDT());
	};

	vectorN areaCosEarth(N);
	{
		for (unsigned int i=0; i<N; i++)
		{
			areaCosEarth[i] = A2Earth.dot(normals[i]);
			if (areaCosEarth[i] < 0)
				areaCosEarth[i] = 0;
		};
	};
	if (areaCosEarth.isNull()) // facet invisible throughout revolution
		return areaCosEarth;

	inertia.calculate_inradiation(A2Sun, dA);
	const vectorN& u4_Lambert = inertia.getInRadiation();       // lightcurve of u^4 (no inertia)

	if (u4_Lambert.isNull()) // facet shadowed throughout revolution
		return u4_Lambert;
	vectorN u_inertia = inertia.uLightCurve(A2Sun, dA);

	const double TSS = ThermalModel::TSS(rAU);
	vectorN result(N);
	for (unsigned int i=0; i<N; i++)
	{
		if (u4_Lambert[i] <= 0) // night side, forget craters
		{
			if (u_inertia[i] < 1e-4)
				result[i]=0;
			else
				result[i] = areaCosEarth[i] / (exp(AsteroidConstants::hck/(TSS*lambdaMu*u_inertia[i] ))-1);
		}
		else
		{
			double temperatureCorrection = u_inertia[i] / sqrt(sqrt(u4_Lambert[i]));
			if (temperatureCorrection > 1.3)
				temperatureCorrection = 1.3;
			if (craterDensity != 0)
			{
				result[i] = craters.fluxModFactors(A2Sun, A2Earth, normals[i], lambdaMu*temperatureCorrection, rAU);
				result[i] *= craterDensity;
				result[i] += (1.-craterDensity)*areaCosEarth[i] / (exp(AsteroidConstants::hck/(TSS*lambdaMu*u_inertia[i] ))-1);
			}
			else
				result[i] = areaCosEarth[i] / (exp(AsteroidConstants::hck/(TSS*lambdaMu*u_inertia[i] ))-1);
		}
	};

	return result;
};
