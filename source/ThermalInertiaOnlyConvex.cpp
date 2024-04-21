// ThermalInertiaOnlyConvex.cpp: implementation of the ThermalInertiaOnlyConvex class.
//
//////////////////////////////////////////////////////////////////////


#include "ThermalInertiaOnlyConvex.h"
using namespace std;


std::vector<double> ThermalInertiaOnlyConvex::fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
															  const std::vector<double>& lambdaMu, double rAU) const

{
	const vector<double>::size_type N = lambdaMu.size();
	double cosEarthArea = dA.dot(A2Earth); 
	if (cosEarthArea < 0)
		return vector<double>(N, 0); // facet invisible

	double TSS = ThermalModel::TSS(rAU);
	double u = this->u(A2Sun, dA);
    if (u<1e-4) // prevent numerical underflow
        return vector<double>(N,0);
	vector<double> result(N);
	double X = AsteroidConstants::hck / (TSS*u);

	for (unsigned int i=0; i<N; i++)
		result[i] = cosEarthArea / (exp (X/lambdaMu[i]) -1);
	return result;
}; // ThermalInertiaOnlyConvex::fluxModFactors



vectorN ThermalInertiaOnlyConvex::ThermalLightCurveModFactors
		(const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		double lambdaMu, double rAU, unsigned int nPointsEach) const
{
	vectorN result(nTime); // ignore nPointsEach
	vectorN areaCosEarth(nTime);
	
	vector3 normal (dA);

	for (int i=0; i<nTime; i++)
	{
		areaCosEarth[i] = normal.dot(A2Earth);
		if (areaCosEarth[i] < 0)
			areaCosEarth[i] = 0;
		normal.rzcw(-dT);
	};

	if (areaCosEarth.isNull())
		return result; // facet invisible throughout entire revolution

	cook(A2Sun, dA);
	
	double TSS = ThermalModel::TSS(rAU);
	for (int i=0; i<nTime; i++)
	{
		if (areaCosEarth[i] == 0 || inRadiation[i] == 0)
			result[i]=0;
		else
		{
			result[i] = areaCosEarth[i] 
							 / (exp (AsteroidConstants::hck /(lambdaMu*TSS* cycle[i]() )) -1);
		};
	};
	return result;
}; // ThermalInertiaOnlyConvex::ThermalLightCurveModFactors



void ThermalInertiaOnlyConvex::solvePDE(const TimeStep &firstGuess) const
{
	// calculate temperature profile at point of interest using firstGuess as input
	cycle[0].calculate (firstGuess, inRadiation[0]);

	// calculate other points on LC
	for (int i=1; i<nTime; i++)
		cycle[i].calculate (cycle[i-1], inRadiation[i]);

	// spin asteroid until temperature at point of interest stable
	double oldU;
	do
	{
		oldU = cycle[0](); // last surface temperature at point of interest
		cycle[0].calculate(cycle[nTime-1], inRadiation[0]);

		for (int i=1; i<nTime; i++)
			cycle[i].calculate (cycle[i-1], inRadiation[i]);
	}
	while (fabs( (cycle[0]()-oldU)/cycle[0]()  ) > ThermalInertia::accuracyGoal);
}; // ThermalInertiaOnlyConvex::solvePDE



inline double ThermalInertiaOnlyConvex::averageU() const
{	
	double energy = 0;
	for (int i=0; i<nTime; i++)
//		energy+=inRadiation[i];		
//	return sqrt(sqrt(energy/nTime));
		energy+=sqrt(sqrt(inRadiation[i]));		
	return energy/nTime;

}; // ThermalInertiaOnlyConvex::averageU



void ThermalInertiaOnlyConvex::calculate_inradiation(const vector3 &A2Sun, const vector3 &dA) const
{
	vector3 unitNormal = dA/dA.modulus();
	for (int i=0; i<nTime; i++)
	{
		inRadiation[i] = unitNormal.dot(A2Sun);
		if (inRadiation[i]<0)
			inRadiation[i] = 0;
		unitNormal.rzcw(-dT); // in astero-centric frame, forward in time means: clockwise.
		                      // here, however, we DON'T rotate the Sun, rather the surface element itself
		                      // ---> have to rotate it counter-clockwise!
	};
}; // ThermalInertiaOnlyConvex::calculate_inradiation



double ThermalInertiaOnlyConvex::u(const vector3 &A2Sun, const vector3 &dA) const
{
	cook(A2Sun, dA);
	return cycle[0]();
}; // ThermalInertiaOnlyConvex::u



vector3 ThermalInertiaOnlyConvex::Yarkovsky_u4dA (const vector3& A2Sun, const vector3& dA) const
{
	double u = this->u(A2Sun, dA);
	return u*u*u*u*dA;
};


/*
double ThermalInertiaOnlyConvex::fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
	double lambdaMu, double rAU) const
{	
	double cosEarthArea = dA.dot(A2Earth); 
	if (cosEarthArea <= 0)
		return 0; // facet invisible

	// fill the std::vector inRadiation with solar flux >= 0
	calculate_inradiation(A2Sun, dA);
	// calculate average temperature
	double T = averageU();
	// make 2nd element constant temperature distribution -- will soon be overwritten
	cycle[1] = T; 
	// iteratively solve the PDE, using the constant distribution as first guess
	solvePDE(cycle[1]);

	double TSS = ThermalModel::TSS(rAU);

	double Planck_exponent = AsteroidConstants::hck / (lambdaMu*TSS*cycle[0]());
	return cosEarthArea / (exp(Planck_exponent)-1);
};
*/
/*
double ThermalInertiaOnlyConvex::fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
	double lambdaMu, double rAU) const
{
	vector<double> dummy(1,lambdaMu);
	//dummy.push_back(lambdaMu);
	vector<double> result = fluxModFactors(A2Sun, A2Earth, dA, dummy, rAU);
	return result[0];
}; // ThermalInertiaOnlyConvex::fluxModFactors 
*/



void ThermalInertiaOnlyConvex::cook(const vector3& A2Sun, const vector3& dA) const
{
	// fill the std::vector inRadiation with solar flux >= 0
	calculate_inradiation(A2Sun, dA);
	// calculate average temperature
	double T = averageU();

	if (T <= 1e-4) // prevent numerical trouble for T==0
	{
		for (int i=0; i<nTime; i++)
			cycle[i]() = 0; // set surface temperatures zero
		return;
	};

	// make 2nd element constant temperature distribution -- will soon be overwritten
	cycle[1] = T; 
	// iteratively solve the PDE, using the constant distribution as first guess
	solvePDE(cycle[1]);

}; // ThermalInertiaOnlyConvex::cook


vectorN ThermalInertiaOnlyConvex::uLightCurve(const vector3& A2Sun, const vector3& dA) const
{
	cook(A2Sun, dA);
	vectorN result(nTime);
	for (int i=0; i<nTime; i++)
		result[i] = cycle[i]();	
	return result;
};
