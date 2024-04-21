#include "InertiaTimesCrater.h"
#include<iostream>
using namespace std;


vector<double> InertiaTimesCrater::fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
			const vector<double>& lambdaMu, double rAU) const 
{
	const vector<double>::size_type N = lambdaMu.size();
	if (dA.dot(A2Earth) < 0) // invisible
		return vector<double> (N,0);
	if (dA.dot(A2Sun) < 0) // night side -- forget about craters
		return inertia.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);
	else
	{
		vector<double> correct = inertia.CorrectionFactor(A2Sun, A2Earth, dA, lambdaMu, rAU);
		vector<double> flux    = craters.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);
		for (unsigned int i=0; i<N; i++)
		{
			if (correct[i]<10)
				flux[i] *= correct[i];
			else
				flux[i] *= 10;
		};
		return flux;
	};
}; // InertiaTimesCrater::fluxModFactors


/*
double InertiaTimesCrater::fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		double lambdaMu, double rAU) const 
{
	if (dA.dot(A2Earth) < 0) // invisible
		return 0;
	if (dA.dot(A2Sun) < 0) // night side -- no craters
		return inertia.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);
	else
		//return inertia.CorrectionFactor(A2Sun, A2Earth, dA, lambdaMu, rAU)
		//* craters.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);
	{
		//double correct = dynamic_cast<const TriangulatedConvex::ThermalModelConvex&>(craters)
			// funnily, VC++ needs that cast to compile - shouldn't be necessary, though
		//	.CorrectionFactor(A2Sun, A2Earth, dA, lambdaMu, rAU);
		//double flux    = inertia.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);

		double correct = inertia.CorrectionFactor(A2Sun, A2Earth, dA, lambdaMu, rAU);
		double flux    = craters.fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);
		if (correct > 10)
		{
			correct = 10;  // correct diverges at the terminator, leading to unphysical results
			//cerr<<"Clipped correction factor at cosSun = "<<A2Sun.dot(dA)/dA.modulus()<<endl;
		};
		// interestingly, it clips more often when the craters' correction factor is used 
		// (tested with thermal parameter = 0.9)

		return flux*correct;
	};
};
*/

