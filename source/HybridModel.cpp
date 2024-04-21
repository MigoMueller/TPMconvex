#include "HybridModel.h"
using namespace std;

vector<double> HybridModel::fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
											const std::vector<double>& lambdaMu, double rAU) const
{
	if (craterDensity < 1e-5) // make boundary somewhat fuzzy
	{
		fluxPreliminaries(A2Sun, A2Earth, dA, lambdaMu, rAU);
		return fluxNoCraters(A2Sun, A2Earth, dA, lambdaMu, rAU);
	};
	if (craterDensity > 0.99999) // make boundary somewhat fuzzy
	{
		fluxPreliminaries(A2Sun, A2Earth, dA, lambdaMu, rAU);
		return fluxCraters(A2Sun, A2Earth, dA, lambdaMu, rAU);
	};

	vectorPair fluxes (this->bothFluxesModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU));
	vector<double>& noCraters = fluxes.first;
	vector<double>& Craters   = fluxes.second;

	const vector<double>::size_type N = lambdaMu.size();
	std::vector<double> result (N,0);
	for (unsigned int i=0; i<N; i++)
		result[i] = (1-craterDensity)*noCraters[i] + craterDensity*Craters[i];
	return result;
};
