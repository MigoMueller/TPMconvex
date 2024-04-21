
#include "HemisphericCrater.h"
#include<sstream>
using namespace std;

HemisphericCrater::HemisphericCrater(double emissivity, double uncorrected_bond, double gammaDeg)
:
gammaDeg(gammaDeg), 
gammaRad(gammaDeg*vector3::DEGREES),
cosGamma2(cos(gammaRad*0.5)),  sinGamma2(sin(gammaRad*0.5)),
cosGamma4(cos(gammaRad*0.25)), sinGamma4(sin(gammaRad*0.25)),
cosGamma (cos(gammaRad)),      sinGamma (sin(gammaRad)),
S(sinGamma4*sinGamma4),
multipleReflection0 ( (1-emissivity)*S*(1-S) / (1 - (1-emissivity)*S)), // improved treatment by myself
//multipleReflection0 ( (1-emissivity)*S*(1-S)), // Lagerros' treatment of multiple reflection
A2Earth(0,0,0), A2Sun(0,0,0)
{
	if (emissivity < 0 || emissivity > 1)
	{
		ostringstream dummy;
		dummy<<"Emissivity must range between zero and one, encountered "<<emissivity<<"!"<<ends;
		throw invalid_argument (dummy.str());
	};
	if (uncorrected_bond < 0 || uncorrected_bond > 1)
	{
		ostringstream dummy;
		dummy<<"Bond albedo must range between zero and one, encountered "<<uncorrected_bond<<'!'<<ends;
		throw invalid_argument (dummy.str());
	};
	if (gammaDeg < 0 || gammaDeg > 180)
	{
		ostringstream dummy;
		dummy<<"Craters' opening angle must range between 0 and 180 deg, encountered "<<gammaDeg<<'!'<<ends;
		throw invalid_argument (dummy.str());
	};
}; // HemisphericCrater::ctor
