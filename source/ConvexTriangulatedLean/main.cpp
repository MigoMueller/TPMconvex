#include "ConvexTriangulated.h"
//#include "NEATMFile.h"
#include "ObjFile.h"
#include<vector>
#include<sstream>
#include "LagerrosConversions.h"
using namespace std;

int main()
{

	const double lambda0     = 4;
	const double lambda1     = 12;
	const int    lambdaSteps = 81;
	vector<double> lambdaMu(lambdaSteps);
	{
		const double lambdaWidth = (lambda1-lambda0)/(lambdaSteps-1);
		for (int i=0; i<lambdaSteps; i++)
			lambdaMu[i] = lambda0 + lambdaWidth*i;
	}; 

	const ObjFile obj("y:/spheres/cube5.obj");
	const SpinState axis(0, 90, 24); // dummy
	const double pV=0.106;     // according to NASA's fact sheet
	const double H =-0.387739; // to set the diameter to 4880km (NASA again)
	const double G = 1.217;    // matches the Bond-albedo

	const ConvexTriangulated aster(obj, axis, pV, H, G);

	const eclipticVector SunMay8 (0, 0, 0.3659876273);
	// r as in May8 1995 -- Earth will be set to (phase, 0, deltaAU)

	const double phases[4] = {0, 45, 94.8, 120};



	const int openingAngle=175; // degrees
	const int thetaBar    =20;  // degrees
	//const double density=0.59985159741381466315514778687818;// --> thetabar = 8deg
	//double thetaBar = LagerrosConversions::HapkeRoughness(density, openingAngle);
	const double density = LagerrosConversions::DensityForThetabar(openingAngle, thetaBar);

	if (density<0 || density>1)
	{
		cerr<<"Density = "<<density<<endl;
		cerr<<"Bail out -- you better change thetabar or the opening angle!"<<endl;
		return -1;
	};

	ofstream out;
	{
		ostringstream dummy;
		dummy<<"D:/mercury/thetabar="<<thetaBar<<".gamma"<<openingAngle<<".dat";
		out.open(dummy.str().c_str());
	};
	if (!out)
		return -1;	
	//ostream& out=cout;


	out<<"# spectrum flat / rough / correction factor\n"
	   <<"# for different phase angles but with identical roughness properties:\n"
	   <<"# craters' opening angle = "<<openingAngle<<", density = "<<density<<'\n'
	   <<"# resulting in a thetabar of "<<LagerrosConversions::HapkeRoughness(density, openingAngle)<<'\n';

	for (int k=0; k<4; k++)
	{ 
		const double phase=phases[k];
		out<<'\n';
		out<<"# phase angle = "<<phase<<'\n';
		const eclipticVector EarMay8 (phase, 0, 0.9104204825);
		cerr<<"Phase:      "<<vector3::angle(SunMay8, EarMay8)<<endl;
		//cerr<<"Elongation: "<<vector3::angle(EarMay8-SunMay8, EarMay8)<<endl;
		{
			//out<<"lambda/mu ; flat/(W/m^2/mu) ; rough/(W/m^2/mu) ; rough/flat"<<endl;
			for (int i=0; i<lambdaSteps; i++)
			{
				double flat = aster.LambertFlux     (SunMay8, EarMay8, 0, lambdaMu[i]);
				double rough= density*aster.LagFluxLookedup (SunMay8, EarMay8, 0, lambdaMu[i], openingAngle)
					+(1-density)*flat;
				
				out<<fixed<<lambdaMu[i]<<'\t'
					<<scientific<<flat<<'\t'
					<<rough<<'\t'<<fixed
					<<rough/flat<<'\n';
			};
		};
	};
	
		
	return 0;
};