
// Checking for conservation of energy from disk-integrated flux

#include "TriangulatedConvex.h"
#include "HemisphericNoInertia.h"
#include "ThermalInertiaOnlyConvex.h"
#include "InertiaTimesCraterLagerros.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<ctime>
#include"jansky.h"
using namespace std;


const double lambda0    = 0.5; 
const double lambda1    = 79.5; 
const int lambdaSteps   = 80; // #points including lambda0 and lambda1
  // must be even!
const double lambdaWidth= (lambda1-lambda0)/(lambdaSteps-1);

const unsigned int nEarths = 500;

const double emissivity = 0.9;
const double pV         = .34;
const double H          = 10;
const double G          = 0.15;
const double periodH     = 6;
const double rAU        = 1;
const double deltaAU    = 1;
const double deltaKM    = deltaAU * AsteroidConstants::AU;
const double gammaDeg   = 0; // opening angle in degrees
const double Inertia    = 500;
const double accuracy_goal = 0.001;


// returns a random number from -1 -- 1
inline double random()
{
	return 2*rand() / (RAND_MAX + 1.) - 1;
}

int main()
{
	srand(static_cast<unsigned>(time(0))); // initialize random number generator
	vector<vector3> Earths;
	Earths.reserve(nEarths);
	{
		double x,y,z;
		for (int i=0; i<nEarths; i++)
		{
			// add new random unit vector to Earths
			x=random();
			y=random();
			z=random();
			Earths.push_back(vector3(x,y,z));
			Earths[i] /= Earths[i].modulus();
		}
	};
	//generated nEarths random orientations

	ConvexFile file ("S:/spheres/cube6.obj.convex");
	TriangulatedConvex aster (file, SpinState(0, 90, periodH), H, G, pV, emissivity);

	const vector3 A2Sun  (1,0,0);

	vector<double> lambdaMu;
	lambdaMu.resize(lambdaSteps);
	{
		double dummy = lambda0;
		for (int i=0; i<lambdaSteps; i++)
		{
			lambdaMu[i] = dummy;
			dummy+=lambdaWidth;
		};
	};

	ThermalInertiaOnlyConvex model (emissivity, aster.getBond(), 
		aster.calculateThermalParameter(Inertia, rAU),
		720, 33, 8, accuracy_goal);
	//InertiaTimesCraterLagerros model (emissivity, aster.getBond(), 
	//	aster.calculateThermalParameter(Inertia, rAU), 
	//	300, 25, 5, accuracy_goal, gammaDeg);
	//TriangulatedConvex::LambertianEmitter model(aster);

	// outmost integral: MC-integration around the asteroid
	// implemented as sum over random orientations, 
	// divided by #orientations thereafter
	double energy = 0;
	const int nLambda = lambdaSteps - 1;
	for (int i=0; i<nEarths; i++)
	{
		const vector3& currentEarth = Earths[i];
		const vector<double> fluxes = aster.ThermalFlux(A2Sun, currentEarth, lambdaMu, rAU,
												deltaKM, model);
		// Simpson integral over fluxes (wrt. d lambda)
		energy += (fluxes[0]+fluxes[nLambda]) * 0.5;
		double odd=0, even=0;
		{
		for (int j=1; j<nLambda; j+=2)
			odd += fluxes[j];
		}
		{
		for (int j=2; j<nLambda; j+=2)
			even += fluxes[j];
		}
		energy += (2*odd+even)*2/3.;
		cerr<<'.';cerr.flush();
	};
	energy *= 4*vector3::PI*deltaKM*deltaKM * lambdaWidth * 1e6/nEarths;

	cerr<<"Attempting to write file ";
	string nameBase = "R:/energy.inertia500.";
	//string nameBase = "R:/energy.Lambertian.";
	ostringstream outname;
	outname<<nameBase<<pV<<'.'<<Inertia
		<<'.'<<gammaDeg
		<<".dat"<<ends;
	ofstream out (outname.str().c_str());
	cerr<<outname.str()<<endl;

	//ofstream& out = cout;

	out<<"Total energy emitted: "<<energy<<" W\n";
	cerr<<"Total energy emitted  "<<energy<<" W/n";
	out<<"Should be:            ";
	cerr<<"Should be             ";

	// no craters:
	double Aeff = aster.getBond(); 

	// with craters:
	//const double S=0.5*(1-cos(gammaDeg*vector3::DEGREES*0.5)); //(1-cos(gamma/2))/2
	//double Aeff = aster.getBond() * (1-S)/(1-S*aster.getBond());

	double radius = aster.getDiameter()/2.;
	double absorbed = (1-Aeff)*AsteroidConstants::solar/(rAU*rAU) * vector3::PI*radius*radius*1e6;
	out<<absorbed<<'\n';
	cerr<<absorbed<<'\n';

	out<<endl
		<<"Parameters used were:"<<endl;
	out<<"Themal inertia + craters, inertia = "<<Inertia
		<<" (SI-units), craters' opening angle "<< gammaDeg<<endl;
	//out<<"Lambertian emitter."<<endl;

	out<<"Emissivity used: "<<emissivity<<", Bond albedo: "
		<<aster.getBond()<<endl;
	out<<"H used: "<<H<<" pV = "<<pV<<", diameter = "
		<<aster.getDiameter()<<" km"<<endl;

	out<<"Number of directions (observers) sampled: "<<nEarths<<endl;
	out<<"Sun @ 1AU, shining on equator."<<endl;
	if (!out)
		cerr<<"Writing to file somehow failed!"<<endl;
	else
		cerr<<"Writing seems to have been successful!"<<endl;


	cin.get();

	return 0;
};