#include "TriangulatedConvex.h"
#include "HemisphericNoInertia.h"
#include "ThermalInertiaOnlyConvex.h"
#include "InertiaTimesCrater.h"
#include "InertiaTimesCraterLagerros.h"
#include<iostream>
#include<fstream>
#include"jansky.h"
using namespace std;


int main()
{
	//ostream& out = cout;
	//ofstream out("E:/C++/crater+inertia.alpha+20.cube6.dat");
	//ofstream out("2002CD.apr4.1015.pv02.poleon.dat");
	//ofstream out("moderateMBA.phase0.dat");
	ofstream out ("P:/sphere.alpha0.thermalinertia500.HIGHaccuracy.dat");

	const char shape_file_name [] = "S:/spheres/cube8.obj.convex";
	const ConvexFile shape (shape_file_name);

	const double ecl_lambda = 0;
	const double ecl_beta   = 90;
	const double periodH    = 6;
	const SpinState axis (ecl_lambda, ecl_beta, periodH);

	//const double H = 20.2;
	const double H = 15.8599; // -- yields D=2km
	const double G = 0.15;
	const double pV= 0.2;
	//const double emissivity = 0.9;
	const double emissivity = 1;

	TriangulatedConvex aster (shape, axis, H, G, pV, emissivity);

	//const double gammaDeg = 180;
	const double gammaDeg = 0;

	//const double rAU = 1.1115846941;;
	//const double deltaAU = 0.1120377794;
	const double deltaAU = 1, rAU = 1;
	//const double phaseDeg = +6.3;
	const double phaseDeg = 0;

	const eclipticVector Sun2Aster   (0,0,rAU);
	const eclipticVector Earth2Aster (phaseDeg, 0, deltaAU);
	
	vector<double> lambdaMu;
	lambdaMu.push_back(5);
	lambdaMu.push_back(8);
	lambdaMu.push_back(10);
	lambdaMu.push_back(12);
	lambdaMu.push_back(20);
	//lambdaMu.push_back(8.8);
	//lambdaMu.push_back(9.7);
	//lambdaMu.push_back(10.5);
	//lambdaMu.push_back(11.7);
	//lambdaMu.push_back(12.4);
	//lambdaMu.push_back(17.7);
	//lambdaMu.push_back(18.8);

	//const double inertia[] = {50, 500, 2500};
	const double inertia = 500;

	//for (double ThermalInertia = 0; ThermalInertia <= 2500; ThermalInertia += 100)
	for (int i=0; i<1; i++)
	{
		const double ThermalInertia = inertia;
		//const double ThermalInertia = inertia[i];
		const int nTime = 1300;
		const int nZ    = 65;
		const double zMax = 16;
		const double accuracyGoal = 0.00001;
		
		ThermalInertiaOnlyConvex inertia (emissivity, aster.getBond(), aster.calculateThermalParameter(ThermalInertia, rAU), 
			nTime, nZ, zMax, accuracyGoal);
		//InertiaTimesCrater inertia_crater (emissivity, aster.getBond(), aster.calculateThermalParameter(ThermalInertia, rAU), 
		//	nTime, nZ, zMax, accuracyGoal, gammaDeg);
		InertiaTimesCraterLagerros inertia_crater_2 (emissivity, aster.getBond(), aster.calculateThermalParameter(ThermalInertia, rAU), 
			nTime, nZ, zMax, accuracyGoal, gammaDeg);

		vector<double> wo_crater = aster.ThermalFlux(Sun2Aster, Earth2Aster, 0, lambdaMu, inertia);
		//vector<double> w__crater = aster.ThermalFlux(Sun2Aster, Earth2Aster, 0, lambdaMu, inertia_crater);
		vector<double> w2_crater = aster.ThermalFlux(Sun2Aster, Earth2Aster, 0, lambdaMu, inertia_crater_2);

		out<<"Thermal inertia = "<<ThermalInertia<<endl;
		out<<"lambda(mu)\tno craters\t(old craters)\t(new craters = Lagerros)"<<endl;
		for (int i=0; i<lambdaMu.size(); i++)
			//out<<lambdaMu[i]<<'\t'<<wo_crater[i]<<'\t'<<w__crater[i]<<endl;
			out<<lambdaMu[i]<<'\t'<<jansky::SI2mJy(wo_crater[i], lambdaMu[i])
		//	<<'\t'<<jansky::SI2mJy(w__crater[i], lambdaMu[i])
			<<'\t'<<jansky::SI2mJy(w2_crater[i], lambdaMu[i])<<endl;
		out<<endl;
	};
	out<<"Parameters were:"<<endl
		<<"H = "<<H<<", G = "<<G<<", pV = "<<pV<<endl
		<<"period = "<<periodH<<'h'<<endl
		<<"craters' opening angle = "<<gammaDeg<<" degrees (full angle, MM notation)"<<endl
		<<"heliocentric distance = "<<rAU<<" AU"<<endl
		<<"geocentric   distance = "<<deltaAU<<" AU"<<endl
		<<"phase angle  = "<<phaseDeg<<" degrees"<<endl
		<<"emissivity = "<<emissivity<<endl;

	cin.get();

	return 0;
};
