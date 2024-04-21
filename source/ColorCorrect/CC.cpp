// Color correction for Spitzer photometers
// Author: Michael.Mueller@dlr.de
// Date:   2006 Jan 19


#include<iostream>
#include<fstream>
#include<string>
#include"../TriangulatedConvex.h"
using namespace std;


// reads in a filter bandpass from filename
// Assumes the file to have two columns:
// Wavelength (microns) and response/photon (*not* per flux)
// The IRAC and PUI bandpasses as trimmed by MM are OK.
vector<pair<double,double> > readpass(const string& filename)
{
	ifstream in(filename.c_str());
	if (!in)
	{
		cerr<<filename<<" appears not to exist!"<<endl;
		throw exception(); // you may want to add some sophistication here
	};
	double lambdaMu, response;
	vector<pair<double,double> > result;

	do
	{
		in>>lambdaMu>>response;
		if (!in)
		{ // you might want to build in some error handling...
		  // (but check for EOF, first!)
		}
		else
		{
			//if (response < 0) response = 0; // might make sense...
			result.push_back(pair<double,double>(lambdaMu, response));
		}
	} while (in);
	return result;
};


// functor giving NEATM fluxes, for use with color correction
// eta, alpha, heliocentric distance, and geometric albedo are specified in ctor
class NEATMftor
{
protected:
	double eta, alpha, rAU, pV;
	const TriangulatedConvex        aster;
	const TriangulatedConvex::NEATM model;
	const eclipticVector Earth2Aster, Sun2Aster;
public:
	// hard wired: H=15, G=0.15, emissivity=0.9, spin axis pointing North, Sun+Earth @ equator
	NEATMftor(const string& sphereFileName, double eta, double alphaDeg, double rAU, double pV)
		: eta(eta), alpha(alphaDeg), rAU(rAU), pV(pV),
		  aster(ConvexFile(sphereFileName.c_str()), SpinState(0,90,1),15, 0.15, pV),
		  model(aster, eta),
		  Earth2Aster(0,0,1), Sun2Aster(alpha, 0, rAU)
	{};
	double operator()(double lambdaMu) const
	{
		return aster.ThermalFlux(Sun2Aster, Earth2Aster, 0, lambdaMu, model);
	};
private:
	NEATMftor();
	NEATMftor(const NEATMftor&);
	NEATMftor operator= (const NEATMftor&);
};


// spectrum proportional to f^exponent or lambda^{-exponent-2}
// (the factor lambda^{-2} is from the differential: df \propto d\lambda/\lambda^2)
class powerLawFtor
{
protected:
	const int exponent;
public:
	double operator() (double lambdaMu) const
		{return pow(lambdaMu, -exponent-2);};
	powerLawFtor(int exponent)
		: exponent(exponent) {};
};


class BlackBodyFtor
{
protected:
	const double T;
public:
	double operator() (double lambdaMu) const
		{return pow(lambdaMu, -5)/(exp(AsteroidConstants::hck/(lambdaMu*T)) - 1.);};
	BlackBodyFtor(double T)
		: T(T) {};
};


// Clumsy way of numerically integrating
// spec(lambda)*lambda*filter.second(lambda)
// over the filter bandpass specified by filter.
// Uses the trapez rule.
template<class ftor>
	double CC_integral(const ftor& spec, const vector<pair<double,double> >& filter)
{
	const vector<double>::size_type N = filter.size();
	pair<double,double> right, left(filter.front());
	double lambda=left.first;
	double right_integrand, left_integrand=spec(lambda)*lambda*left.second;
	double result=0;
	for (unsigned int i=1; i<N; i++)
	{
		right=filter.at(i);
		lambda = right.first;
		right_integrand=spec(lambda)*lambda*right.second;
		if (right_integrand < 0)
			right_integrand = 0; // discard unphysical filter transmission < 0 (given for PUI)
		result += (right_integrand+left_integrand)*(lambda-left.first);
		left=right; left_integrand=right_integrand;
	};
	result *= 0.5;
	return result;
};



// Calculates the color correction factor from reference spectrum ref
// to target spectrum target for the filter and nominal wavelength specified.
//
// For color correction, *divide* by this factor!!!
template <class refFtor, class targetFtor>
	double CC (
		   const refFtor& ref, const targetFtor& target,
		   const vector<pair<double,double> >& filter, 
		   double nominalLambdaMu)
{
	return 1./(CC_integral(ref, filter)/CC_integral(target, filter) 
		* target(nominalLambdaMu)/ref(nominalLambdaMu));
};


int main()
{
	const string bandpassFile =
	//	"c:/Documents and Settings/Migo/Desktop/SST/SST documentation/spectral response/PUI/redPU.trimmed.txt";
	////const double nominalLambdaMu = 22; // "official" wavelength in pipeline 15.3.0
	////const double nominalLambdaMu = 22.319; // effective lambda from this routine; matches that given in PU_fluxcal_memo from Jan 2008
	//const double nominalLambdaMu = 22.3; 
	//	"c:/Documents and Settings/Migo/Desktop/SST/SST documentation/spectral response/PUI/bluePU.trimmed.txt";
	////const double nominalLambdaMu = 15.6; // "official" wavelength since pipeline 15.3.0
	////const double nominalLambdaMu = 15.5;
	////const double nominalLambdaMu = 15.7971; // effective lambda from this routine; matches that given in PU_fluxcal_memo from Jan 2008
	//const double nominalLambdaMu = 15.8;
		"c:/Documents and Settings/Migo/Desktop/SST/SST documentation/spectral response/IRAC/ch4.trimmed.txt";
	const double nominalLambdaMu = 7.872;
	//	"c:/Documents and Settings/Migo/Desktop/SST/SST documentation/spectral response/IRAC/ch3.trimmed.txt";
	//const double nominalLambdaMu = 5.731;
	//	"c:/Documents and Settings/Migo/Desktop/SST/SST documentation/spectral response/IRAC/ch2.trimmed.txt";
	//const double nominalLambdaMu = 4.493;
	//	"c:/Documents and Settings/Migo/Desktop/SST/SST documentation/spectral response/IRAC/ch1.trimmed.txt";
	//const double nominalLambdaMu = 3.550;
//	    "W:/results/MIPS/ColorCorrection/Transmit24_divided.txt";
//	const double nominalLambdaMu = 23.6748;
//#define MIPS
//	    "W:/results/MIPS/ColorCorrection/Transmit70_divided.txt";
//	const double nominalLambdaMu = 71.4194;
//#define MIPS
//	    "W:/results/MIPS/ColorCorrection/Transmit160_divided.txt";
//	const double nominalLambdaMu = 155.9;
//#define MIPS

	vector<pair<double,double> > filter = readpass(bandpassFile);

#ifndef MIPS
	powerLawFtor  ref_spectrum(-1);  // for IRAC or PUI
#else
	BlackBodyFtor ref_spectrum(10000);  // for MIPS
#endif

 // Start NEATM
 //   // Centaurs:
	//const double eta = 1;
 //   const double pV = 0.1;
 //   //const double rAU = 8.7;
	//const double rAU = 15;
	////const double rAU = 17;
	////const double rAU = 26;
 //   const double alphaDeg = 2;

	//// Karins:
	//const double eta = 1;
	//const double pV = 0.2;
	//const double rAU = 3;
	//const double alphaDeg = 20;

	//// PH5:
	//const double eta = 2;
	//const double pV = 0.2;
	//const double rAU = 1;
	//const double alphaDeg = 60;

	////ML:
	//// paper as first submitted:
	////const double eta = 2.31;
	////const double pV = 0.41;
	////const double rAU = 1.27;
	////const double alphaDeg = 52.33;
	////// 2nd revision after Marco's report:
	//const double eta = 2.48;
	////const double pV = 0.35;
	//const double rAU = 1.27;
	//const double alphaDeg = 52.33;

	// Spitzer 'warm NEOs'
	// part 1: 6037
	 //those were randomly chosen
	const double pV = 0.1; 
	const double eta = 2.5;
	//// first fit:
	//const double pV = 0.3;
	//const double eta= 1.6; 
	const double alphaDeg = 46.97;
	const double rAU = 1.2424;



	//// Jovian and Saturnian satellites
	//const double eta = 1;
	//const double pV = 0.2;
	// // Jupiter
	//const double rAU=5.2;
	//const double alphaDeg=10;
	////// Saturn
	////const double rAU = 9.6;
	////const double alphaDeg = 5;

	const string sphereFileName = "s:/spheres/6144facets.obj.convex";
	NEATMftor target_spectrum(sphereFileName, eta, alphaDeg, rAU, pV);
	{
 // end NEATM



////// start power law
////	for (int i=-4;i<=4;i++)
////	{
////		powerLawFtor target_spectrum(i);
////		cout<<"alpha = "<<i<<"\t";
////// end power law



//// start black body
//	vector<double> temperatures;
// ////IRAC:	
//	//temperatures.push_back(5000);
//	//temperatures.push_back(2000);
//	//temperatures.push_back(1500);
//	//temperatures.push_back(1000);
//	//temperatures.push_back(800);
//	//temperatures.push_back(600);
//	//temperatures.push_back(400);
//	//temperatures.push_back(300);
//	//temperatures.push_back(250);
//	//temperatures.push_back(220);
//	//temperatures.push_back(200);
//	//temperatures.push_back(180);
//	//temperatures.push_back(100);
// ////PUI:
//	////temperatures.push_back(20);
//	////temperatures.push_back(25);
//	////temperatures.push_back(30);
//	////temperatures.push_back(35);
//	////temperatures.push_back(40);
//	//temperatures.push_back(50);
//	//temperatures.push_back(60);
//	//temperatures.push_back(70);
//	//temperatures.push_back(80);
//	//temperatures.push_back(100);
//	//temperatures.push_back(120);
//	//temperatures.push_back(160);
//	//temperatures.push_back(200);
//	//temperatures.push_back(240);
//	//temperatures.push_back(320);
//	//temperatures.push_back(640);
//	//temperatures.push_back(1280);
//	//temperatures.push_back(2560);
//	//temperatures.push_back(5120);
//	//temperatures.push_back(10240);
//// MIPS
//	temperatures.push_back(10000);
//	temperatures.push_back(5000);
//	temperatures.push_back(1000);
//	temperatures.push_back(500);
//	temperatures.push_back(300);
//	temperatures.push_back(200);
//	temperatures.push_back(150);
//	temperatures.push_back(100);
//	temperatures.push_back(70);
//	temperatures.push_back(50);
//	temperatures.push_back(30);
//	temperatures.push_back(20);
//
//	for (unsigned int i=0; i<temperatures.size(); i++)
//	{
//		cout<<"T = "<<temperatures[i]<<"K\t";
//		BlackBodyFtor target_spectrum(temperatures[i]);
//// end black body

		double CCfactor = CC(ref_spectrum, target_spectrum, filter, nominalLambdaMu);
		cout<<"divide by "<<CCfactor<<" or multiply by "<<1./CCfactor<<endl;
	}

	// nominal wavelength according to IRAC datahandbook 3.0 p47
	powerLawFtor oben(-1); // --> int (lambda^{-1} * lambda*R*dlambda)
	powerLawFtor unten(0); // --> int (lambda^{-2} * lambda*R*dlambda)
	cout<<"Nominal wavelength in microns: "
		<<CC_integral(oben, filter)/CC_integral(unten,filter)<<endl;

/* just testing something else
	BlackBodyFtor Sun(5800), Vega(10000);
	double V=0.545, L=3.6;
	double ratioSun = Sun(L)/Sun(V);
	double ratioVega = Vega(L)/Vega(V);
	cout<<"From V to L Sun is redder than Vega by ~"<<ratioSun/ratioVega
		<<" or ~"<<-2.5*log10(ratioSun/ratioVega)<<"mag."<<endl;
*/

	return 0;
};
