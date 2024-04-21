// Color correction 
// Author: Michael.Mueller@dlr.de
// Date:   2006 Jan 19

// Batch mode for IRAC
// Michael.Mueller@as.arizona.edu
// 2008 May 27


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
  const ConvexFile cf;
  const TriangulatedConvex        aster;
  const TriangulatedConvex::NEATM model;
  const eclipticVector Earth2Aster, Sun2Aster;
public:
  // hard wired: H=15, G=0.15, emissivity=0.9, spin axis pointing North, Sun+Earth @ equator
  NEATMftor(const string& sphereFileName, double eta, double alphaDeg, double rAU, double pV)
    : eta(eta), alpha(alphaDeg), rAU(rAU), pV(pV),
      cf(sphereFileName.c_str()),
      aster(cf, SpinState(0,90,1),15, 0.15, pV),
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


int main(int nargs, char* argsv[])
{
  const string homeDir = getenv("HOME");
  const string inputDir = homeDir + "//tools/miri/";
  const string bandpassFiles [] = {
    "MIRIM_F560W.txt",
    "MIRIM_F770W.txt",
    "MIRIM_F1000W.txt", 
    "MIRIM_F1130W.txt",
    "MIRIM_F1280W.txt",
    "MIRIM_F1500W.txt",
    "MIRIM_F1800W.txt",
    "MIRIM_F2100W.txt",
    "MIRIM_F2550W.txt",
    "MIRIM_F2550WR.txt",
    "MIRIM_F1065C.txt",
    "MIRIM_F1140C.txt",
    "MIRIM_F1550C.txt",
    "MIRIM_F2300C.txt"
  };
  const double nominalLambdaMus [] = {
    // MM 2017/09/15: waiting for updates from Maca!
    5.6,
    7.7,
    10.0,
    11.3,
    12.8,
    15.0,
    18.0,
    21.0,
    25.5,
    25.5,
    10.65,
    11.4,
    15.5,
    23.0
  };
  
  powerLawFtor  ref_spectrum(-1);
  const string sphereFileName = homeDir+"/tools/6144facets.obj.convex";

  //const int exponents [] = {-2, -1, 0, 1, 2};
  const int exponents [] = {1, 0, -1, -2, -3};
  
  cout<<"# MIRI color-correction factors for exponential target spectra"<<endl;
  cout<<"# exponent, CC factors (divide by) for filters"<<endl;
  cout<<"# ";
  for (int j=0; j<sizeof(nominalLambdaMus)/sizeof(nominalLambdaMus[0]); j++)
    cout<<bandpassFiles[j]<<'\t';
  cout<<'\n';

  for (int j=0; j< sizeof(exponents)/sizeof(exponents[0]); j++)
    {
      powerLawFtor target_spectrum(exponents[j]);
      cout<<exponents[j]<<'\t';
      for (int i=0; i<sizeof(nominalLambdaMus)/sizeof(nominalLambdaMus[0]); i++)
	{
	  const string& bandpassFile = inputDir+bandpassFiles[i];
	  const double nominalLambdaMu = nominalLambdaMus[i];
	  vector<pair<double,double> > filter = readpass(bandpassFile);
	  double CCfactor = CC(ref_spectrum, target_spectrum, filter, nominalLambdaMu);
	  cout<<CCfactor<<"\t";
	}; // end for loop
      cout<<endl;
    }

  cout<<'#'<<endl;
  cout<<"# Blackbody spectra now"<<endl;
  //const int temperatures[] = {5000, 2000, 1500, 1000, 800, 600, 400, 200};
  const int temperatures[] = {100, 200, 300, 400, 500, 600, 700, 800};
  for (int j=0; j<sizeof(temperatures)/sizeof(temperatures[0]); j++)
    {
      BlackBodyFtor target_spectrum(temperatures[j]);
      cout<<temperatures[j]<<"K\t";
      for (int i=0; i<sizeof(nominalLambdaMus)/sizeof(nominalLambdaMus[0]); i++)
	{
	  const string& bandpassFile = inputDir+bandpassFiles[i];
	  const double nominalLambdaMu = nominalLambdaMus[i];
	  vector<pair<double,double> > filter = readpass(bandpassFile);
	  double CCfactor = CC(ref_spectrum, target_spectrum, filter, nominalLambdaMu);
	  cout<<CCfactor<<"\t";
	}; // end for loop
      cout<<endl;
    }

  return 0;
};
