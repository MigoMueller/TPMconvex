// Color correction for WISE
// Author: Michael.Mueller@oca.eu
// 2011/06/10

// This routine reproduces CC factors given in Wright et al., AJ 140, 1868-1881 (2010)
// In-band fluxes must be divided by CC factors to obtain monochromatic fluxes.
// (This is the notation used by Wright et al., even if it's stated somewhat unclearly)
//
// Caveat: ref spectrum used here is a T=10,000K black body, while Wright et al. use Vega.
// Vega is somewhat more cumbersome to include here...
// Nevertheless, agreement between CC factors calculated here and those given in Wright et al.
// is usually much better than 1%, except in channel 1 for low-temperature black bodies.
// The worst match is for the lowest temperature, T=100K, where "my" result is ~2% too large.

// Based on routines for MIPS and IRAC developed by MM since 2006/01/19


#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include"../TriangulatedConvex.h"
using namespace std;


// BAUSTELLE: WISE ref spectrum is Vega, here I only use a 
// T=10,000K black body for testing.
// Should be a fair approximation.


// reads in a filter bandpass from filename
// Assumes the file to have three columns:
// Wavelength (microns) and response/photon (*not* per flux), #3 is ignored.
// No comment lines are allowed.
vector<pair<double,double> > readpass(const string& filename)
{
  ifstream in(filename.c_str());
  if (!in)
    {
      cerr<<filename<<" appears not to exist!"<<endl;
      throw exception(); // you may want to add some sophistication here
    };
  double lambdaMu, response, dummy;
  vector<pair<double,double> > result;
  
  do
    {
      in>>lambdaMu>>response>>dummy;
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
  const string bandpassFiles [] = {
    homeDir+"/tools/wise/RSR-W1_trimmed.asc",
    homeDir+"/tools/wise/RSR-W2_trimmed.asc",
    homeDir+"/tools/wise/RSR-W3_trimmed.asc",
    homeDir+"/tools/wise/RSR-W4_trimmed.asc"
  };
  const double nominalLambdaMus [] = { 3.35, 4.60, 11.56, 22.08};
  
  //BAUSTELLE
  //  powerLawFtor  ref_spectrum(-1);
  BlackBodyFtor ref_spectrum(10000);

  const string sphereFileName = homeDir+"/tools/6144facets.obj.convex";

  const int exponents [] = {3, 2, 1, 0, -1, -2, -3, -4};
  
  cout<<"# WISE color-correction factors for exponential target spectra"<<endl;
  cout<<"# exponent, CC factors (divide by) for channels 1--4"<<endl;
  for (int j=0; j< sizeof(exponents)/sizeof(exponents[0]); j++)
    {
      powerLawFtor target_spectrum(exponents[j]);
      cout<<exponents[j]<<'\t';
      for (int i=0; i<4; i++)
	{
	  const string& bandpassFile = bandpassFiles[i];
	  const double nominalLambdaMu = nominalLambdaMus[i];
	  vector<pair<double,double> > filter = readpass(bandpassFile);
	  double CCfactor = CC(ref_spectrum, target_spectrum, filter, nominalLambdaMu);
	  cout<<CCfactor<<"\t";
	}; // end for loop
      cout<<endl;
    }

  cout<<'#'<<endl;
  cout<<"# Blackbody spectra now"<<endl;
  const int temperatures[] = {100, 141, 200, 283, 400, 566, 800, 1131, 10000};
  for (int j=0; j<sizeof(temperatures)/sizeof(temperatures[0]); j++)
    {
      BlackBodyFtor target_spectrum(temperatures[j]);
      cout<<temperatures[j]<<"K\t";
      for (int i=0; i<4; i++)
	{
	  const string& bandpassFile = bandpassFiles[i];
	  const double nominalLambdaMu = nominalLambdaMus[i];
	  vector<pair<double,double> > filter = readpass(bandpassFile);
	  double CCfactor = CC(ref_spectrum, target_spectrum, filter, nominalLambdaMu);
	  cout<<CCfactor<<"\t";
	}; // end for loop
      cout<<endl;
    }

  return 0;
};
