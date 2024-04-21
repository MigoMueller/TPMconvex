// Color correction 
// Author: Michael.Mueller@dlr.de
// Date:   2006 Jan 19

// Batch mode for MIPS
// Michael.Mueller@as.arizona.edu
// 2008 May 27


#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
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
  if (nargs != 5)
    {
      cerr<<"CC_MIPS, calculates MIPS color-correction factors for NEATM-like spectra."<<endl;
      cerr<<"Please provide four arguments: eta, pv, r (AU), and phase angle (deg)!"<<endl;
      return -1;
    }
  const double eta = atof(argsv[1]);
  const double pV  = atof(argsv[2]);
  const double rAU = atof(argsv[3]);
  const double alphaDeg = atof(argsv[4]);
  
  const string homeDir = getenv("HOME");
  const string bandpassFiles [] = {
    homeDir+"/__work/results/MIPS/ColorCorrection/Transmit24_divided.txt",
    homeDir+"/__work/results/MIPS/ColorCorrection/Transmit70_divided.txt",
    homeDir+"/__work/results/MIPS/ColorCorrection/Transmit160_divided.txt"
  };
  const double nominalLambdaMus [] = { 23.6748, 71.4194, 155.9};
  
  BlackBodyFtor  ref_spectrum(10000); // 10,000K black body
  const string sphereFileName = homeDir+"/tools/6144facets.obj.convex";
  NEATMftor target_spectrum(sphereFileName, eta, alphaDeg, rAU, pV);
	
  cerr<<"# MIPS color-correction factors for NEATM spectrum with following parameters:"<<endl;
  cerr<<"# [eta, pV, r (AU), phase angle (deg)]:"<<endl
      <<"# "<<eta<<'\t'<<pV<<'\t'<<rAU<<'\t'<<alphaDeg<<'\t'<<endl;
  cerr<<"# Channel\tdivide by\tmultiply by"<<endl;
  for (int i=0; i<3; i++)
    {
      const string& bandpassFile = bandpassFiles[i];
      const double nominalLambdaMu = nominalLambdaMus[i];
      vector<pair<double,double> > filter = readpass(bandpassFile);
      cout<<i+1<<"\t\t";
      double CCfactor = CC(ref_spectrum, target_spectrum, filter, nominalLambdaMu);
      cout<<CCfactor<<"\t\t"<<1./CCfactor<<endl;
    }; // end for loop
  return 0;
};
