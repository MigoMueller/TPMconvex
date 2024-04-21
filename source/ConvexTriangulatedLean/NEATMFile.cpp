// NEATMFile.cpp: implementation of the NEATMFile class.
//
//////////////////////////////////////////////////////////////////////

#include<sstream>
#include<iostream>
#include "NEATMFile.h"
using namespace std;


const string         NEATMFile::sphereFile = "y:/Spheres/cube6.obj";
const ObjFile::proxy NEATMFile::obj(sphereFile.c_str());

// start values for fit
const double NEATMFile::etaMin       = 0.5;  // no smaller eta can be fitted to!
const double NEATMFile::etaStep      = 0.3;  // initial step width for 'boxing' the minimum
// initial guess for pV
const double NEATMFile::pvDefault    = 0.15;
// largest possible values (fit will throw invalid_argument when they're exceeded!)
const double NEATMFile::etaMAX       = 4.5;
const double NEATMFile::pvMAX        = 1.5;

// How many refinements are made after minimum is boxed
// Final accuracy of eta is etaStep/2^nRefinements
const int    NEATMFile::nRefinements = 5;



NEATMFile::NEATMFile(const char*const fileName)
:H(0), G(0), alphaDeg(0), rAU(0), deltaAU(0),
 Sun2Aster(0,0), Earth2Aster(0,0)
{
	// all const'ness is casted away inside the ctor!
	double& H = const_cast<double&> (this->H);
	double& G = const_cast<double&> (this->G);
	double& alphaDeg = const_cast<double&> (this->alphaDeg);
	double& rAU      = const_cast<double&> (this->rAU);
	double& deltaAU  = const_cast<double&> (this->deltaAU);
	fluxMap& data    = const_cast<fluxMap&> (this->data);

	std::ifstream in(fileName);
	if (!in)
	{
		ostringstream dummy;
		dummy<<"Couldn't open "<<fileName<<" for reading!"<<ends;
		throw invalid_argument (dummy.str().c_str());
	};

	skipComments(in);
	in>>H;
	skipComments(in);
	in>>G;	
	skipComments(in);
	in>>alphaDeg;
	skipComments(in);
	in>>rAU;	
	skipComments(in);
	in>>deltaAU;
	unsigned int N;
	skipComments(in);
	in>>N;
	eclipticVector& Sun2Aster   = const_cast<eclipticVector&> (this->Sun2Aster);
	eclipticVector& Earth2Aster = const_cast<eclipticVector&> (this->Earth2Aster);
	Sun2Aster   = eclipticVector(0,0, rAU);
	Earth2Aster = eclipticVector(alphaDeg, 0, deltaAU);

	if (!in)
	{
		ostringstream dummy;
		dummy<<fileName<<" malformatted: Header must contain H, G, pV, r(AU), delta(AU) and N\n"
			 <<"where N is #points to follow!"<<ends;
		throw invalid_argument(dummy.str().c_str());
	};

	double lambdaMu, flux, error;
	while (!in.eof())
	{
		skipComments(in);
		in>>lambdaMu;
		skipComments(in);
		in>>flux;
		skipComments(in);
		in>>error;
		if (!in)
		{
			if (!in.eof())
			{
				ostringstream dummy;
				dummy<<fileName<<" is malformatted in data line "<<data.size()<<'!'<<ends;
				throw invalid_argument(dummy.str().c_str());
			};
		}
		else
		{			
			if (lambdaMu<=0 || flux<=0 || error<=0)
			{
				ostringstream dummy;
				dummy<<fileName<<" is malformatted in data line "<<data.size()<<'!'<<ends;
				throw invalid_argument(dummy.str().c_str());
			};
			//data.insert(pair<double, ValueWithError> (lambdaMu, ValueWithError(flux, error)));
			data.insert(dataLine(lambdaMu, flux, error));
		}; 	
	}; // while (!in.eof())
	
	if (data.size() != N)
	{
		cerr<<"Inconsistency in NEATM-file"<<fileName<<":"<<endl;
		cerr<<"#points given doesn't match the actual #points found later on!"<<endl;
		cerr<<"Continue using the #points found."<<endl;
	};
}; // NEATMFile::ctor



void NEATMFile::skipComments(std::ifstream &in) const
{
	while (in.peek() == '#')
	{
		char dummy[513];
		in.getline(dummy, 512);
	};
}; // NEATMFile::skipComments



NEATMFile::etaPV NEATMFile::fit(double accuracy) const
{
	const SpinState axis(0, 90, 24); // spin axis is z-axis
	ConvexTriangulated asteroid(obj, axis, pvDefault, H, G); 

	// box the minimum by letting eta run (pV is updated 'on the fly')
	// use a coarser accuracy-goal to speed up

	const double coarseAccuracy=5*accuracy;

	etaPV leftParms(etaMin, pvDefault);
	double leftChi2 = fitPV(asteroid, leftParms, coarseAccuracy);
	etaPV midParms(etaMin+etaStep, leftParms.pV);
	double midChi2 = fitPV(asteroid, midParms, coarseAccuracy);
	etaPV rightParms (etaMin+2*etaStep, midParms.pV);
	double rightChi2 = fitPV(asteroid, rightParms, coarseAccuracy);

	while (rightChi2<midChi2)
	{ // forget about the left values, proceed to higher eta
		leftParms=midParms;   leftChi2=midChi2;
		midParms =rightParms; midChi2 =rightChi2;

		rightParms.eta += etaStep; // let pV constant (1st guess)
		if (rightParms.eta > etaMAX)
		{
			ostringstream dummy;
			dummy<<"NEATMFile::fit: need an eta higher than "<<etaMAX<<", which the programmer deems unphysical!"<<ends;
			throw invalid_argument (dummy.str());
		};
		rightChi2 = fitPV(asteroid, rightParms, coarseAccuracy);
	}; // minimum boxed

	// update Chi^2 to higher accuracy (just in case)
	leftChi2  = fitPV(asteroid, leftParms,  accuracy);
	midChi2   = fitPV(asteroid, midParms,   accuracy);
	rightChi2 = fitPV(asteroid, rightParms, accuracy);

	if (leftChi2<midChi2) // may have happened due to numerical instability
	{
		rightParms=midParms;  rightChi2=midChi2;
		midParms  =leftParms; midChi2  =leftChi2;
		leftParms.eta -= etaStep;
		leftChi2      = fitPV(asteroid, leftParms, accuracy);

		if (leftChi2<midChi2) // shouldn't happen, anyway!
			int a=12;
	}
	else if (rightChi2<midChi2)
	{
		leftParms=midParms;   leftChi2=midChi2;
		midParms =rightParms; midChi2 =rightChi2;
		rightParms.eta+=etaStep;
		rightChi2     = fitPV(asteroid, rightParms, accuracy);

		if (rightChi2<midChi2)
			int a=12;
	};

	// now do bisections
	for (int i=0; i<nRefinements; i++)
	{
		if (leftChi2<rightChi2)
		{ // most likely minimum not on the 'right hand' side
			etaPV newParms = leftParms.average(midParms); // go to middle of lhs interval
			double newChi2 = fitPV(asteroid, newParms, .01);
			if (newChi2 < midChi2)
			{ // new minimum on lhs
				rightParms=midParms; rightChi2=midChi2;
				midParms  =newParms; midChi2  =newChi2;
			}
			else
			{
				etaPV verynewParms = midParms.average(rightParms);
				double verynewChi2 = fitPV(asteroid, verynewParms, .01);
				if (verynewChi2 < midChi2)
				{ // guessed wrong in the first place
					leftParms=midParms;     leftChi2=midChi2;
					midParms =verynewParms; midChi2 =verynewChi2;
				}
				else
				{ // old estimate remains best one
					leftParms =newParms;     leftChi2 =newChi2;
					rightParms=verynewParms; rightChi2=verynewChi2;
				};
			};
		}
		else
		{ // most likely minimum not on lhs
			etaPV newParms = midParms.average(rightParms); // go to middle of rhs interval
			double newChi2 = fitPV(asteroid, newParms, .01);
			if (newChi2 < midChi2)
			{ // new minimum on rhs
				leftParms=midParms;  leftChi2=midChi2;
				midParms =newParms;  midChi2 =newChi2;
			}
			else
			{
				etaPV verynewParms = midParms.average(leftParms);
				double verynewChi2 = fitPV(asteroid, verynewParms, .01);
				if (verynewChi2 < midChi2)
				{
					rightParms=midParms;     rightChi2=midChi2;
					midParms  =verynewParms; midChi2=verynewChi2;
				}
				else
				{
					rightParms=newParms;     rightChi2= newChi2;
					leftParms =verynewParms; leftChi2 = verynewChi2;
				};
			};
		};
	}; // for (i<nRefinements)
	return midParms;
}; // NEATMFile::fit



double NEATMFile::chi2(ConvexTriangulated &asteroid, double pV, double eta) const
{
	asteroid.setPV(pV);
	fluxMap::const_iterator it=data.begin(), last=data.end();
	double chi2=0, summand, lambda, modelFlux;
	while (it!=last)
	{
		lambda=it->first;
		modelFlux=asteroid.NEATMFlux(Sun2Aster, Earth2Aster, 0, lambda, eta);
		const int nOccurances = data.count(lambda);
		for (int i=0; i<nOccurances; ++i, ++it)
		{
			const ValueWithError& fluxAndSuch=it->second;
			summand = (modelFlux-fluxAndSuch.value)/fluxAndSuch.sigma;
			summand*=summand;
			chi2 += summand;
		};
	};
	return chi2;
}; // NEATMFile::chi2



// Algorithm: 
// To first order, changing pV affects model fluxes like: model ~ 1/pV
// Chi^2 = \sum_i [ (correction*model_i - flux_i)/sigma_i ]^2
// --> best fit correction =  ( \sum_i model_i*flux_i/sigma_i^2 ) / (\sum_i model_i^2/sigma_i^2)
// (by simple derivating Chi^2 w.r.t. correction and equation the derivative to 0)
double NEATMFile::fitPV(ConvexTriangulated& asteroid, double eta, double &pV, double accuracy) const
{
	double correction;
	// sum over 1. (model_i/sigma_i)^2 ; 2. model_i*data_i/(sigma_i)^2 ; 3. (data_i/sigma_i)^2
	double modelmodel, modeldata, datadata;
	do
	{
		modelmodel=modeldata=datadata=0; 
		asteroid.setPV(pV);

		fluxMap::const_iterator it=data.begin(), end=data.end();
		while (it!=end)
		{
			const double lambda=it->first;
			const double model =asteroid.NEATMFlux(Sun2Aster, Earth2Aster, 0, lambda, eta);
			const int nOccurances = data.count(lambda); // how many data points with lambda?
			for (int i=0; i<nOccurances; i++, ++it)
			{
				const double flux         =     it->second.value;
				const double sigmaSquared = pow(it->second.sigma, 2);
				modelmodel += model*model/sigmaSquared;
				modeldata  += model*flux /sigmaSquared;
				datadata   += flux*flux  /sigmaSquared;
			};
		};
		correction=modeldata/modelmodel;
		pV/=correction;
		if (pV>pvMAX)
		{
			ostringstream dummy;
			dummy<<"NEATMFile::fitPV: Encountered a geometric albedo > "<<pvMAX<<", which the developer deems unrealistic!"<<ends;
			throw invalid_argument (dummy.str());
		};
	} // do
	while ( fabs(correction-1) > accuracy );

	// apply last correction to estimate for chi^2 (don't be fooled, this is only a binomial performed!)
	return modelmodel*correction*correction - 2*correction*modeldata + datadata;
}; // NEATMFile::fitPV
