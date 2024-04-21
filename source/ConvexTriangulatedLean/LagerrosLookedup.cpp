// LagerrosLookedup.cpp: implementation of the LagerrosLookedup class.
//
//////////////////////////////////////////////////////////////////////

#include "LagerrosLookedup.h"
#include "XFile.h"
#include <fstream>
#include <sstream>
#include <strstream>
#include <string>
#include <direct.h>
#include <cmath>
#include "vector3.h"
using namespace std;


// use the (repaired) lookup tables with an accuracy of 0.001;
// must end with '/'!
const string gammaini::directoryRoot="E:/C++/LagerrosLookup/1e-3/";


arrays::pointerArray<XFileManager*> LagerrosLookedup::managers; // trickily initialized below...
const gammaini    LagerrosLookedup::gammaList   ("gamma.ini");
const bluenessini LagerrosLookedup::bluenessList("blueness.ini");
const int LagerrosLookedup::xxxxxxx = LagerrosLookedup::initializeManagers();
// static members of LagerrosLookedup

const int LagerrosLookedup::initializeManagers()
{
	const int nGamma = gammaList.getN();
	managers.makeSize(nGamma);
	for (int i=0; i<nGamma; i++)
		managers[i] = new XFileManager(gammaList.gammas[i]);
	return 0;
};



LagerrosLookedup::LagerrosLookedup (int gamma, double blueness)
:
gamma(gamma), blue(blueness),
bluenessInfo(manager(gamma).getFiles(blueness)),
X0(bluenessInfo.Xlow),
X1(bluenessInfo.Xhigh),
blueRemainder (bluenessInfo.interpol),
sunStep   (X0.getSunStep()),
earthStep (X0.getEarthStep()),
aziStep   (X0.getAziStep())
{
	if (sunStep != X1.getSunStep())
		throw invalid_argument("");
	if (earthStep != X1.getEarthStep())
		throw invalid_argument("");
	if (aziStep != X1.getAziStep())
		throw invalid_argument("");
} // LagerrosLookedup::ctor



LagerrosLookedup::values::values(double value, double first, double last, double step)
{
	if(value<first)
	{
		if ( (first-value)/step  > 1e-3)
			throw logic_error("Value out off bounds!");
		value=first;
	};
	if (value>last)
	{
		if ( (value-last)/step > 1e-3)
			throw logic_error("Value out off bounds!");
		value=last;
	};
	remainder = (value-first)/step;
	iLow = floor(remainder);
	remainder -= iLow;
	iHigh = iLow+1;
	//if (value == last) 
	//	iHigh=iLow;
	//else
	//	iHigh=iLow+1;
	if (iHigh > (last-first)/step)
		iHigh=iLow;
};


// 1D-interpolation
inline double LagerrosLookedup::interpol(double low, double high, double remainder)
	{return low + remainder*(high-low);};
// 2D-interpolation
inline double LagerrosLookedup::interpol (double lowlow, double lowhigh,
						double highlow, double highhigh,
						double remainder1, double remainder2)
{return interpol( interpol(lowlow, lowhigh, remainder2), interpol(highlow, highhigh, remainder2), remainder1);};


double LagerrosLookedup::interpolate(double cosSun, double cosEarth, double cosAzi) const
{
	values Sun   (cosSun,   0, 1, sunStep);
	values Earth (cosEarth, 0, 1, earthStep);
	values azi   (cosAzi,  -1, 1, aziStep);

	double lowerBlue = interpol // 3D-interpolation based on X0 (lower blueness)
		( interpol(X0(Sun.iLow, Earth.iLow, azi.iLow),  X0(Sun.iLow, Earth.iLow, azi.iHigh),
				   X0(Sun.iLow, Earth.iHigh, azi.iLow), X0(Sun.iLow, Earth.iHigh, azi.iHigh),
				   Earth.remainder, azi.remainder),
		  interpol(X0(Sun.iHigh, Earth.iLow, azi.iLow),  X0(Sun.iHigh, Earth.iLow, azi.iHigh),
				   X0(Sun.iHigh, Earth.iHigh, azi.iLow), X0(Sun.iHigh, Earth.iHigh, azi.iHigh),
				   Earth.remainder, azi.remainder),
		  Sun.remainder);
	double upperBlue = interpol // 3D-interpolation at higher blueness (X1)
		( interpol(X1(Sun.iLow, Earth.iLow, azi.iLow),  X1(Sun.iLow, Earth.iLow, azi.iHigh),
				   X1(Sun.iLow, Earth.iHigh, azi.iLow), X1(Sun.iLow, Earth.iHigh, azi.iHigh),
				   Earth.remainder, azi.remainder),
		  interpol(X1(Sun.iHigh, Earth.iLow, azi.iLow),  X1(Sun.iHigh, Earth.iLow, azi.iHigh),
				   X1(Sun.iHigh, Earth.iHigh, azi.iLow), X1(Sun.iHigh, Earth.iHigh, azi.iHigh),
				   Earth.remainder, azi.remainder),
		  Sun.remainder);
	return interpol(lowerBlue, upperBlue, blueRemainder);
}



gammaini::gammaini(const string& gammaName)
{
	int result = chdir(directoryRoot.c_str());
	if (result != 0)
	{
		stringstream dummy;
		dummy<<"Couldn't change to directory "<<directoryRoot<<"!";
		throw logic_error(dummy.str());
	};

	ifstream inFile (gammaName.c_str());
	if (!inFile)
	{
		stringstream dummy;
		dummy<<"Couldn't open the initialisation file "<<gammaName<<" in "<<directoryRoot<<"!";
		throw invalid_argument(dummy.str());
	};

	stringstream in;
	in<<inFile.rdbuf();

	while (in.peek() == '#')
	{
		char dummy[513];
		in.getline(dummy, 512);
		if (!in)
			throw invalid_argument ("Initialisation file seems to be corrupted!");
	};
	
	in>>n;
	if (!in)
	{
		stringstream dummy;
		dummy<<"The initialisation file "<<gammaName<<" in "<<directoryRoot<<" seems to be corrupted!";
		throw invalid_argument(dummy.str());
	};

	gammas.makeSize(n);
	names.makeSize(n);

	for (int i=0; i<n; i++)
		in>>gammas[i]>>names[i];
	if (!in)
	{
		if (! in.eof())
		{
			stringstream dummy;
			dummy<<"The initialisation file "<<gammaName<<" in "<<directoryRoot<<" seems to be corrupted!";
			throw invalid_argument(dummy.str());
		};
	};
};  // gammaini::ctor



bluenessini::bluenessini(const string& fileName)
{
	ifstream inFile(fileName.c_str());
	if (!inFile)
	{
		stringstream dummy;
		dummy<<"Couldn't open initialisation file "<<fileName<<"!";
		throw invalid_argument(dummy.str());
	};
	stringstream in;
	in<<inFile.rdbuf(); // read in entire file at once -- should speed up!

	while (in.peek() == '#')
	{
		char dummy [513];
		in.getline(dummy, 512);
		if (!in)
		{
			stringstream dummy;
			dummy<<"Initialisation file "<<fileName<<" seems to be corrupted!";
			throw invalid_argument(dummy.str());
		};
	};

	in>>n;
	if (!in)
	{
		stringstream dummy;
		dummy<<"Initialisation file "<<fileName<<" seems to be corrupted!";
		throw invalid_argument(dummy.str());
	};

	blue.makeSize(n);
	names.makeSize(n);

	for (int i=0; i<n; i++)
		in>>blue[i]>>names[i];
	if (!in)
	{
		if (!in.eof())
		{
			stringstream dummy;
			dummy<<"Initialisation file "<<fileName<<" seems to be corrupted!";
			throw invalid_argument(dummy.str());
		}
	};

	first = blue[0];
	last  = blue[n-1];
	step  = (last-first)/(n-1);
}; // bluenessini::ctor



const string gammaini::dir(int gamma) const
{
	try
	{
		int i=0;
		while (gammas[i] != gamma)
			i++; // go through entire list -- if we get beyond, array throws an exception
		ostringstream dummy;
		dummy<<gammaini::directoryRoot<<names[i];
		return dummy.str();
	}
	catch(invalid_argument& )
	{
		stringstream dummy;
		dummy<<"There is no opening angle of "<<gamma<<" in our list!";
		throw invalid_argument(dummy.str());
	};
}; // gammaini::dir



const string& bluenessini::filename (int blue) const
{
	try
	{
		int i=0;
		while (this->blue[i] != blue)
			i++;
		return names[i];
	}
	catch (invalid_argument& )
	{
		stringstream dummy;
		dummy<<"There is no blueness of "<<blue<<" in our list!";
		throw invalid_argument(dummy.str());
	};
}; // bluenessini::filename



inline XFileManager& LagerrosLookedup::manager(int gamma) const
{
	for (int i=0; i<gammaList.n; i++)
	{
		if (gammaList.gammas[i] == gamma)
			return *managers[i];
	}
	stringstream dummy;
	dummy<<"There is no opening angle of "<<gamma<<" degrees in our list!";
	throw invalid_argument (dummy.str().c_str());
} // LagerrosLookedup::manager




double LagerrosLookedup::CorrectionFactor(const vector3 &dA, const vector3 &A2Earth, const vector3 &A2Sun)
{
	vector3 unitNormal (dA);
	unitNormal/=dA.modulus();

	double mu0= unitNormal.dot(A2Sun);
	if (mu0<=0)
		return 0;
	double mu = unitNormal.dot(A2Earth);
	if (mu<=0)
		return 0;

	double cosAzimuth;
	if (mu==1 || mu0==1)
		cosAzimuth = 0;
	else
		cosAzimuth = (A2Earth.dot(A2Sun) - mu*mu0) / sqrt( (1-mu*mu)*(1-mu0*mu0) );
	return interpolate(mu0, mu, cosAzimuth);		
} // LagerrosLookedup::CorrectionFactor
