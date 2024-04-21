#include "shape.h"
#include "jansky.h"
#include<iostream>
#include<fstream>
using namespace std;


shape::shape(double pV, 
		  const SpinState& spinAxis, // length is angular velocity in h^{-1}
		  double H, double G) // G-> 0.15 by default
		  throw(invalid_argument)
		  : spinAxis(spinAxis), H(H), G(G),
			period(2*vector3::PI/spinAxis.getLength())
{
	emissivity=0.9;
	setPV(pV);
}; // ctor(...eclipticVector...)


vector3 shape::bodyFix(const eclipticVector& vec, double JD0) const throw()
	{return bodyFix(vec.getLatitude(), vec.getLongitude(), JD0);};


// Compared my implementations of bodyfix with Mikko's (cf. email!)
// The first one is exactly what he does! (however, stated pretty differently...)
vector3 shape::bodyFix(double ecl_latitude, double ecl_longitude, double JD0) const throw()
{
	vector3 dummy(ecl_longitude - spinAxis.getLongitude(),
				 ecl_latitude); // rotate to spin axis' longitude
	dummy.rotate(2, 90.-spinAxis.getLatitude() ); // rotate to spin axis' latitude
	//dummy.rycw((90.-spinAxis.getLatitude())*vector3::DEGREES); // equivalent but less transparent
	if (JD0!=0)
		spin(dummy, 24*(JD0-spinAxis.getJD0()));
	return dummy;
}; // shape::bodyFix(double^2)
// Stick to 1st implementation to start with
/*
{
	// try new implementation
	double beta    = spinAxis.getLatitude()*vector3::DEGREES;
	double lambda  = (spinAxis.getLongitude()-ecl_longitude)*vector3::DEGREES;
	double betax   = ecl_latitude*vector3::DEGREES;
	double cosbeta   = cos(beta);
	double sinbeta   = sin(beta);
	double cosbetax  = cos(betax);
	double sinbetax  = sin(betax);
	double coslambda = cos(lambda);

	return vector3(cosbetax * sin(lambda),
		          sinbeta * cosbetax * coslambda - cosbeta * sinbetax,
				  cosbeta * cosbetax * coslambda + sinbeta * sinbetax);
}
*/
// equivalent to 1st implementation up to a phase shift of PI/2

// Q: What is standard?????????????????????????????????????????

// Difference among the two (apart from cosmetics):
// 1st: Rotate about z until lambda=0,  then rotate about y
// 2nd: Rotate about z until lambda=90, then rotate about x

// Problem: Transformation is not well defined, leaves you with a phase shift


// So, Mikko's convention is R(z'') R(y') R(z), which is pretty much Eulerian angles 
// (but also those exist in an zxz-convention...)



const unsigned int shape::fluxRecord::read(fluxRecord **& records, const char inName[])
{
	ifstream in (inName);
	if (!in)
		throw invalid_argument("fluxRecord::read: Cannot open file for reading!");

	unsigned int n;
	in>>n;
	if (!in)
		throw invalid_argument("fluxRecord::read: Couldn't read #records from file! (must be first)");
	try
	{
		records = new fluxRecord* [n];
	}
	catch (...)
	{
		cerr<<"fluxRecord::read: Ran out of memory (too many records)"<<endl;
		throw bad_alloc();
	};

	double JD, hLong, hLat, obsLong, obsLat, rAU, deltaAU, dummy, lambda, flux, sigma;
	char c;
	int i;
	for (i=0; i<n; i++)
		records[i]=0;
	try
	{
		for (i=0; i<n; i++)
		{
			in>>JD;
			in>>c; // 1st comma
			in>>c; 
			if (c!=',')
			{
				in>>c;
				if (c!=',')
					throw invalid_argument ("fluxRecord::read: Malformatted input file (three commas required after JD)");
			};
			in>>c; 
			if (c!=',')
			{
				in>>c;
				if (c!=',')
					throw invalid_argument ("fluxRecord::read: Malformatted input file (three commas required after JD)");
			};

			in>>hLong>>c; // eat one more comma
			in>>hLat>>c;
			in>>rAU>>c;
			in>>dummy>>c; // rdot, not interesting for our purposes (but given via the web-interface, nevertheless)
			in>>deltaAU>>c;
			in>>dummy>>c; // deltadot, who cares?
			in>>obsLong>>c;
			in>>obsLat>>c; // leave last comma from ephemeris -- now come the fluxes (no commas anymore!)
			in>>lambda;
			in>>flux; // in mJy
			in>>sigma;
			if (!in)
				throw invalid_argument("fluxRecord::read: Malformatted input file!");
			records[i]=new fluxRecord(JD, lambda, jansky::mJy2SI(flux, lambda), jansky::mJy2SI(sigma, lambda),
				hLong, hLat, rAU, obsLong, obsLat, deltaAU);
		}
	}
	catch (...)
	{
		for (int j=0; j<i; j++)
			delete [] records;
		delete [] records;
		// i still known!
		throw;
	};
	return n;
}; // fluxRecord::read



shape::fitFile::~fitFile()
{
	for (int i=0; i<nRecords; i++)
		delete records[i];
	delete [] records;
};