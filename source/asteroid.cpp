#include "asteroid.h"
#include "jansky.h"
#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;


const double asteroid::ThermalModel::DefaultAccuracyGoal = 0.01;


asteroid::asteroid(double pV, 
		  const SpinState& spinAxis, // length is angular velocity in h^{-1}
		  double H, double G, double emissivity) // G-> 0.15 by default
		  throw(invalid_argument)
		  : spinAxis(spinAxis), H(H), G(G),
			period(2*vector3::PI/spinAxis.getLength()),
			emissivity(emissivity)
{
	setPV(pV);
}; 



vector3 asteroid::bodyFix(const eclipticVector& vec, double JD0) const throw()
	{return bodyFix(vec.getLatitude(), vec.getLongitude(), JD0);};



// Compared my implementations of bodyfix with Mikko's (cf. email!)
// This is exactly what he does! (however, stated pretty differently...)
vector3 asteroid::bodyFix(double ecl_latitude, double ecl_longitude, double JD0) const throw()
{
	vector3 dummy(ecl_longitude - spinAxis.getLongitude(),
				 ecl_latitude); // rotate to spin axis' longitude
	dummy.rotate(2, 90.-spinAxis.getLatitude() ); // rotate to spin axis' latitude
    dummy.rotate(3, spinAxis.getPhi0()); // rotate cw by phi0
										 // this is *not* the radar definition of phi0
	                                     // (rather, 90deg must be added to it)
	                                     
	if (JD0!=0)
		spin(dummy, 24*(JD0-spinAxis.getJD0()));
	return dummy;
}; // asteroid::bodyFix(double^2)



const unsigned int asteroid::fluxRecord::readSI(fluxRecord **& records, const char inName[])
{
	ifstream in (inName);
	if (!in)
		throw invalid_argument("fluxRecord::read: Cannot open file for reading!");
	while(in.peek() == '#')
	{
		char dummy[512];
		in.getline(dummy, 511);
	};
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
	unsigned int i;
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
			in>>obsLat>>c; // eat last comma from ephemeris -- now come the fluxes (no more commas!)
			in>>lambda;
			in>>flux; // in mJy
			in>>sigma;
			if (!in)
			{
				ostringstream dummy;
				dummy<<"fluxRecord::read: Parse error in line "<<i<<"!";
				throw invalid_argument(dummy.str());
			};
			records[i]=new fluxRecord(JD, lambda, flux, sigma,
				hLong, hLat, rAU, obsLong, obsLat, deltaAU);
		}
	}
	catch (...)
	{
		for (unsigned int j=0; j<i; j++)
			delete [] records;
		delete [] records;
		// i still known!
		throw;
	};
	return n;
}; // fluxRecord::readSI



asteroid::fitFileSI::~fitFileSI()
{
	for (unsigned int i=0; i<nRecords; i++)
		delete records[i];
	delete [] records;
};



asteroid::fitFileMJy::fitFileMJy(const char inName[])
: asteroid::fitFileSI(inName)
{
	// file read in as if fluxes were in W/m^2/mu -- convert them in here
	for (unsigned int i=0; i<nRecords; i++)
	{
		fluxRecord& currentRecord = *records[i];
		currentRecord.fluxSI  = jansky::mJy2SI(currentRecord.fluxSI, currentRecord.lambdaMu);
		currentRecord.sigmaSI = jansky::mJy2SI(currentRecord.sigmaSI, currentRecord.lambdaMu);
	};
};



asteroid::ThermalModel::ThermalModel (double emissivity, double bond)
			:emissivity(emissivity), bond(bond),
			accuracyGoal(DefaultAccuracyGoal)
{
	check();
};


asteroid::ThermalModel::ThermalModel (const asteroid& shape)
			:emissivity(shape.getEmissivity()), bond(shape.getBond()),
			accuracyGoal(DefaultAccuracyGoal)
{
	check();
};



inline void asteroid::ThermalModel::check()
{
	if (emissivity < 0 || emissivity > 1)
	{
		ostringstream dummy;
		dummy<<"Emissivity must range between zero and one, encountered "<<emissivity<<"!"<<ends;
		throw invalid_argument (dummy.str());
	};
	if (bond < 0 || bond > 1)
	{
		ostringstream dummy;
		dummy<<"Bond albedo must range between zero and one, encountered "<<bond<<'!'<<ends;
		throw invalid_argument (dummy.str());
	};
}; // asteroid::ThermalModel::check



void asteroid::ThermalModel::setAccuracy(double accuracyGoal) const
{
	if (accuracyGoal <= 0)
	{
		ostringstream dummy;
		dummy<<"You attempted to set the accuracy goal to "<<accuracyGoal<<", must be nonnegative, though!"<<ends;
		throw logic_error(dummy.str());
	};
	if (accuracyGoal >= 1)
	{
		ostringstream dummy;
		dummy<<"You attempted to set the accuracy goal to "<<accuracyGoal<<", it obviously should be smaller than 1!"<<ends;
		throw logic_error(dummy.str());
	};
	this->accuracyGoal = accuracyGoal;
} // asteroid::ThermalModel::setAccuracy



double asteroid::calculateThermalParameter(double ThermalInertiaSI, double rAU) const
{
	// TSS^3
	double T3 = pow((1-bond)*AsteroidConstants::solar_over_stefan_boltzmann/(emissivity*rAU*rAU), 0.75);
	double omega = spinAxis.getOmegaH()/3600.; // in 1/s
	return ThermalInertiaSI * sqrt(omega) / (emissivity*AsteroidConstants::sigma*T3);
}; // asteroid::calculateThermalParameter



double asteroid::calculateThermalInertia(double ThermalParameter, double rAU) const
{
	// TSS^3
	double T3 = pow((1-bond)*AsteroidConstants::solar_over_stefan_boltzmann/(emissivity*rAU*rAU), 0.75);
	double omega = spinAxis.getOmegaH()/3600.; // in 1/s
	return ThermalParameter * emissivity * AsteroidConstants::sigma*T3 / sqrt(omega);
}; // asteroid::calculateThermalInertia
