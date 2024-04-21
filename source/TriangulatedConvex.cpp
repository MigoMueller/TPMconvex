
#include "TriangulatedConvex.h"
#include<ostream>
#include"vectorN.h"
#include"jansky.h"
using namespace std;


TriangulatedConvex::TriangulatedConvex (const ConvexFile& file, const SpinState& spinAxis,
		double H, double G, double pV, double emissivity)
		: asteroid(pV, spinAxis, H, G, emissivity),
		  dA(file.dA),
		  intrinsicDiameter(file.intrinsicDiameter)
{
	setScaleFactor();
}; // TriangulatedConvex::ctor



double TriangulatedConvex::ThermalFlux (const vector3& A2Sun, const vector3& A2Earth,
		double lambdaMu, double rAU, double deltaKM, const ThermalModelConvex& model) const
{
	vector<double> lambda(1, lambdaMu);
	vector<double> result = ThermalFlux(A2Sun, A2Earth, lambda, rAU, deltaKM, model);
	return result[0];	
}; // TriangulatedConvex::ThermalFlux



vector<double> TriangulatedConvex::ThermalFlux (const vector3& A2Sun, const vector3& A2Earth,
		const vector<double>& lambdaMu, double rAU, double deltaKM, const ThermalModelConvex& model) const
{
	const vector<double>::size_type N = lambdaMu.size();
	double factor = 2*AsteroidConstants::hc2*emissivity/pow(deltaKM,2);
	factor *= pow(scaleFactor, 2);
	vector<double> prefactor(N);
	{
		for (unsigned int i=0; i<N; i++)
			prefactor[i] = factor/pow(lambdaMu[i], 5);
	};

	vector<double> result(N,0);
	vector<vector3>::const_iterator it, end=dA.end();
	// counter: only for debugging purposes
//	int counter=0;
	for (it=dA.begin(); it!=end; it++)
	{
//		if (counter == 3131)
//			double x=3;
//		counter++;
		vector<double> summand = model.fluxModFactors(A2Sun, A2Earth, *it, lambdaMu, rAU);
		//flux += model.fluxModFactors (A2Sun, A2Earth, *it, lambdaMu, rAU);
		vector<double>::iterator result_it, sum_last = result.end(), summand_it=summand.begin();
		for (result_it = result.begin(); result_it != sum_last; result_it++, summand_it++)
		{
			if (*summand_it < 0)
				throw logic_error("Negative flux encountered!");
			if (*summand_it == *summand_it)
				*result_it += *summand_it;
			else // NaN
				double x=3; // make breakpoint here!
		};
	};
	//return result*factor;
	for (unsigned int i=0; i<N; i++)
	{
		result[i] *= prefactor[i];
/*		if (result[i] == result[i])
		{}
		else // result[i] is NaN
		{
			double x=3; // make a breakpoint here!
		};
*/
	};
	return result;
}; // TriangulatedConvex::ThermalFlux 



double TriangulatedConvex::ThermalModelConvexNoInertia::fluxModFactors 
	(const vector3& A2Sun, const vector3& A2Earth, const vector3& dA, double lambdaMu, double rAU) const
{
	static const double EPS = 1e-4;
	static const double OnePlusEps  = 1+EPS;
	static const double OneMinusEps = 1-EPS;
	
	const double area = dA.modulus();

	double cosSun   = dA.dot(A2Sun) / area;
	if (cosSun < EPS)
		return 0;
	if (cosSun >= OneMinusEps)
	{
		if (cosSun >= OnePlusEps)
			throw invalid_argument ("TriangulatedConvex::ThermalModelConvexNoInertia::fluxModFactors: A2Sun must be of length 1!");
		cosSun = OneMinusEps;
	};

	double cosEarth = dA.dot(A2Earth) / area;
	if (cosEarth < EPS)
		return 0;
	if (cosEarth >= OneMinusEps)
	{
		if (cosEarth >= OnePlusEps)
			throw invalid_argument ("TriangulatedConvex::ThermalModelConvexNoInertia::fluxModFactors: A2Earth must be of length 1!");
		cosEarth = OneMinusEps;
	};

	double cosAzi = (A2Sun.dot(A2Earth) - cosSun*cosEarth) / sqrt( (1-cosSun*cosSun) * (1-cosEarth*cosEarth));

	return area*fluxModFactorsPerArea(cosSun, cosEarth, cosAzi, lambdaMu, rAU);
}; // TriangulatedConvex::ThermalModelConvexNoInertia::fluxModFactors


void TriangulatedConvex::ThermalLightCurve(std::ostream& out, 
			const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster, 
			double lambdaMu, double JD0, const ThermalModelConvex& model, int nPointsEach) const
{
	// makes two revolutions with nPointsEach points each (only one will be calculated, however)
	// (if thermal inertia is involved: nTime points each - nPointsEach will simply be ignored)
	out<<"# Thermal light curve\n#\n";
	out<<"# Format: Time(JD)  /  flux(W/m^2/mu)  /  flux(mJy)\n";
	out<<"#\n# Plotting two entire periods at a wavelength of "<<fixed<<lambdaMu<<" microns.\n#"<<endl;

	vector3 A2Sun   = bodyFix(Sun2Aster, JD0);
	A2Sun *= -1; // unit vector pointing away from asteroid at zero time
	vector3 A2Earth = bodyFix(Earth2Aster, JD0);
	A2Earth *= -1;	
	double rAU= Sun2Aster.getLength();
	double deltaKM=Earth2Aster.getLength()*AsteroidConstants::AU;

	double factor = 2*AsteroidConstants::hc2*emissivity/pow(deltaKM,2);
	factor *= pow(scaleFactor, 2)/pow(lambdaMu, 5);

	const vector<double>::size_type N=dA.size();
	vectorN lightcurve = model.ThermalLightCurveModFactors 
		(A2Sun, A2Earth, dA[0], lambdaMu, rAU, nPointsEach);
	for (unsigned int i=1; i<N; i++)
		lightcurve += model.ThermalLightCurveModFactors
							(A2Sun, A2Earth, dA[i], lambdaMu, rAU, nPointsEach);
	lightcurve *= factor;

	const vector<double>::size_type nPoints = lightcurve.size(); // need not be nPointsEach!
	const double dT = period/(24*nPoints); // spacing between two points in units of days

	double lightcurve_average = 0;
	for (unsigned int i=0; i<nPoints; i++)
		lightcurve_average += lightcurve[i];
	lightcurve_average /= nPoints;
	out<<"# Lightcurve average: "
	   <<scientific<<lightcurve_average
	   <<fixed<<" W/m^2/mu\t"<<jansky::SI2mJy(lightcurve_average, lambdaMu)
	   <<" mJy\n#\n";

	for (unsigned int i=0; i<nPoints; i++)
	{
		out<<fixed<<JD0<<'\t'
			<<scientific<<lightcurve[i]<<'\t'
			<<fixed<<jansky::SI2mJy(lightcurve[i], lambdaMu)
			<<'\n';
		JD0 += dT;
	}; // first revolution plotted
	for (unsigned int i=0; i<nPoints; i++)
	{
		out<<fixed<<JD0<<'\t'
			<<scientific<<lightcurve[i]<<'\t'
			<<fixed<<jansky::SI2mJy(lightcurve[i], lambdaMu)
			<<'\n';
		JD0 += dT;
	}; // second revolution plotted
	out.flush();
};



void TriangulatedConvex::ThermalLightCurve(std::ostream& out, 
			const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster, 
			double lambdaMu, double JD0, 
			const ThermalModelConvex& model1, const ThermalModelConvex& model2, const double density1,
			int nPointsEach) const
{
	// makes two revolutions with nPointsEach points each (only one will be calculated, however)
	// (if thermal inertia is involved: nTime points each - nPointsEach will simply be ignored)
	out<<"# Thermal light curve\n#\n";
	out<<"# Format: Time(JD)  /  flux(W/m^2/mu)  /  flux(mJy)\n";
	out<<"#\n# Plotting two entire periods at a wavelength of "<<fixed<<lambdaMu<<" microns.\n#"<<endl;

	vector3 A2Sun   = bodyFix(Sun2Aster, JD0);
	A2Sun *= -1; // unit vector pointing away from asteroid at zero time
	vector3 A2Earth = bodyFix(Earth2Aster, JD0);
	A2Earth *= -1;	
	double rAU= Sun2Aster.getLength();
	double deltaKM=Earth2Aster.getLength()*AsteroidConstants::AU;

	double factor = 2*AsteroidConstants::hc2*emissivity/pow(deltaKM,2);
	factor *= pow(scaleFactor, 2)/pow(lambdaMu, 5);

	const vector<double>::size_type N=dA.size();
	vectorN lightcurve1 = model1.ThermalLightCurveModFactors 
		(A2Sun, A2Earth, dA[0], lambdaMu, rAU, nPointsEach);
	vectorN lightcurve2 = model2.ThermalLightCurveModFactors 
		(A2Sun, A2Earth, dA[0], lambdaMu, rAU, nPointsEach);
	for (unsigned int i=1; i<N; i++)
	{
		lightcurve1 += model1.ThermalLightCurveModFactors
							(A2Sun, A2Earth, dA[i], lambdaMu, rAU, nPointsEach);
		lightcurve2 += model2.ThermalLightCurveModFactors
			(A2Sun, A2Earth, dA[i], lambdaMu, rAU, nPointsEach);
	};

	// merge both lightcurves into one:
	lightcurve1*=density1;
	lightcurve2*=(1-density1);
	vectorN lightcurve(lightcurve1);
	lightcurve+=lightcurve2;
	lightcurve *= factor;

	const vector<double>::size_type nPoints = lightcurve.size(); // need not be nPointsEach!
	const double dT = period/(24*nPoints); // spacing between two points in units of days

	double lightcurve_average = 0;
	for (unsigned int i=0; i<nPoints; i++)
		lightcurve_average += lightcurve[i];
	lightcurve_average /= nPoints;
	out<<"# Lightcurve average: "
	   <<scientific<<lightcurve_average
	   <<fixed<<" W/m^2/mu\t"<<jansky::SI2mJy(lightcurve_average, lambdaMu)
	   <<" mJy\n#\n";

	for (unsigned int i=0; i<nPoints; i++)
	{
		out<<fixed<<JD0<<'\t'
			<<scientific<<lightcurve[i]<<'\t'
			<<fixed<<jansky::SI2mJy(lightcurve[i], lambdaMu)
			<<'\n';
		JD0 += dT;
	}; // first revolution plotted
	for (unsigned int i=0; i<nPoints; i++)
	{
		out<<fixed<<JD0<<'\t'
			<<scientific<<lightcurve[i]<<'\t'
			<<fixed<<jansky::SI2mJy(lightcurve[i], lambdaMu)
			<<'\n';
		JD0 += dT;
	}; // second revolution plotted
	out.flush();
};


vectorN TriangulatedConvex::ThermalModelConvexNoInertia::ThermalLightCurveModFactors 
		(const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
		double lambdaMu, double rAU, unsigned int nPointsEach) const
{
	vectorN result (nPointsEach);
	const double dPhi = vector3::PI2/nPointsEach;
	vector3 normal (dA);
	for (unsigned int i=0; i<nPointsEach; i++)
	{
		result[i] = this->fluxModFactors(A2Sun, A2Earth, normal, lambdaMu, rAU);
		normal.rzcw(-dPhi); // spin forward in time (asteroid moves ccw)
	};
	return result;
};
