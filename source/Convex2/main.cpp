#include "../TriangulatedConvex.h"
//#include "../HemisphericNoInertia.h"
//#include "../ThermalInertiaOnlyConvex.h"
#include "../InertiaTimesCraterLagerros.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include <list>
#include"../jansky.h"
using namespace std;


// procedure to output a file containing chi^2 
// as a function of crater density and albedo

// ready to be plotted by 'gnuplot'
// file format: pV \t rho \t chi^2

// crater density runs from 0 thru 1, using nRho points,
// pV runs from given value (given in aster) / pV_factor thru given value * pV_factor, using nPV points

void plot_chi2(const TriangulatedConvex& aster,
			   const list<ThermalInertiaOnlyConvex*>& no_craters,
			   const list<InertiaTimesCraterLagerros*>& craters,
			   const list<asteroid::fitFileSI*>& data,
			   ostream& out,
			   unsigned int nPV = 100,
			   unsigned int nRho = 101,
			   double pV_factor = 4,
			   ostream& err = cerr)
{
	double cn, n, m; // cn: craters-nocraters, n: craters, m: measured -- all /= sigma
	double cncn=0, ncn=0, mcn=0, nn=0, mn=0, mm=0; // partial sums

	list<asteroid::fitFileSI*>::const_iterator it_fitfiles = data.begin();
	list<InertiaTimesCraterLagerros*>::const_iterator it_craters = craters.begin();
	list<ThermalInertiaOnlyConvex*>::const_iterator it_no_craters = no_craters.begin();

	for ( ; it_fitfiles != data.end(); it_fitfiles++, it_craters++, it_no_craters++)
	{
		const asteroid::fitFileSI&        data       = **it_fitfiles;
		const InertiaTimesCraterLagerros& craters    = **it_craters;
		const ThermalInertiaOnlyConvex&   no_craters = **it_no_craters;

		for (int i=0; i<data.nRecords; i++)
		{
			const asteroid::fluxRecord& record  = data[i];
			m  = record.fluxSI/record.sigmaSI;
			double nnnn = aster.ThermalFlux(record, no_craters);
			err<<'.';
			n  =  nnnn / record.sigmaSI;
			double cccc = aster.ThermalFlux(record, craters);
			err<<'.';
			cn =  (cccc / record.sigmaSI) - n;
			cncn += cn*cn;
			ncn  += n*cn;
			mcn  += m*cn;
			nn   += n*n;
			mn   += m*n;
			mm   += m*m;
		}; // all sums performed (1 data file)
	}; // all sums performed (all data files)

	const double pV = aster.getPV();
	const double pV_min = pV/pV_factor;
	const double pV_max = pV*pV_factor;
	const double pV_step = (pV_max-pV_min)/(nPV-1);

	const double rho_step = 1./(nRho-1);
	vector<double> rho (nRho);
	for (int i=0; i<nRho; i++)
		rho[i] = rho_step*i;
	rho[nRho-1] = 1;
	vector<double> rho2 (nRho);
	for (int i=0; i<nRho; i++)
		rho2[i] = rho[i]*rho[i];

	double currentPV = pV_min;
	for (int i=0; i<nPV; i++)
	{
		const double kappa = pV / currentPV;
		const double kappa2 =  kappa*kappa;

		for (int j=0; j<nRho; j++)
		{
			out<<currentPV<<'\t'
				<<rho[j]<<'\t'
				<<kappa2 * rho2[j] * cncn
				    + 2*kappa2*rho[j]*ncn
					+ kappa2 * nn
					- 2 * kappa * rho[j] * mcn
					- 2 * kappa * mn
					+ mm
				<<'\n';
		};
		currentPV += pV_step;
	};
};




// returns the best-fit chi^2
// store furthermore best-fit pV in newPV + best fit crater density in newDensity

// Method: linear regression
// minimizes Chi^2=\sum_k [ (modelflux_k-measuredflux_k)/sigma_k ]^2,
// where modelflux = pv(used)/pv(bestfit) * [ (1-rho)*nocraters + rho*craters ]
// assuming flux ~ 1/pV (pretty good approximation, handle with care, nevertheless!)

double fit(double& newPV, double& sigmaPV, 
		   double& newDensity, double& sigmaDensity,
		   double& correlation,
		   const TriangulatedConvex& aster,
		   const list<ThermalInertiaOnlyConvex*>& no_craters,
		   const list<InertiaTimesCraterLagerros*>& craters,
		   const list<asteroid::fitFileSI*>& data,
		   ostream& important_output,
		   ostream& diagnostic_output = cerr)
{
	ostream& out = important_output;
	ostream& err = diagnostic_output;

	//	out<<"# Fitting crater density and albedo, calculating "<<2*data.nRecords<<" model fluxes.\n";
	out<<"# JD\tlambda\tThermal inertia only\tThermal inertia + craters\tmeasured\tsigma"<<endl;

	double cn, n, m; // cn: craters-nocraters, n: craters, m: measured -- all /= sigma
	double cncn=0, ncn=0, mcn=0, nn=0, mn=0, mm=0; // partial sums

	list<asteroid::fitFileSI*>::const_iterator it_fitfiles = data.begin();
	list<InertiaTimesCraterLagerros*>::const_iterator it_craters = craters.begin();
	list<ThermalInertiaOnlyConvex*>::const_iterator it_no_craters = no_craters.begin();

	for ( ; it_fitfiles != data.end(); it_fitfiles++, it_craters++, it_no_craters++)
	{
		const asteroid::fitFileSI&        data       = **it_fitfiles;
		const InertiaTimesCraterLagerros& craters    = **it_craters;
		const ThermalInertiaOnlyConvex&   no_craters = **it_no_craters;

		for (int i=0; i<data.nRecords; i++)
		{
			const asteroid::fluxRecord& record  = data[i];
			out<<fixed;
			out<<"#\t"<<record.JD<<'\t'<<record.lambdaMu<<'\t';
			out.flush();

			m  = record.fluxSI/record.sigmaSI;

			double nnnn = aster.ThermalFlux(record, no_craters);
			out<<scientific<<nnnn<<'\t'; out.flush(); err<<'.';
			n  =  nnnn / record.sigmaSI;

			double cccc = aster.ThermalFlux(record, craters);
			out<<cccc
				<<'\t'
				<<record.fluxSI<<'\t'
				<<record.sigmaSI<<endl; 
			err<<'.';
			cn =  (cccc / record.sigmaSI) - n;

			cncn += cn*cn;
			ncn  += n*cn;
			mcn  += m*cn;
			nn   += n*n;
			mn   += m*n;
			mm   += m*m;
		}; // all sums performed (1 data file)
	}; // all sums performed (all data files)
	out<<'#'<<endl;
	err<<endl;

	// this solution minimizes Chi^2 (partial derivatives vanish)
	double kappa = (ncn*mcn - mn*cncn) / (ncn*ncn - nn*cncn);
	out<<"#Note: The best fit is achieved when all the model fluxes are multiplied by "<<kappa<<endl
		<<'#'<<endl;

	double sigma_kappa = sqrt( cncn/(nn*cncn - ncn*ncn) );
	double rel_kappa   = sigma_kappa / kappa;
	
	//double delta_rho = (-NCN^2+NN*CNCN)*sqrt(MCN^2*NN+MN^2*CNCN-2*NCN*MN*MCN)/((-NCN*MCN+MN*CNCN)^2);

	newPV        = aster.getPV() / kappa;
	sigmaPV      = newPV*rel_kappa;
	newDensity   = (mn*ncn - nn*mcn) / (ncn*mcn - mn*cncn);
	sigmaDensity = (nn*cncn-ncn*ncn)*sqrt(mcn*mcn*nn + mn*mn*cncn - 2*ncn*mn*mcn)/((mn*cncn-ncn*mcn)*(mn*cncn-ncn*mcn));

	correlation  = -mcn*(-ncn*ncn+nn*cncn)/((ncn*mcn-mn*cncn)*(ncn*mcn-mn*cncn)*sigma_kappa*sigmaDensity);

	if (newDensity >= 0 && newDensity <= 1)
		return newDensity*newDensity * kappa*kappa * cncn
		+ 2 * newDensity * kappa*kappa * ncn
		+ kappa*kappa * nn
		+ mm
		- 2 * kappa * mn
		- 2 * newDensity * kappa * mcn;

	// else: rho < 0 or rho > 1, check which one's 'better':
	double kappa0 = mn/nn;
	double newPV0 = aster.getPV() / kappa;
	double chi2_0 = kappa0*kappa0*nn + mm - 2*kappa0*mn;

	double kappa1 = (mcn+mn)/(nn+2*ncn+cncn);
	double newPV1 = aster.getPV() / kappa1;
	double chi2_1 = kappa1*kappa1*cncn 
			+ 2*kappa1*kappa1*ncn
			+ kappa1*kappa1*nn 
			+ mm 
			- 2*kappa1*mn
			- 2*kappa1*mcn;

	if (chi2_0 < chi2_1)
	{
		kappa = kappa0;
		newPV = newPV0;
		newDensity = 0;
		return chi2_0;
	}
	else
	{
		kappa = kappa1;
		newPV = newPV1;
		newDensity = 1;
		return chi2_1;
	};
};
/*
double fit(double& newPV, double& sigmaPV, 
		   double& newDensity, double& sigmaDensity,
		   const TriangulatedConvex& aster,
		   const ThermalInertiaOnlyConvex& no_craters,
		   const InertiaTimesCraterLagerros& craters,
		   const asteroid::fitFileSI& data,
		   ostream& important_output,
		   ostream& diagnostic_output = cerr)
{
	ostream& out = important_output;
	ostream& err = diagnostic_output;

	out<<"# Fitting crater density and albedo, calculating "<<2*data.nRecords<<" model fluxes.\n";
	out<<"\t# JD\tlambda\tThermal inertia only\tThermal inertia + craters\tmeasured\tsigma"<<endl;
	
	double cn, n, m; // cn: craters-nocraters, n: craters, m: measured -- all /= sigma
	double cncn=0, ncn=0, mcn=0, nn=0, mn=0, mm=0; // partial sums
	for (int i=0; i<data.nRecords; i++)
	{
		const asteroid::fluxRecord& record  = data[i];
		out<<fixed;
		out<<"#\t"<<record.JD<<'\t'<<record.lambdaMu<<'\t';
		out.flush();
		
		m  = record.fluxSI/record.sigmaSI;
		
		double nnnn = aster.ThermalFlux(record, no_craters);
		out<<scientific<<nnnn<<'\t'; out.flush(); err<<'.';
		n  =  nnnn / record.sigmaSI;
		
		double cccc = aster.ThermalFlux(record, craters);
		out<<cccc
		   <<'\t'
		   <<record.fluxSI<<'\t'
		   <<record.sigmaSI<<endl; 
		err<<'.';
		cn =  (cccc / record.sigmaSI) - n;
		
		cncn += cn*cn;
		ncn  += n*cn;
		mcn  += m*cn;
		nn   += n*n;
		mn   += m*n;
		mm   += m*m;
	}; // all sums performed
	
	out<<endl;
	err<<endl;

	// this solution minimizes Chi^2 (partial derivatives vanish)
	double kappa = (ncn*mcn - mn*cncn) / (ncn*ncn - nn*cncn);
	//newPV      = aster.getPV() * (ncn*ncn - nn*cncn) / (ncn*mcn - mn*cncn);
	out<<"#Note: The best fit is achieved when all the model fluxes are multiplied by "<<kappa<<endl<<endl;
	newPV        = aster.getPV() / kappa;
	newDensity = (mn*ncn - nn*mcn) / (ncn*mcn - mn*cncn);

	if (newDensity >= 0 && newDensity <= 1)
		return newDensity*newDensity * kappa*kappa * cncn
		+ 2 * newDensity * kappa*kappa * ncn
		+ kappa*kappa * nn
		+ mm
		- 2 * kappa * mn
		- 2 * newDensity * kappa * mcn;
	// rho < 0 or rho > 1, check which one's 'better':

	double kappa0 = mn/nn;
	double newPV0 = aster.getPV() / kappa;
	double chi2_0 = kappa0*kappa0*nn + mm - 2*kappa0*mn;

	double kappa1 = (mcn+mn)/(nn+2*ncn+cncn);
	double newPV1 = aster.getPV() / kappa1;
	double chi2_1 = kappa1*kappa1*cncn 
			+ 2*kappa1*kappa1*ncn
			+ kappa1*kappa1*nn 
			+ mm 
			- 2*kappa1*mn
			- 2*kappa1*mcn;
	if (chi2_0 < chi2_1)
	{
		kappa = kappa0;
		newPV = newPV0;
		return chi2_0;
	}
	else
	{
		kappa = kappa1;
		newPV = newPV1;
		return chi2_1;
	};

	if (newDensity < 0) // crop density to >= 0 / output negative value nevertheless
	{
		// best-fit kappa for density==0
		kappa = mn/nn;
		newPV = aster.getPV() / kappa;
		return kappa*kappa*nn + mm - 2*kappa*mn;
	};
	if (newDensity > 1) // crop density to <= 1 / output high value nevertheless
	{
		kappa = (mcn+mn)/(nn+2*ncn+cncn);
		newPV = aster.getPV() / kappa;
		return kappa*kappa*cncn 
			+ 2*kappa*kappa*ncn
			+ kappa*kappa*nn 
			+ mm 
			- 2*kappa*mn
			- 2*kappa*mcn;
	};
	return newDensity*newDensity * kappa*kappa * cncn
		  + 2 * newDensity * kappa*kappa * ncn
		  + kappa*kappa * nn
		  + mm
		  - 2 * kappa * mn
		  - 2 * newDensity * kappa * mcn;
};
*/



// returns the best-fit chi^2
// stores best-fit crater density in newDensity
// keeps albedo constant
// Method and terminology: as above
double fit(double& newDensity, double& sigmaDensity, const TriangulatedConvex& aster,
		 const ThermalInertiaOnlyConvex& no_craters,
		 const InertiaTimesCraterLagerros& craters,
		 const asteroid::fitFileSI& data,
		 ostream& important_output,
		 ostream& diagnostic_output = cerr) 
{
	ostream& out = important_output;
	ostream& err = diagnostic_output;
	out<<"# Fitting crater density with fixed albedo, calculating "<<2*data.nRecords<<" model fluxes.\n";
	out<<"#\tJD\tlambda\tThermal inertia only\tThermal inertia + craters\tmeasured fluxes\tsigma"<<endl;
	
	double cn, nm; // cn: craters-nocraters, nm: nocraters-measured -- all /= sigma
	double c,n,m; // dummy variables
	double cncn=0, cnnm=0, nmnm=0; // partial sums
	for (int i=0; i<data.nRecords; i++)
	{
		const asteroid::fluxRecord& record  = data[i];
		out<<"#\t"<<fixed<<record.JD<<'\t'<<record.lambdaMu<<'\t'; out.flush();

		m  = record.fluxSI;// / record.sigmaSI;

		n = aster.ThermalFlux(record, no_craters);
		out<<scientific<<n<<'\t'; out.flush(); err<<'.';
		
		c = aster.ThermalFlux(record, craters);
		out<<c<<'\t'
			<<m<<'\t'
			<<record.sigmaSI<<endl; 
		err<<'.';
		
		cn = (c-n) / record.sigmaSI;
		nm = (n-m) / record.sigmaSI;
		
		cncn += cn*cn;
		cnnm += cn*nm;
		nmnm += nm*nm;
	}; // all sums performed
	
	out<<endl;
	err<<endl;
	// this solution minimizes Chi^2 (partial derivative vanishes)
	newDensity   = -cnnm / cncn;
	sigmaDensity = sqrt(2.3 / cncn);

	if (newDensity >= 0 && newDensity <= 1)
		return newDensity*cnnm + nmnm;
	double chi2_rho0 = nmnm;                 // Chi^2 for rho = 0
	double chi2_rho1 = cncn + 2*cnnm + nmnm; // Chi^2 for rho = 1
	if (chi2_rho0 < chi2_rho1)
		return chi2_rho0;
	else
		return chi2_rho1;
/*
	if (newDensity<0)
		return nmnm; // crop density to >= 0 but output neg. value nevertheless (diagnostic)
	if (newDensity<=1)
		return newDensity*cnnm + nmnm;
	else // crop density to <= 1 but output high value nevertheless (diagnostic)
		return cncn + 2*cnnm + nmnm;
*/
};


// returns the best-fit chi^2
// stores best-fit geometric albedo in newPV
// keeps cratering constant (only one model as input)
// Method and terminology: as above
template<class model_type>
double fit(double& newPV, double& sigmaPV, const TriangulatedConvex& aster,
		   const list<model_type*>& models,
		   const list<asteroid::fitFileSI*>& data,
		   ostream& important_output,
		   ostream& diagnostic_output = cerr) 
{
	ostream& out = important_output;
	ostream& err = diagnostic_output;
	out<<"# Fitting albedo with fixed cratering\n";
	out<<"#\tJD\tlambda\tModel fluxes\tmeasured fluxes\tsigma"<<endl;
	
	double sm=0, ss=0, mm=0; //partial sums 
	                         // ss: synthetic*measured, ss: synthetic^2, mm: measured^2 
	                         // all sums /= sigma
	double s,m; // dummy variables
	
	list<asteroid::fitFileSI*>::const_iterator it_fitfiles = data.begin();
	list<model_type*>::const_iterator it_models = models.begin();
	
	for ( ; it_fitfiles != data.end(); it_fitfiles++, it_models++)
	{
		const asteroid::fitFileSI&                    data    = **it_fitfiles;
		const TriangulatedConvex::ThermalModelConvex& model   = **it_models;

		for (int i=0; i<data.nRecords; i++)
		{
			const asteroid::fluxRecord& record  = data[i];
			out<<"#\t"<<fixed<<record.JD<<'\t'<<record.lambdaMu<<'\t'; out.flush();

			m = record.fluxSI;// / record.sigmaSI;
			s = aster.ThermalFlux(record, model);
			out<<scientific;
			out<<s<<'\t'
				<<m<<'\t'
				<<record.sigmaSI<<endl; 
			err<<'.';

			s/=record.sigmaSI;
			m/=record.sigmaSI;

			sm += s*m;
			mm += m*m;
			ss += s*s;
		};
	}; // all sums performed
	err<<endl;
	// this solution minimizes Chi^2 (partial derivative vanishes)
	double kappa = sm/ss;
	out<<"#Note: The best fit is achieved when all the model fluxes are multiplied by "<<kappa<<endl
		<<'#'<<endl;

	newPV   = aster.getPV()/kappa;
	sigmaPV = newPV/(sqrt(ss)*kappa);
	return kappa*kappa*ss - 2*kappa*sm + mm;
};


// Thermal-inertia-fitting
int main(int nargs, char* argsv[])
{
	
	if (nargs != 2)
	{
		cerr<<"Specify output file!"<<endl;
		return -1;
	}
	ofstream out(argsv[1]);
	if (!out)
	{
		cerr<<"Couldn't open "<<argsv[1]<<" for writing, bail out!"<<endl;
		return -1;
	};
	
	//ostream& out = cout; // for debugging / profiling
	
try
{
	list<asteroid::fitFileSI*> fitfiles;

/*
// 1998 WT24
	fitfiles.push_back(new asteroid::fitFileMJy ("R:/1998_WT24/TPM/ESO20011204.dat"));
	//fitfiles.push_back(new asteroid::fitFileMJy ("R:/1998_WT24/TPM/IRTF20011218.dat"));
	//fitfiles.push_back(new asteroid::fitFileMJy ("R:/1998_WT24/TPM/IRTF20011218_Harris2006.dat"));
	fitfiles.push_back(new asteroid::fitFileMJy ("R:/1998_WT24/TPM/IRTF20011219.dat"));
	fitfiles.push_back(new asteroid::fitFileMJy ("R:/1998_WT24/TPM/IRTF20011221.dat"));
	//ConvexFile shape ("R:/1998_WT24/shapes/elli2_1.25.obj.convex");
	ConvexFile shape ("S:/spheres/cube5.obj.convex");

	// to run on epsi (tweak paths)
	//fitfiles.push_back(new asteroid::fitFileMJy ("fitfiles/ESO20011204.dat"));
	//fitfiles.push_back(new asteroid::fitFileMJy ("R:/1998_WT24/TPM/IRTF20011218.dat"));
	//fitfiles.push_back(new asteroid::fitFileMJy ("fitfiles/IRTF20011218_Harris2006.dat"));
	//fitfiles.push_back(new asteroid::fitFileMJy ("fitfiles/IRTF20011219.dat"));
	//fitfiles.push_back(new asteroid::fitFileMJy ("fitfiles/IRTF20011221.dat"));
	//ConvexFile shape ("shape/cube5.obj.convex");

	// spin axis:
	// assumed to be normal to viewing plane on 2001 Dec 16, 18.00 UT:
	// prograde: (170, 45) degrees J2000 ecliptic longitude/latitude
	// retro:    (350, -45)
	// best fit JD0 (manually determined, subject to optimization):
	// prograde: 2452262.03
	// retro:    2452262.016

	// prograde:
    //double ecl_lambda = 170;
	//double ecl_beta   = 45;
	//double JD0        = 2452262.03;

	// retrograde:
	//double ecl_lambda = 170+180;
	//double ecl_beta   = -45;
	//double JD0        = 2452262.016;

	// Nominal axis, Harris 2006 (or 2007)
	double ecl_lambda = 355;
	double ecl_beta   = -52;
	double JD0 = 0; // irrelevant since we're dealing with a sphere

	// Nominal axis tilted by 30 degrees towards the Sun;
	// cf. the maple worksheet in folder R:\1998_WT24\fits\sphere, new pole + 30deg off
	//double ecl_lambda = 38.16069179;
	//double ecl_beta   = -43.03439951;
	//double JD0 = 0; // see above


	// old, wrong, spin axis, used in June 2005:
    // (result from bug in ecl...)
	// prograde: (4.425, 49.3256) degrees J2000 ecliptic longitude/latitude
	// retro:    (184.425, -49.3256)
	// best fit JD0 (manually determined, subject to optimization):
	// prograde: 2452261.993
	// retro:    2452261.963

	// prograde:
	//double ecl_lambda = 4.425;
	//double ecl_beta   = 49.3256;
	//double JD0        = 2452261.993;

	// retrograde:
	//double ecl_lambda = 184.425;
	//double ecl_beta   = -49.3256;
	//double JD0        = 2452261.963;


	// first try:
	// prograde: perpendicular to viewing geometry of 2001 Dec 21
	// JD0 was (loosely) fitted to optical lightcurves by P. Pravec
	//double JD0        = 2452262.0;
	//double ecl_lambda = 1.468;
	//double ecl_beta   = 60.63;
	
	// retrograde: (as above)
	// again, JD0 was (loosely = manually) fitted to lightcurve data by PP
	//double JD0          = 2452261.97;
	//double ecl_lambda   = 181.468;
	//double ecl_beta     = -60.63;

//	double period_h   = 3.6996; // by P. Pravec
    double period_h   = 3.697; // Krugly et al. 2002

    double H = 18.5;  // Harris et al. 2006
    double G = 0.15; // maybe try 0.4???

	//double H          = 17.9; // EARN
//	double H          = 18.39; // Kiselev et al., ACM 2002
//	double G          = 0.15;  // Kiselev et al. give a slope in mag/deg --> G???
	double emissivity = 0.9;
	double pV         = 0.6;
    const double phi0 = 0;
*/  // WT24



// (5381) Sekhmet
/*
	fitfiles.push_back(new asteroid::fitFileSI("R:/sekhmet/tpm_fitfiles/sekhmet.20030512.dat"));
	fitfiles.push_back(new asteroid::fitFileSI("R:/sekhmet/tpm_fitfiles/sekhmet.20030513.dat"));
	fitfiles.push_back(new asteroid::fitFileSI("R:/sekhmet/tpm_fitfiles/sekhmet.20030514.dat"));
	fitfiles.push_back(new asteroid::fitFileSI("R:/sekhmet/tpm_fitfiles/sekhmet.20030515.dat"));
	fitfiles.push_back(new asteroid::fitFileSI("R:/sekhmet/tpm_fitfiles/sekhmet.20030516.dat"));

	ConvexFile shape ("S:/spheres/cube5.obj.convex");
	// prograde: equatorial view on 2003 May 14, 14 UT (radar?)
	double ecl_lambda = 12.5218;
	double ecl_beta   = 69.6029;

	// retrograde:
	//double ecl_lambda = 192.5218;
	//double ecl_beta   = -69.6029;

	double period_h   = 2.7;
	double JD0        = 0; // sphere --> JD0 doesn't matter

	double H = 16.5;
	double G = 0.15;
	double emissivity = 0.9;
	double pV = 0.3;
*/

/*
// (22) Kalliope
	fitfiles.push_back(new asteroid::fitFileMJy("R:/Hestroffer, Marchis, HST/22 Kalliope/fitfile_IRAS22_12+25.dat"));
	
	ConvexFile shape ("S:/22 Kalliope/Kalliope.ver.mikko.convex");
	double ecl_lambda =  20.3108817;
	double ecl_beta   =  90 - 112.885576;
	double period_h   = 4.14819973;
	double JD0        = 2436258.63627;

	double H = 6.45; // Horizons
	double G = 0.21; // Horizons
	double emissivity = 0.9;
	double pV = 0.15;
*/
/*
// (21) Lutetia
//	fitfiles.push_back(new asteroid::fitFileMJy("R:/Lutetia/IRAS.fitfile.dat"));
//	fitfiles.push_back(new asteroid::fitFileMJy("R:/Lutetia/IRAS60.fitfile.dat"));
//	fitfiles.push_back(new asteroid::fitFileSI ("R:/Lutetia/MIRSI.fitfile.dat"));
//	fitfiles.push_back(new asteroid::fitFileSI ("R:/Lutetia/MIRSI.fitfile.larger_error.dat"));
	fitfiles.push_back(new asteroid::fitFileSI ("R:/Lutetia/MIRSI.fitfile.June2005.dat"));
	ConvexFile shape ("S:/Lutetia/Lutetia.ver.convex");

	double ecl_lambda = 38.5329802;
	double ecl_beta   = 90-86.9990486;
	double period_h   = 8.16545586;
	double JD0        = 2444822.35116;

	//double H   = 7.284;  // used all the time -- don't recall why
	//double H = 7.297; // AstDys
	double H = 7.35; // IRAS
	double G   = 0.11;
	double emissivity = 0.9;
	double pV  = 0.2;
//	const double gammaDeg = 0;
//	const double gammaDeg = 180;
//	const double gammaDeg = 150;
*/


// Betulia:
//	const char datafilename[] = "R:/betulia/fitfile.all+lebofsky.dat";
//	fitfiles.push_back(new asteroid::fitFileMJy ("R:/betulia/fitfile.paper.dat"));
/*	fitfiles.push_back(
		new asteroid::fitFileMJy("r:/betulia/fitfile.all+lebofsky_less_error_1976+recalibrated+lcaverage.dat"));
	*/
//	const char datafilename[] = "R:/betulia/fitfile.lebofsky.dat";
//	const char datafilename[] = "R:/betulia/fitfile.paper+lebofsky_less_error_1976.dat";
//	const char datafilename[] = "R:/betulia/fitfile.paper+lebofsky_less_error_1976+recalibrated.dat";
//	fitfiles.push_back(new asteroid::fitFileMJy(
//		"R:/betulia/fitfile.all+lebofsky_less_error_1976+recalibrated+lcaverage.dat"));
// Beta, lambda, and period (hours)
//68.077881 133.275794 6.1383641
// epoch of zero time t=0, fi0 then
//2442916.96183, 0.
/*
	double ecl_lambda = 133.275794;
	double ecl_beta = 90-68.077881; // = 21.922119
	//double period_h = 6.1383641;   // Mikko's value
	//double period_h = 6.138356107; // corrected by MM to fit optical data taken by Yan Fernandez on 2 June 2002
	double period_h = 6.1383602379;  // corrected by MM to 'better' fit the same optical data in March 2005 (???)
	double JD0      = 2442916.96183;
	
	double H = 15.1;
	double G = 0.09;
	double emissivity = 0.9;
	double pV = 0.07325; 
	const double gammaDeg = 180;	

	ConvexFile shape ("S:/betulia/betulia.mikko.convex");
*/



//Itokawa:
//	const char datafilename[] = "R:/itokawa/measurements/ito.japan+marco+myself.no5mu.dat";	
//	const char datafilename[] = "R:/itokawa/measurements/ito.japan+marco+myself.dat";	
//	fitfiles.push_back(
//		new asteroid::fitFileSI("R:/itokawa/measurements/ito.japan+marco+myself.larger_error.dat"));	
	fitfiles.push_back(
		new asteroid::fitFileSI("R:/itokawa/measurements/ito.ThM.dat")); // 9 data points by Th. M.
	fitfiles.push_back(
		new asteroid::fitFileSI("R:/itokawa/measurements/ito.sekionly.dat")); // 1 N-band data-point by Sekiguchi, reanalyzed by Th. M.
	fitfiles.push_back(
		new asteroid::fitFileSI("R:/itokawa/measurements/ito.Delbo_Apr8.dat")); // 5 data points by Delbo, reanalyzed by Th. M.
	//fitfiles.push_back(
	//	new asteroid::fitFileSI("R:/itokawa/measurements/jdflux.ito.2004jul10.dat")); // 5 data points by MM, error probably underestimated
	fitfiles.push_back(
		new asteroid::fitFileSI("r:/itokawa/measurements/ito.IRTF.as_in_paper.dat")); // 5 data points by MM, error bars like in M. Mueller et al. (submitted 2004)
	//fitfiles.push_back(
	//	new asteroid::fitFileSI("R:/itokawa/measurements/ito.M_only.dat")); // 1 data point by Ishiguro et al.
//	const char datafilename[] = "R:/itokawa/measurements/ito.japan+marco+myself.noM.larger_error.dat";	
//	const char datafilename[] = "R:/itokawa/measurements/jdflux.ito.2004jul10.dat";
//	const char datafilename[] = "R:/itokawa/measurements/ito.sekionly.dat";
//	const char datafilename[] ="R:/itokawa/measurements/jdflux.ito.2001mar11+14.dat";
//	const char datafilename[] ="R:/itokawa/measurements/jdflux.ito.2001apr8+9.dat";


	//zero time t0, and fi0 then
	//2451933.95456,0.
	//beta, lambda, and period (hrs):
	//185.367481  176.543541  12.1320138
	//double ecl_lambda = 176.543541; 
	//double ecl_beta   = 90 - 185.367481; // Mikko counts beta from the north pole
//	double ecl_lambda = 176.543541+180;
//	double ecl_beta   = -180-(90-185.367481); // change coordinates such that -90<beta<90 (else bodyfix fails)

	double ecl_beta = 90-179.393675; // values for the 'new' model (after Hayabusa meeting)
	double ecl_lambda = 330.957231;

	double period_h   = 12.13237; 
	double H = 19.5; // M. Mueller et al.
	double G = 0.15;
	
	//double H = 19.9; // used by Th. Mueller
	//double G = 0.21; // dto.

	double emissivity = .9;
	//double pV = 0.218416; // radar diameter = 0.358 km
	//double pV = 0.389486;   // best fit pV (iteration 1)
	//double pV = 0.338609;   // best fit PV (iteration 2)
	double pV = 0.34;         // may be iteration 3
	//double pV = 0.1625;     // radar diameter as updated @ meeting: 415m
	//double pV = .183;

	//ConvexFile shape ("S:/25143 Itokawa/mikko/2001_model/itokawa.ver.convex");
	ConvexFile shape ("S:/25143 Itokawa/mikko/itonew.ver.obj.convex"); // use NEW pole position!
	double JD0        = 2451933.95456;
	//ConvexFile shape ("S:/25143 Itokawa/itokawa.radar.big.convex");
	//double JD0          = 2451933.82; // radar! determined 'experimentally' from optical data

    const double phi0=0;


/*
// 2000 PH5
	//ConvexFile shape ("S:/54509 2000 PH5/vmod.obj.convex");
	//ConvexFile shape ("S:/54509 2000 PH5/vmod2.wf.convex"); // new shape, identical spin axis
	// to run on a remote machine:
	//ConvexFile shape ("input/vmod.obj.convex");
	ConvexFile shape ("input/vmod2.wf.convex");
	double ecl_beta = -81;
	double ecl_lambda = 173;
	double phi0 = -38.8101327056 + 90.; // 3rd Euler angle to define spin state at JD0
									    // 90 must be added to conform my definition of Euler angles
	//const double period_h   = 0.202901;
	const double deg_day = 42582.39651;
	const double period_h = 24/(deg_day/360.);
	double JD0 = 2452117.5; // 2001/07/27.0 UT
	JD0 -= 0.5 * 0.00021352 * 1483*1483/deg_day; // to account for YORP between JD0 and JD 2453600

	
	double emissivity = 0.9;

	double H = 22.562; // NEO-Dys 2006/06/15
	double G = 0.15; 

	// Radar volume: 0.001270 km^3 --> D ~ 134.4 m
	//double pV = 0.0923531;
	//double pV = 0.23; // as suggested by my data
	double pV = 0.3; // for TI = 0

    //fitfiles.push_back(
//        new asteroid::fitFileMJy("R:/PH5/IRAC/fitfile_IRAC4.dat"));
   // fitfiles.push_back(
     //   new asteroid::fitFileMJy("R:/PH5/PUI/fitfile_PUIr_IDLphot814.dat"));
	//fitfiles.push_back(
      //  new asteroid::fitFileMJy("R:/PH5/PUI/fitfile_PUIb_IDLphot814.dat"));

	// to run on a remote machine!
	fitfiles.push_back(
        new asteroid::fitFileMJy("input/fitfile_IRAC4.dat"));
    fitfiles.push_back(
        new asteroid::fitFileMJy("input/fitfile_PUIr_IDLphot814.dat"));
	fitfiles.push_back(
        new asteroid::fitFileMJy("input/fitfile_PUIb_IDLphot814.dat"));
*/


//// Eros
//    fitfiles.push_back(
//	   new asteroid::fitFileMJy ("R:/eros/harris_davies.1998/fluxes.N.dat"));
//
//	double ecl_lambda = 17.22; // used astro.pro to convert RA & Dec 11.350deg & 17.216deg into ecliptic coordinates
//	double ecl_beta   = 11.35;
//	double period_h   = 5.2702573080758831260625733693225; // (from Thomas' shape models)
//	double JD0        = 2451545.19887;
//    const double phi0 = 0;
//	ConvexFile shape ("S:/433 Eros/eros001708.tab.obj.convex");
//
//	double H = 10.82;//10.31; // 10.31: literature; 10.82: diameter from volume = 2535km^3, pV=0.29
//	double emissivity = 0.9;
//	double pV = 0.29;
//	double G = 0.181; // leads to Bond albedo = 0.12 (Domingue 2002, Icarus)



	//ConvexFile shape ("S:/spheres/cube6.obj.convex");

// 
// GENERAL PART - good for any old asteroid
//

    SpinState axis (ecl_lambda, ecl_beta, period_h, JD0, phi0);
	TriangulatedConvex aster(shape, axis, H, G, pV, emissivity);

	int nTime = 300;
	int nZ = 25;
	double zMax = 6;
	double accuracyGoal = 0.005;


	// high roughness: rho = 1, f = 1
	const double gammaDeg = 151.75834;
	const double craterDensity = 1;

	//// default roughness (Mueller 99): rho = 0.7, f=0.6
	//const double gammaDeg = 144.59046;  
	//const double craterDensity = 0.6; 

	//// low roughness (Mueller 99): rho = 0.4, f=0.4
	//const double gammaDeg = 117.70116;
	//const double craterDensity = 0.4;

	//// no roughness at all:
	//const double gammaDeg = 0;
	//const double craterDensity = 0;



/*
	// plotting a Chi^2 - shape:
	const double ThermalInertia = 350;
	//const double gammaDeg = 

	list<ThermalInertiaOnlyConvex*> inertias;
	list<InertiaTimesCraterLagerros*> inertia_craters;
	list<double> ThermalParameters;

	for (list<asteroid::fitFileSI*>::const_iterator it = fitfiles.begin(); 
		it != fitfiles.end(); it++)
	{
		ThermalParameters.push_back(
			aster.calculateThermalParameter(ThermalInertia, (**it)[0].Sun2Aster.getLength()));
		inertias.push_back( new ThermalInertiaOnlyConvex(
			emissivity, aster.getBond(), ThermalParameters.back(), nTime, nZ, zMax, accuracyGoal));
		inertia_craters.push_back( new InertiaTimesCraterLagerros(
			emissivity, aster.getBond(), ThermalParameters.back(), nTime, nZ, zMax, accuracyGoal, gammaDeg));
	};
	//ofstream out ("test.dat");
	plot_chi2(aster, inertias, inertia_craters, fitfiles, out, 201, 101, 2);//, 3, 3, 2);
*/

///*
	for (double ThermalInertia = 1200.01; ThermalInertia >= 0; ThermalInertia -= 12.5)
	//for (double ThermalInertia = 0; ThermalInertia < 2501; ThermalInertia += 50)
	//for (double ThermalParameter = 0; ThermalParameter < 35; ThermalParameter += .2)
	{		
		//list<ThermalInertiaOnlyConvex*> inertias;
		list<InertiaTimesCraterLagerros*> inertia_craters;
		list<double> ThermalParameters;

		for (list<asteroid::fitFileSI*>::const_iterator it = fitfiles.begin(); 
			it != fitfiles.end(); it++)
		{
			ThermalParameters.push_back(
				aster.calculateThermalParameter(ThermalInertia, (**it)[0].Sun2Aster.getLength()));
			//inertias.push_back( new ThermalInertiaOnlyConvex(
			//	emissivity, aster.getBond(), ThermalParameters.back(), nTime, nZ, zMax, accuracyGoal));
			inertia_craters.push_back( new InertiaTimesCraterLagerros(
				aster, ThermalParameters.back(), nTime, nZ, zMax, accuracyGoal, gammaDeg
				, craterDensity
				));
		};

		//double ThermalParameter = aster.calculateThermalParameter(ThermalInertia, fit_file[0].Sun2Aster.getLength());
		//ThermalInertiaOnlyConvex inertia          (emissivity, aster.getBond(), ThermalParameter, nTime, nZ, zMax, 
		//									accuracyGoal);
		//InertiaTimesCraterLagerros inertia_crater (emissivity, aster.getBond(), ThermalParameter, nTime, nZ, zMax, 
		//									 accuracyGoal, gammaDeg);

		//out<<"# Using thermal parameter: "<<ThermalParameter<<'\n';
		out<<"# Using a thermal inertia of "
			//<<aster.calculateThermalInertia(ThermalParameter, fit_file[0].Sun2Aster.getLength())
			<<ThermalInertia
			<<"."<<endl;
		//double rho; double sigmaRho=-12;
		double newPV; double sigmaPV=-12; 
		//double correlation = -12;
		//double chi2 = fit (newPV, sigmaPV, rho, sigmaRho, correlation, 
		//	aster, inertias, inertia_craters, fitfiles, out);
		//double chi2 = fit (rho, sigmaRho, aster, inertia, inertia_crater, fit_file, out);
		double chi2 = fit (newPV, sigmaPV, aster, inertia_craters, fitfiles, out);
		
		out<<fixed;
		out<<"# Best fit geometric albedo (used "<<aster.getPV()<<" as start parameter): "<<newPV<<"+-"<<sigmaPV<<'\n';
		/*
		out<<"# Best fit crater density: "<<rho<<"+-"<<sigmaRho<<endl;
		out<<"# corresponding to a thetabar of ";
		double g2rad = vector3::DEGREES*gammaDeg/2;
		if // make 0<rho<1
			(rho < 0) rho = 0;
		else 
			if (rho > 1) rho = 1;
		out<<vector3::RAD*atan( 2*rho/vector3::PI* ( log( (1+sin(g2rad))/cos(g2rad) ) -sin(g2rad)) / (1-cos(g2rad)) );
		out<<" (Lagerros), ";
		out<<vector3::RAD*atan( 2*rho/vector3::PI * (g2rad-sin(g2rad)*cos(g2rad)) / (1-cos(g2rad)*cos(g2rad)) );
		out<<" (MM)\n";
		*/
		out<<"# with a Chi^2 of "<<chi2<<endl;
		//out<<"# Formal correlation between crater density and scale parameter (pV^{-1}): "<<correlation<<endl;
		//out<<aster.calculateThermalInertia(ThermalParameter, fit_file[0].Sun2Aster.getLength())
		out<<ThermalInertia
			<<'\t'<<chi2
			<<'\t'<<newPV<<'\t'<<sigmaPV
			//<<'\t'<<rho<<'\t'<<sigmaRho
			//<<'\t'<<correlation
			<<'\n'
			<<'#'<<endl;

        if (newPV > 2.5) newPV = 2.5; // just in case; to prevent prog from bugging
        aster.setPV(newPV); // new feature, added by MM on 2006/07/21
        // this increases the fit stability if pV varies significantly with inertia
        // pV is updated after each fit, so for a 'new' value of thermal inertia,
        // the albedo is used which best fit the previous thermal inertia and so on.
        // Comparison with its precursor shows this method to bring down Chi^2,
        // also even for wildly 'bad' start values of pV, 
        // best-fit albedos for the 4th or 5th thermal inertia are independent of the start value.


		// free memory:
		//list<ThermalInertiaOnlyConvex*>::iterator   it_inertias = inertias.begin();
		list<InertiaTimesCraterLagerros*>::iterator it_craters  = inertia_craters.begin();
//		for ( ; it_inertias != inertias.end(); it_inertias ++, it_craters++)
		for ( ; it_craters  != inertia_craters.end(); it_craters++)
//		for ( ; it_inertias  != inertias.end(); it_inertias++)
		{
			//delete *it_inertias; *it_inertias=0;
			delete *it_craters;  *it_craters =0;
		};
	};
	out<<endl;
	out<<"# Geometric albedo:       "<<pV<<endl;
	out<<"# Craters' opening angle: "<<gammaDeg<<endl;
	//out<<"# Data file used was: "<<datafilename<<endl;
//*/	
	
/*
	ostream& out = cout;
	out<<"# Comparison of model vs. measured data, \n";
	out<<"# using a thermal parameter of "<<ThermalParameter<<" and an opening angle of "<<gammaDeg<<" degrees.\n";
	out<<"# Format: JD / lambda(micron) / phase angle (deg) / fluxes (W/m^2/mu): measured / no craters / 100% craters"<<endl;
	for (int i=0; i<fit_file.nRecords; i++)
	{
		const asteroid::fluxRecord& record = fit_file[i];
		out<<fixed;
		out<<record.JD<<'\t'<<record.lambdaMu<<'\t';
		out<<eclipticVector::phaseDegrees(record.Earth2Aster, record.Sun2Aster)<<'\t';
		out<<scientific<<record.fluxSI<<'\t';
		out<<aster.ThermalFlux(record, inertia)<<'\t';
		out.flush();
		out<<aster.ThermalFlux(record, inertia_crater)<<endl;
	};
*/
	
	//vector3 dA (0.174, 0.985, 1); // local time = 80deg from noon = 40 minutes before sunset / after sunrise, 
									// evening side: y>0, morning side: y<0
/*	vector3 dA    (1,0,0);
	vector3 A2Sun (1,0,0);
	vector3 A2Earth (A2Sun);

	double lambdaMu = 10;
	double rAU = 1;

	ThermalInertiaOnlyConvex inert (emissivity, bond, ThermalParameter, nTime, nZ, zMax, accuracyGoal);
	double factor = inert.CorrectionFactor (A2Sun, A2Earth, dA, lambdaMu, rAU);
	cout<<"Correction factor: "<<factor<<endl;
*/	
	//ConvexFile file ("E:/C++/NewThermalModel/ClassHierarchyFortran/obj/itokawa.radar.big.convex");
//	ConvexFile file ("Y:/betulia.mikko.convex");

//	const double H=15;
//	const double G=0.15;
//	const double pv=0.2;
//	const double eps=0.9;
//	const eclipticVector Sun2Aster   (0,0,1);
//	const eclipticVector Earth2Aster (0,0,.1);
//	const double JD = 0;
//	const double lambdaMu = 20;
//
//	const TriangulatedConvex    asteroid1 (file, SpinState(0, 90, 24), H, G, pv, eps);
//	const HemisphericNoInertia  craters   (asteroid1, 180);
//
//	cout<<asteroid1.ThermalFlux(Sun2Aster, Earth2Aster, JD, lambdaMu, craters);	
//	cout<<craters.CorrectionFactor(1,1,0,10,1)<<endl;

	cout<<"Done!"<<endl;
	return 0;
}
catch (exception& exc)
{
	cerr<<exc.what()<<endl;
	return -1;
}
catch (...)
{
	cerr<<"Caught some non-exception (?)"<<endl;
	return -2;
};
};
