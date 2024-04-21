#include<iostream>
#include<fstream>
#include<list>
#include "Parmfile.h"
#include "../asteroid.h"
#include "../TriangulatedConvex.h"
#include "../InertiaTimesCraterLagerros.h"
#include<sstream>
#include "../jansky.h"

using namespace std;


// returns the best-fit chi^2
// stores best-fit geometric albedo in newPV
// keeps cratering constant (only one model as input)
// Method: linear regression
// minimizes Chi^2=\sum_k [ (modelflux_k-measuredflux_k)/sigma_k ]^2,
// assuming flux ~ 1/pV (pretty good approximation, handle with care, nevertheless!)
// (see also discussion in dissertation)
template<class model_type>
double fit(double& newPV, double& sigmaPV, const TriangulatedConvex& aster,
	   const list<model_type*>& models,
	   const list<asteroid::fitFileSI*>& data,
	   ostream& important_output,
	   ostream& diagnostic_output = cerr) 
{
  ostream& out = important_output;
  ostream& err = diagnostic_output;
  //out<<"# Fitting albedo with fixed cratering\n";
  out<<"#\tJD\tlambda\tModel flux\tmeasured flux\tsigma"<<endl;
  
  double sm=0, ss=0, mm=0; //partial sums 
  // ss: synthetic*measured, ss: synthetic^2, mm: measured^2 
  // all sums /= sigma
  double s,m; // dummy variables
  
  list<asteroid::fitFileSI*>::const_iterator it_fitfiles = data.begin();
  typename list<model_type*>::const_iterator it_models = models.begin();
  
  for ( ; it_fitfiles != data.end(); it_fitfiles++, it_models++)
    {
      const asteroid::fitFileSI&                    data    = **it_fitfiles;
      const TriangulatedConvex::ThermalModelConvex& model   = **it_models;
      
      for (unsigned int i=0; i<data.nRecords; i++)
	{
	  const asteroid::fluxRecord& record  = data[i];
	  out<<"#\t"<<fixed<<record.JD<<'\t'<<record.lambdaMu<<'\t'; out.flush();
	  
	  m = record.fluxSI;
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
}; // fit function




int main(int nargs, char* argsv[])
{
  if (nargs != 2)
    {
      cerr<<"Specify parameter file!"<<endl;
      return -1;
    }


  try{
    cerr<<"Reading in parameter file "<<argsv[1]<<":\n";
    const Parmfile parms (argsv[1]);
    cerr<<"Done."<<endl;
    
    cerr<<"Opening output file "<<parms.outName<<":\n";
    ofstream out(parms.outName.c_str());
    if (!out)
      {
	string dummy = "Fatal error: Couldn't open "+parms.outName+" for writing!\n";
	throw runtime_error (dummy);
      };
    cerr<<"Done.\n";

    if (parms.data_in_mJy.size() > 1)
      cerr<<"Opening "<<parms.dataName.size()<<" data files:\n";
    else
      cerr<<"Opening 1 data file:\n";
    list<asteroid::fitFileSI*> fitfiles;
    for (vector<double>::size_type i=0; i<parms.data_in_mJy.size(); i++)
      {
	cerr<<"\t"<<parms.dataName[i]<<endl;
	if (parms.data_in_mJy[i])
	  fitfiles.push_back(new asteroid::fitFileMJy(parms.dataName[i].c_str()));
	else
	  fitfiles.push_back(new asteroid::fitFileSI(parms.dataName[i].c_str()));
      };
    cerr<<"Done.\n";
    
    cerr<<"Opening shape file "<<parms.shapeName<<":\n";
    ConvexFile shape (parms.shapeName.c_str());
    cerr<<"Done.\n";

    const double ecl_lambda = parms.ecl_lambda,
      ecl_beta   = parms.ecl_beta,
      period_h   = parms.period_h,
      JD0        = parms.JD0,
      phi0       = parms.phi0, 
      
      H          = parms.H,
      G          = parms.G,
      pV         = parms.pV,
      emissivity = parms.emissivity,
      
      TImax = parms.TImax, 
      TImin = parms.TImin, 
      TIstep = parms.TIstep, 
      
      accuracyGoalCrater = parms.accuracyGoalCrater,
      accuracyGoalTI = parms.accuracyGoalTI,
      
      gammaDeg = parms.gammaDeg,
      craterDensity = parms.craterDensity,
      
      zMax = parms.zMax;
    const unsigned int
      nTime = parms.nTime, 
      nZ = parms.nZ;
    
    cerr<<"Generating asteroid model:\n";
    SpinState axis (ecl_lambda, ecl_beta, period_h, JD0, phi0);
    TriangulatedConvex aster(shape, axis, H, G, pV, emissivity);
    cerr<<"Done.\n";
    
    out<<"# Generated using FitFromParmfile\n"
       <<"# Contains TPM fits to thermal-IR asteroid data\n"
       <<"# Parameter file used:  "<<argsv[1]<<"\n#\n#\n#\n";
    
    for (double ThermalInertia = TImax; ThermalInertia >= TImin; ThermalInertia -= TIstep)
      {
	cerr<<"Start fitting to TI = "<<ThermalInertia<<":\n";
	
	list<InertiaTimesCraterLagerros*> inertia_craters;
	list<double> ThermalParameters;
	for (list<asteroid::fitFileSI*>::const_iterator it = fitfiles.begin(); 
	     it != fitfiles.end(); it++)
	  {
	    ThermalParameters.push_back(
					aster.calculateThermalParameter(ThermalInertia, (**it)[0].Sun2Aster.getLength()));
	    inertia_craters.push_back( new InertiaTimesCraterLagerros(
			  aster, ThermalParameters.back(), nTime, nZ, zMax, accuracyGoalTI, gammaDeg
			  , craterDensity
		     ));
	    inertia_craters.back()->setAccuracies(accuracyGoalTI, accuracyGoalCrater);
	  }
	out<<"# Using a thermal inertia of "
	   <<ThermalInertia
	   <<":"<<endl;
	double newPV; double sigmaPV=-12; 
	double chi2 = fit (newPV, sigmaPV, aster, inertia_craters, fitfiles, out);
	out<<fixed;
	out<<"# Best fit geometric albedo (used "<<aster.getPV()<<" as start parameter): "<<newPV<<"+-"<<sigmaPV<<'\n';
	out<<"# with a Chi^2 of "<<chi2<<endl;
	double diameter = 1329 * pow(10., -parms.H/5.)/sqrt(newPV),
	  diameterUnc = diameter * 0.5 * sigmaPV/newPV;
	out<<ThermalInertia
	   <<'\t'<<chi2
	   <<'\t'<<newPV<<'\t'<<sigmaPV
	   <<'\t'<<diameter<<'\t'<<diameterUnc
	   <<'\n'
	   <<'#'<<endl;
	if (newPV > 2.5) newPV = 2.5; // just in case; to prevent prog from bugging
	aster.setPV(newPV); 
	// release memory:
	list<InertiaTimesCraterLagerros*>::iterator it_craters  = inertia_craters.begin();
	for ( ; it_craters  != inertia_craters.end(); it_craters++)
	  {
	    delete *it_craters;  
	    *it_craters =0;
	  }
	cerr<<"Done.\n";
      } // for-loop (ThermalInertia)
    out<<endl;
    out<<"# Geometric albedo:       "<<pV<<endl;
    out<<"# Craters' opening angle: "<<gammaDeg<<endl;

    cout<<"Program finishes regularly. Bye bye!"<<endl;
    return 0;
  } // end of try block
  catch(exception& exc)
    {
      cerr<<exc.what();
      return -2;
    }
  catch(...)
    {
      cerr<<"Unknown fatal error occurred, exit."<<endl;
      return -3;
    };
}; // main
