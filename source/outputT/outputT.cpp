// outputT.cpp
// M.Mueller@sron.nl, 2014/05/29
// Read in shape file, output temperatures per facet for each point in rotational lightcurve.
// Goal: combine T maps, shape model to generate videos using POVray (using script)

#include<iostream>
#include<fstream>
#include<list>
#include<cstdlib>
#include "../TriangulatedConvex.h"
#include "../ThermalInertiaOnlyConvex.h"
#include<sstream>
#include<iomanip>  // for setprecision

using namespace std;





int main(int nargs, char* argsv[])
{
  if (nargs != 19)
    {
      cerr<<argsv[0]<<endl;
      cerr<<"Please provide 11 command-line parameters:"<<endl;
      cerr<<"1  output file name"<<endl;
      cerr<<"2  shape model"<<endl;
      cerr<<"3  spinEclLon"<<endl;
      cerr<<"4  spinEclLat"<<endl;
      cerr<<"5  periodH"<<endl;
      cerr<<"6  JD0"<<endl;
      cerr<<"7  JD"<<endl;
      cerr<<"8  eclLonSun"<<endl;
      cerr<<"9  eclLatSun"<<endl;
      cerr<<"10 eclLonObs"<<endl;
      cerr<<"11 eclLatObs"<<endl;
      cerr<<"12 rAU"<<endl;
      cerr<<"13 TI"<<endl;
      cerr<<"14 pV"<<endl;
      cerr<<"15 nTimes (try 300)"<<endl;
      cerr<<"16 zMax (try 6)"<<endl;
      cerr<<"17 nDepthSteps (try 25)"<<endl;
      cerr<<"18 accuracyGoal (try 5e-3)"<<endl;
      return -1;
    }

  try{
    cerr<<"Opening shape file "<<argsv[2]<<endl;
    ConvexFile shape(argsv[2]);
    cerr<<"Opening output file "<<argsv[1]<<endl;
    ofstream out(argsv[1]);
    if (!out)
      {
	string dummy = "Fatal error: Couldn't open output file for writing!\n";
	throw runtime_error (dummy);
      };
    cerr<<"Parsing spin-axis parameters"<<endl;
    const double spinEclLon = atof(argsv[3]), 
      spinEclLat = atof(argsv[4]),
      periodH = atof(argsv[5]),
      JD0 = atof(argsv[6]);
    cerr<<"Instantiating spin axis"<<endl;
    SpinState axis (spinEclLon, spinEclLat, periodH, JD0, 0); // hard-wire phi0 to be 0
    cerr<<"Parsing obs parameters"<<endl;
    const double JD = atof(argsv[7]),
      eclLonSun=atof(argsv[8]),
      eclLatSun=atof(argsv[9]),
      eclLonObs=atof(argsv[10]),
      eclLatObs=atof(argsv[11]),
      rAU=atof(argsv[12]),
      TI=atof(argsv[13]),
      pV=atof(argsv[14]);
    cerr<<"Instantiating asteroid"<<endl;
    TriangulatedConvex aster(shape, axis, 15, 0.15, pV, 0.9);
    cerr<<"Instantiating ecliptic vectors"<<endl;
    const eclipticVector Sun2Aster(eclLonSun, eclLatSun, rAU),
      Obs2Aster(eclLonObs, eclLatObs, 1); // hard-wire deltaAU to be 1
    //    const unsigned int nTimes = atoi(argsv[3]);
    cerr<<"Parsing thermal-model parameters"<<endl;
    const unsigned int nTimes=atoi(argsv[15]);
    const double zMax = atof(argsv[16]);
    const unsigned int nDepthSteps=atoi(argsv[17]);
    const double accuracyGoal = atof(argsv[18]);
    cerr<<"Instantiating thermal model"<<endl;
    ThermalInertiaOnlyConvex model(aster, aster.calculateThermalParameter(TI, rAU), 
				   nTimes, nDepthSteps, zMax, accuracyGoal);
    cerr<<"Starting calculations..."<<endl;
    vector<vector3> A2Sun, A2Obs; // will be set by getULightCurves
    vector<double> JDs;
    vector<vectorN> temperatures = aster.getULightCurves(Sun2Aster, Obs2Aster, JD, model, A2Sun, A2Obs, JDs);
    cerr<<"Starting output"<<endl;
    out<<"# Output of "<<argsv[0]<<'\n';
    out<<"# Format: an explanatory line starting in '#', followed by actual value(s)\n";
    out<<"# nTimes\n"<<nTimes<<endl;
    out<<"# shape file name\n"<<argsv[2]<<endl;
    out<<"# number of facets\n"<<shape.dA.size()<<endl;
    out<<"# Cartesian vectors (in asteroid-centric system) of direction to Sun: x,y,z; x,y,z; x,y,z; ...\n";
    out<<fixed<<setprecision(4)<<A2Sun[0].x<<", "<<A2Sun[0].y<<", "<<A2Sun[0].z;
    for (int i=1; i<nTimes; i++)
      {
	out<<"; "<<A2Sun[i].x<<", "<<A2Sun[i].y<<", "<<A2Sun[i].z;
      };
    out<<endl;
    out<<"# Cartesian vectors (in asteroid-centric system) of direction to Observer (as above)\n";
    out<<fixed<<setprecision(4)<<A2Obs[0].x<<", "<<A2Obs[0].y<<", "<<A2Obs[0].z;
    for (int i=1; i<nTimes; i++)
      {
	out<<"; "<<A2Obs[i].x<<", "<<A2Obs[i].y<<", "<<A2Obs[i].z;
      };
    out<<endl;
    out<<"# JD of lightcurve points\n";
    out<<fixed<<setprecision(6)<<JDs[0];
    for (int i=1; i<nTimes; i++)
      {
	out<<"; "<<JDs[i];
      };
    out<<endl;
    out<<"# sub-solar temperature TSS (K)\n"<<model.getTSS(rAU)<<endl;
    out<<"# dimensionless temperatures (T/TSS) over lightcurve: one line per facet: facet # (zero-based); T0, T1, ...\n";
    for (int i=0; i<temperatures.size(); i++)
      {
	out<<i<<"; "<<fixed<<setprecision(5)<<temperatures[i][0];
	for (int j=1; j<temperatures[i].size(); j++)
	  {
	    out<<", "<<temperatures[i][j];
	  };
	out<<endl;
      };

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



