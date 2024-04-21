#include "TriangulatedConvex.h"
#include "HemisphericNoInertia.h"
#include "ThermalInertiaOnlyConvex.h"
#include "InertiaTimesCraterLagerros.h"
#include<iostream>
#include<fstream>
#include"jansky.h"
using namespace std;


int main()
{
//	const char datafilename[] = "R:/itokawa/measurements/ito.japan+marco+myself.no5mu.dat";	
//	const char datafilename[] = "R:/itokawa/measurements/ito.japan+marco+myself.dat";	
//	const char datafilename[] = "R:/itokawa/measurements/ito.japan+marco+myself.larger_error.dat";	
//	const char datafilename[] = "R:/itokawa/measurements/ito.japan+marco+myself.noM.larger_error.dat";	
//	const char datafilename[] = "R:/itokawa/measurements/jdflux.ito.2004jul10.dat";
//	const char datafilename[] = "R:/itokawa/measurements/ito.sekionly.dat";
	//asteroid::fitFileSI  fit_file(datafilename);
//	const char datafilename[] ="R:/itokawa/measurements/jdflux.ito.2001mar11+14.dat";
//	const char datafilename[] ="R:/itokawa/measurements/jdflux.ito.2001apr8+9.dat";


	//const char datafilename[] = "R:/betulia/fitfile.all+lebofsky.dat";
//	const char datafilename[] = "R:/betulia/fitfile.paper.dat";
//	const char datafilename[] = "R:/betulia/fitfile.lebofsky.dat";
//	const char datafilename[] = "R:/betulia/fitfile.paper+lebofsky_less_error_1976.dat";
//	const char datafilename[] = "R:/betulia/fitfile.paper+lebofsky_less_error_1976+recalibrated.dat";
//	const char datafilename[] = "R:/betulia/fitfile.all+lebofsky_less_error_1976+recalibrated+lcaverage.dat";
//	const char datafilename[] = "e:/results/backup/eros/harris_davies.1998/fluxes.N.dat";
	//asteroid::fitFileMJy fit_file(datafilename);
	//const unsigned int which_line = 13; // geometry of which data line is to be used? Zero-based!
	const unsigned int which_line = 0;
	//const unsigned int which_line = 27;



	//ConvexFile shape ("S:/spheres/cube6.obj.convex");


//Itokawa:
	//zero time t0, and fi0 then
	//2451933.95456,0.
	//beta, lambda, and period (hrs):
	//185.367481  176.543541  12.1320138
	//double ecl_lambda = 176.543541; 
	//double ecl_beta   = 90 - 185.367481; // Mikko counts beta from the north pole
//	double ecl_lambda = 176.543541+180;
//	double ecl_beta   = -180-(90-185.367481); // change coordinates such that -90<beta<90 (else bodyfix fails)
/*	
	//ConvexFile shape ("S:/itokawa/mikko/2001_model/itokawa.ver.convex");
	//ConvexFile shape ("S:/itokawa/mikko/itonew.ver.obj.convex"); // use NEW pole position!
	//double JD0        = 2451933.95456;
	//ConvexFile shape ("S:/itokawa/itokawa.radar.big.convex");
	//double JD0          = 2451933.82; // radar! determined 'experimentally' from optical data

	double ecl_beta = 90-179.393675; // values for the 'new' model (after Hayabusa meeting)
	double ecl_lambda = 330.957231;

	double period_h   = 12.13237; 
	//double H = 19.5;
	//double G = 0.15;
	
	double H = 19.9; // used by Th. Mueller
	double G = 0.21; // dto.

	double emissivity = .9;
	//double pV = 0.218416; // radar diameter = 0.358 km
	//double pV = 0.389486;   // best fit pV (iteration 1)
	//double pV = 0.338609;   // best fit PV (iteration 2)
	//double pV = 0.34;         // may be iteration 3
	//double pV = 0.1625;     // radar diameter as updated @ meeting: 415m
	double pV = .183;
	const double f = 0.6;
	double th_inertia = 1000; // compare with Th. Mueller
*/
// Betulia:
// Beta, lambda, and period (hours)
//68.077881 133.275794 6.1383641
// epoch of zero time t=0, fi0 then
//2442916.96183, 0.
/*
	//ConvexFile shape ("S:/betulia/betulia.mikko.convex");
	double ecl_lambda = 133.275794;
	double ecl_beta = 90-68.077881; // = 21.922119
	//double period_h = 6.1383641;   // Mikko's value
	//double period_h = 6.138356107; // corrected by MM to fit optical data taken by Yan Fernandez on 2 June 2002
	double period_h = 6.1383602379;  // corrected by MM to 'better' fit the same optical data in March 2005 (???)

	double JD0      = 2442916.96183;
	
	double H = 15.1;
	double G = 0.09;
	double emissivity = 0.9;
	//double pV = 0.0735; 
	//double pV = 0.066;
	//double pV = 0.08;
	//double pV = 0.114;
	double pV = 0.077;
*/
/*
// Eros
	ConvexFile shape ("S:/eros/eros001708.tab.obj.convex");
	double ecl_lambda = 17.22; // used astro.pro to convert RA & Dec 11.350deg & 17.216deg into ecliptic coordinates
	double ecl_beta   = 11.35;
	double period_h   = 5.2702573080758831260625733693225; // (from Thomas' shape models)
	double JD0        = 2451545.19887;

	double H = 10.82;//10.31; // 10.31: literature; 10.82: diameter from volume = 2535km^3, pV=0.29
	double emissivity = 0.9;
	double pV = 0.29;
	double G = 0.181; // leads to Bond albedo = 0.12 (Domingue 2002, Icarus)
*/

// Lutetia
	ConvexFile shape ("S:/Lutetia/Lutetia.ver.convex");

	double ecl_lambda = 38.5329802;
	double ecl_beta   = 90-86.9990486;
	double period_h   = 8.16545586;
	double JD0        = 2444822.35116;
	const double phi0 = 0;

	//double H   = 7.284;  // used all the time -- don't recall why
	//double H = 7.297; // AstDys
	double H = 7.35; // IRAS
	double G   = 0.11;
	double emissivity = 0.9;
	double pV = 0.208; // M. Mueller et al. (2006)
	asteroid::fitFileSI fit_file ("R:/Lutetia/MIRSI.fitfile.larger_error.dat");
	ofstream out ("R:/Lutetia/__10mulc_ThM_defrough.dat");
	const double JD = 2453181;



/*
// 2000 PH5
	ConvexFile shape ("S:/54509 2000 PH5/vmod.obj.convex");
	double ecl_beta = -81;
	double ecl_lambda = 173;
	double JD0 = 2452117.5; // 2001/07/27.0 UT
	double phi0 = -38.8101327056; // 3rd Euler angle to define spin state at JD0
	double period_h = 0.202901;	
	double emissivity = 0.9;
	double H = 22.562; // NEO-Dys 2006/06/15
	double G = 0.15; 
	// Radar volume: 0.001270 km^3 --> D ~ 134.4 m
	double pV = 0.0923531;

    //asteroid::fitFileMJy fit_file("R:/PH5/IRAC/fitfile_IRAC4.dat");
    asteroid::fitFileMJy fit_file("R:/PH5/PUI/fitfile_PUIr_IDLphot814.dat");
*/



// general part, good for any old asteroid
	SpinState axis (ecl_lambda, ecl_beta, period_h, JD0,phi0);
	TriangulatedConvex aster(shape, axis, H, G, pV, emissivity);
	//const double gammaDeg = 0;
	//const double gammaDeg = 144.59; // leading to rho = 0.7 (with f=0.6)

//	TriangulatedConvex::LambertianEmitter model(aster);


	const unsigned int nTime = 300;
	const unsigned int nZ    = 25;
	const double       zMax  = 6;
	const double accuracyGoal = 0.001;
	//const double inertia      = 270;
	const double inertia = 50;
	//const double inertia = 2430;


	const double gammaDeg = //150;
		                    144.59046; 
	const double craterDensity = //0.46;
		                    0.6; // (default roughness)
	ThermalInertiaOnlyConvex model1 (aster, aster.calculateThermalParameter(inertia, fit_file[which_line].Sun2Aster.getLength()),
		nTime, nZ, zMax, accuracyGoal);
	InertiaTimesCraterLagerros model2 (aster.getEmissivity(), aster.getBond(),
		aster.calculateThermalParameter(inertia, fit_file[which_line].Sun2Aster.getLength()),
		nTime, nZ, zMax, accuracyGoal, gammaDeg);

	
/*	
	ThermalInertiaOnlyConvex model (aster, 
		aster.calculateThermalParameter(inertia, fit_file[which_line].Sun2Aster.getLength()),
		nTime, nZ, zMax, accuracyGoal);
*/


/*
	//const double gammaDeg = 150;
	const double gammaDeg = 125.63;
	InertiaTimesCraterLagerros model(aster.getEmissivity(), aster.getBond(),
		aster.calculateThermalParameter(inertia, fit_file[which_line].Sun2Aster.getLength()),
		nTime, nZ, zMax, accuracyGoal, gammaDeg);
*/

	//const double gammaDeg = 125.63;
	//HemisphericNoInertia model(aster, gammaDeg);


	const eclipticVector& Sun2Aster   = fit_file[which_line].Sun2Aster;
	const eclipticVector& Earth2Aster = fit_file[which_line].Earth2Aster;
	//const double JD                   = fit_file[which_line].JD;
	//const double lambdaMu             = 7.872;
    //const double lambdaMu             = 22.3272;
	const double lambdaMu             = 10;
	const unsigned int nPointsEach    = 300;

	//ofstream out ("betulia2002_1stnight@11.7mu_TI260_pV.062_gamma15100perc.dat");
	//ofstream out ("betulia1976@20mu.1stnight.dat");
	//ofstream out ("betulia2002@11.7mu.1stnight.fit76+02.dat");
	//ofstream out ("betulia1976_1stnight@10.6mu.2002-parameters_TI260_pV.062_gamma15_100perc.dat");
	//ofstream out ("betulia1976_1stnight@10.6mu.fitboth_TI240_pV.066_gamma150_30perc.dat");
	//ofstream out ("betulia2002_1stnight@7.91mu.fitboth_TI240_pV.66_gamma150_30perc.dat");
	//ofstream out ("betulia2002_1stnight@10.27mu_TI180_pv.077_gamma150_100perc.dat");
	//ofstream out ("betulia1976_1stnight@10.6mu_TI180_pV.077_gamma150_100perc.dat");
	//ofstream out ("lutetia2004.11.6mu_TI50_pV020_gamma150_100perc.dat");
	//ofstream out ("lutetia2004.8.7mu_TI0_pV024_gamma125.63.dat");

    //ofstream out ("PH5.IRAC4.TI0_gamma125.63.dat");
    //ofstream out ("PH5.IRAC4.TI1500_no_roughness.dat");
    //ofstream out ("PH5.PUIr.TI50_no_roughness.dat");

	//aster.ThermalLightCurve(out, Sun2Aster, Earth2Aster, lambdaMu, JD, model, nPointsEach);
	aster.ThermalLightCurve(out, Sun2Aster, Earth2Aster, lambdaMu, JD, model1, model2, 1.-craterDensity, nZ);
	
	cerr<<"Done!"<<endl;
	return 0;
}; // main
