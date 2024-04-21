#include "ParmfileLC.h"
#include<string>
using namespace std;


ParmfileLC::ParmfileLC (const char filename[])
{
	in.open(filename);
    {
        string dummy = "Opening parameter file "+ (string)filename;
        checkIn(dummy);
    };
	{
		ostringstream dummy;
		dummy<<filename<<".out.dat"<<ends;
		this->outName = dummy.str();
	};


	this->readItem(shapeName, "shape file name");

	this->readItem(H, "absolute optical magnitude H");
	this->readItem(G, "slope parameter G");
	this->readItem(emissivity, "emissivity (try 0.9)");
	this->readItem(pV, "geometric albedo pV");

	this->readItem(axisLambda, "spin axis orientation: ecl. longitude (deg)");
	this->readItem(axisBeta, "spin axis orientation: ecl. latitude (deg)");
	this->readItem(periodH, "rotation period (hours)");
	this->readItem(JD0, "JD0");
	this->readItem(phi0, "rotational phase at JD0 (deg; try zero)");

	this->readItem(nTime, "number of time steps");
	this->readItem(nZ, "number of depth steps");
	this->readItem(zMax, "maximum depth (in skin depths)");
	this->readItem(accuracyGoalTI, "fractional accuracy goal for thermal conduction");
	this->readItem(accuracyGoalCrater, "fractional accuracy goal for crater integral");

	this->readItem(gammaDeg, "crater opening angle (deg, 0--180)");
	this->readItem(craterDensity, "crater density (0--1)");

	this->readItem(lambdaMu, "wavelength (micron)");
	this->readItem(this->TI, "thermal inertia (SI units)");
	this->readItem(JD, "JD of desired lightcurve");

	this->readItem(SunLambda, "ecl. vector Sun->Asteroid: longitude (deg)");
	this->readItem(SunBeta, "ecl. vector Sun->Asteroid: latitude (deg, -90--90)");
	this->readItem(rAU, "heliocentric distance (AU)");
	this->readItem(ObsLambda, "ecl. vector Observer->Asteroid: longitude (deg)");
	this->readItem(ObsBeta, "ecl. vector Observer->Asteroid: latitude (deg, -90--90)");
	this->readItem(deltaAU, "observer-centric distance (AU)");
	

} // ParmfileLC::ctor




// If in is OK: Adds "Action" + success note to log. Else: adds to log, throws exception with log
void ParmfileLC::checkIn(const string& Action)
{
    log<<Action<<": ";
    if (!in)
    {
        log<<"Unsuccessful, exit!\n";
        throwexception();
    }
    else
        log<<"OK.\n";
}
