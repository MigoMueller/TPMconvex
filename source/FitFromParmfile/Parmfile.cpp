#include "Parmfile.h"
//#include<sstream>
using namespace std;

Parmfile::Parmfile (const char inName[]) 
{
    //// Lines starting with one of the commentDelims will be ignored
    //static const string commentDelims = "#";

    in.open(inName);
    {
        string dummy = "Opening parameter file "+ (string)inName;
        checkIn(dummy);
    };

	{
		ostringstream dummy;
		dummy<<inName<<".out.dat"<<ends;
		this->outName = dummy.str();
	};

    readItem(shapeName, "shape file name");
    readItem(H, "H value");
    readItem(G, "G value");
    readItem(emissivity, "emissivity");
    readItem(pV, "initial pV");
    readItem(ecl_lambda, "spin axis ecliptic lambda (deg)");
    readItem(this->ecl_beta, "spin axis ecliptic beta (deg)");
    readItem(this->period_h, "spin period (hours)");
    this->readItem(this->JD0, "JD_0");
    this->readItem(this->phi0, "phi_0");

    this->readItem(this->nTime, "# time steps");
    this->readItem(this->nZ,    "# depth steps");
    this->readItem(this->zMax,  "max depth (in skin lengths)");
    this->readItem(this->accuracyGoalTI, "fract. accuracy goal (conduction)");
    this->readItem(this->accuracyGoalCrater, "fract. accuracy goal (crater)");
    this->readItem(this->gammaDeg, "crater opening angle (deg; between 0 and 180)");
    this->readItem(this->craterDensity, "crater density (between 0 and 1)");
    this->readItem(this->TImin, "minimum thermal inertia (SI units)");
    this->readItem(this->TImax, "maximum thermal inertia (SI units)");
    this->readItem(this->TIstep, "step width in thermal inertia (SI units)");

    unsigned int nData;
    this->readItem(nData, "# data files");
    this->dataName.resize(nData);
    this->data_in_mJy.resize(nData);
    for (unsigned int i=0; i<nData; i++)
    {
        ostringstream dummy;
        dummy<<"data file #"<<(i+1)<<" of "<<nData;
        //string dummy = "data file #" + (string)(i+1) + " of " 
        this->readItem(dataName[i], dummy.str());
        string SImJy;
        this->readItem(SImJy, "flag for SI or mJy");
        if (SImJy == "SI")
            this->data_in_mJy[i] = false;
        else
        {
            if (SImJy == "mJy")
                this->data_in_mJy[i] = true;
            else
            {
                log<<"Exit: Flag is '"<<SImJy<<"'; must be either 'SI' of 'mJy' (case sensitive)\n";
                throwexception();
            };
        }
    } // (for (i=0..nData)


    in.close();
    this->checkIn("Closing parameter file");
    //ofstream outLog("ParmfileLog.dat");
    //outLog<<log.str();
};


// If in is OK: Adds "Action" + success note to log. Else: adds to log, throws exception with log
void Parmfile::checkIn(const string& Action)
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
