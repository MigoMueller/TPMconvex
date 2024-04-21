#ifndef __ParmfileLC_H_INCLUDED
#define __ParmfileLC_H_INCLUDED

#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<stdexcept>
#include "../../MM_textutils.h"

class ParmfileLC
{
public:
	ParmfileLC (const char filename[]);
	~ParmfileLC(void){};
	
	std::string outName, shapeName;
	double lambdaMu;
	double JD,
		   SunLambda, SunBeta, rAU,
		   ObsLambda, ObsBeta, deltaAU,
		   axisLambda, axisBeta, periodH, JD0, phi0,
		   H, G, pV, emissivity;
	double TI;
	double accuracyGoalCrater, accuracyGoalTI;
	double gammaDeg, craterDensity;
	double zMax;
	unsigned int nTime, nZ;

protected:
    std::ifstream in;
    std::ostringstream log;
private:
	ParmfileLC(void);
	ParmfileLC (const ParmfileLC&);
	ParmfileLC& operator= (const ParmfileLC&);
protected:

    // If in is OK: Adds "Action" + success note to log. Else: adds to log, throws exception with log
    void checkIn(const std::string& Action);

    void checkReading(const std::string& item)
    {
        std::string dummy = "Reading line containing " + item;
        checkIn(dummy);
    }

    void readItem(std::string& stringvar, const std::string& name)
    {
        skipCommentLines();
        getline(in, stringvar);
        checkReading(name);
    }

    template<class numbertype>
    void readItem(numbertype& number, const std::string& name)
    {
        skipCommentLines();
        in>>number;
        checkReading(name);
    }

    void skipCommentLines(const std::string& delims = "# \n")
        {MM_textutils::skipCommentLines(in, delims);};

    void throwexception()
    {
        const std::string dummy = 
            "Fatal error reading in parameter file. See log for details:\n"
             + log.str();
	//        throw std::exception(dummy.c_str());
        throw std::runtime_error(dummy);
    };

};


#endif
