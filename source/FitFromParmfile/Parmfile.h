#ifndef __Parmfile_h_included
#define __Parmfile_h_included

#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<stdexcept>
#include "../../MM_textutils.h"


class Parmfile
{
public:
    virtual ~Parmfile(void){};
    Parmfile (const char inName[]);

    std::string shapeName, outName;
    std::vector<std::string> dataName;
    std::vector<bool>        data_in_mJy;
    double H, G, emissivity, pV;
    double ecl_beta, ecl_lambda, period_h, JD0, phi0;
    int nTime, nZ;
    double zMax, accuracyGoalTI, accuracyGoalCrater;
    double gammaDeg, craterDensity;
    double TImax, TImin, TIstep;
protected:
    std::ifstream in;
    std::ostringstream log;
private:
    Parmfile(void);
    Parmfile(const Parmfile&);
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
        throw std::runtime_error(dummy.c_str());
    };

};

#endif
