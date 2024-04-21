// checkNEATM.cpp: implementation of the checkNEATM class.
//
//////////////////////////////////////////////////////////////////////

#include<fstream>
#include<direct.h>
#include "checkNEATM.h"
#include "NEATMFile.h"
using namespace std;

checkNEATM::checkNEATM(const string &directoryRoot)
	: filenames(directoryRoot)
{
	const string outName = nameserver::verifyDir(directoryRoot)+"NEATMfit.dat";
	ofstream out(outName.c_str());
	//ostream out=cout;
	if (!out)
		throw logic_error ("Couldn't open " + outName + " for writing!");

	out<<"# NEATM-fits \n";
	out<<"# Format: \n";
	out<<"# 1st line: # craters' opening angle (deg), crater density (percent) \n"
	   <<"# next lines: phase angle (degrees), best-fit eta, best-fit pV \n"
	   <<"# and so forth until no more craters are available :) \n";
	out<<"#\n";

	for (multimap<int,int>::const_iterator it=filenames.GammaDensity.begin(); it!=filenames.GammaDensity.end(); it++)
	{
		const int openingAngle   = it->first;
		const int densityPercent = it->second;
		out<<"# "<<openingAngle<<'\t'<<densityPercent<<'\n';
		for (list<int>::const_iterator it=filenames.alphaList.begin(); it!= filenames.alphaList.end(); it++)
		{
			const int alphaDeg = *it;
			NEATMFile NEATM ( filenames(alphaDeg, openingAngle, densityPercent).c_str());
			out<<alphaDeg<<'\t';
			NEATMFile::etaPV fit = NEATM.fit(.02);
			out<<fit.eta<<'\t'<<fit.pV<<'\n';
		};
		out<<endl;
	}; // done with GammaDensity
}; // checkNEATM::ctor



checkNEATM::~checkNEATM()
{CleanUp();}; // checkNEATM::dtor


inline void checkNEATM::CleanUp()
{
//	delete NEATM;
//	NEATM=0;
}; //checkNEATM::CleanUp



checkNEATM::nameserver::nameserver(const string& directoryRoot)
	: directoryRoot(verifyDir(directoryRoot)),
	  alpha(this->directoryRoot+"alpha.ini"), 
	  craters(this->directoryRoot+"craters.ini"),
	  alphaList (alpha.getAlphas()),
	  GammaDensity (craters.getGammaDensity())
{}; // checkNEATM::nameserver::ctor



const string checkNEATM::nameserver::operator () 
	(int phaseDeg, int openingAngleDeg, int craterDensityPercent) const
{
	return directoryRoot + // has '/' at the end
		   verifyDir(alpha.dirName(phaseDeg)) +  // now has '/' at the end, too
		   craters.fileName(openingAngleDeg, craterDensityPercent);
}; // checkNEATM::nameserver::operator ()



inline const string checkNEATM::nameserver::verifyDir(const string &rhs) 
{
	if (rhs.empty())
		return rhs;
	// there IS a last letter since rhs is not empty --> foolproof
	if (*(rhs.end()-1) == '/')
		return rhs;
	else
		return rhs+'/';
}; // checkNEATM::nameserver::verifyDir



const string& checkNEATM::nameserver::cratersIni::fileName(int openingAngleDeg, int craterDensity) const
{
	const_iteratorO itO = data.find(openingAngleDeg);
	if (itO == data.end())
	{
		ostringstream dummy;
		dummy<<"There is no opening angle of "<<openingAngleDeg<<" degrees in our list!"<<ends;
		throw invalid_argument(dummy.str());
	}; // opening angle successfully looked up
	
	const_iteratorI itI = itO->second.find(craterDensity);
	if (itI == itO->second.end())
	{
		ostringstream dummy;
		dummy<<"There is no crater density of "<<craterDensity<<"% in our list!"<<ends;
		throw invalid_argument(dummy.str());
	};
	return itI->second;
}; // checkNEATM::cratersIni::fileName



const string& checkNEATM::nameserver::alphaIni::dirName(int alpha) const
{
   const const_iterator it=data.find(alpha);
   if (it==data.end()) // alpha not found
   {
	   ostringstream dummy;
	   dummy<<"Couldn't find "<<alpha<<" in list of phase angles!"<<ends;
	   throw invalid_argument(dummy.str());
   };
   return it->second;
}; // checkNEATM::nameserver::alphaIni::dirName



checkNEATM::nameserver::alphaIni::alphaIni(const string& fileName)
{
	ifstream in(fileName.c_str());
	if (!in)
	{
		stringstream dummy;
		dummy<<"Couldn't open "<<fileName<<" for reading!"<<ends;
		throw invalid_argument(dummy.str().c_str());
	};
	in>>H>>G>>pV>>rAU>>deltaAU;
	if (!in)
	{
		stringstream dummy;
		dummy<<"Error in header of "<<fileName<<".\nMust contain: H, G, pV, r(AU), delta(AU)!"<<ends;
		throw invalid_argument(dummy.str().c_str());
	};
	int alpha;
	string dirName;
	while (!in.eof())
	{
		in>>alpha>>dirName;
		if (!in)
		{ // something went wrong or eof
			if (!in.eof())
			{
				stringstream dummy;
				dummy<<"Error in "<<fileName<<": malformatted data line "<<data.size()+1<<"!"<<ends;
				throw invalid_argument(dummy.str().c_str());
			};
		}
		else
		{
			resultType dummy = data.insert(dataType(alpha, dirName));
			if (!dummy.second) // insert failed
			{
				stringstream dummy;
				dummy<<"Error in "<<fileName<<", data line "<<data.size()+1<<": Each alpha must occur only once!"<<ends;
				throw invalid_argument(dummy.str().c_str());
			};
		};
	}; // while (!in.eof())
	if (data.empty())
	{
		stringstream dummy;
		dummy<<fileName<<" doesn't contain data lines!"<<ends;
		throw invalid_argument(dummy.str().c_str());
	};
}; // checkNEATM::alphaIni::ctor


 
checkNEATM::nameserver::cratersIni::cratersIni(const string& fileName)
{
	ifstream in(fileName.c_str());
	if (!in)
	{
		stringstream dummy;
		dummy<<"Couldn't open "<<fileName<<" for reading!"<<ends;
		throw invalid_argument (dummy.str().c_str());
	};
	int angle, density;
	string name;
	while (!in.eof())
	{
		in>>angle>>density>>name;
		if (in)
		{
			resultTypeO dummyO = data.insert(dataTypeO (angle, innerMap()));
			resultTypeI dummyI = dummyO.first->second.insert (dataTypeI (density, name));
			if (!dummyI.second)
			{
				ostringstream dummy;
				dummy<<"Inserting opening angle of "<<angle<<" and crater density of "<<density<<" failed!"<<ends;
				throw invalid_argument(dummy.str().c_str());
			};
			// everything worked fine
		}
		else if (!in.eof())
		{
			stringstream dummy;
			dummy<<fileName<<" is malformatted, bail out!"<<ends;
			throw invalid_argument (dummy.str().c_str());
		};
	}; // while (!in.eof())
}; // checkNEATM::cratersIni::ctor



std::list<int> checkNEATM::nameserver::alphaIni::getAlphas() const
{
	std::list<int> dummy; 
	for (const_iterator it=data.begin(); it!=data.end(); it++)
		dummy.push_back(it->first);
	return dummy;
};



multimap<int,int> checkNEATM::nameserver::cratersIni::getGammaDensity() const
{
	multimap<int,int> dummy;
	for (const_iteratorO it=data.begin(); it!=data.end(); it++)
	{
		const int gamma = it->first;
		const innerMap& densities=it->second;
		for (const_iteratorI it=densities.begin(); it!=densities.end(); it++)
			dummy.insert(pair<int,int>(gamma, it->first));
	};
	return dummy;
};
