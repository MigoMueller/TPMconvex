// LagerrosLookedup:
//
// Reproduces correction factors for hemispherical craters as computed by LookupTables
// and interpolates them
//
//////////////////////////////////////////////////////////////////////

// Author: michael.mueller@dlr.de
// Date:   2004 Mar 16

// Caveat:
// Relies on a well formatted 'gamma.ini' at a position hard-wired in LagerrosLookedup.cpp!
// From this, the further file-and-directory-structure is determined.


// Not supposed to be a base-class for anything, non-virtual dtor!!!


#if !defined(AFX_LAGERROSLOOKEDUP_H__21B45DA7_2FFF_4914_ACFC_044AE91C77E0__INCLUDED_)
#define AFX_LAGERROSLOOKEDUP_H__21B45DA7_2FFF_4914_ACFC_044AE91C77E0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
//using namespace std;
#include "arrays.h"
#include "XFile.h"
#include "XFileManager.h"

class vector3; // defined in vector3.h

class gammaini;		// defined below LagerrosLookedup
class bluenessini;	// dto.


class LagerrosLookedup  
{
private:
	class values; // defined further below
	
protected:
	static const gammaini gammaList;
	static const bluenessini bluenessList;
	static arrays::pointerArray<XFileManager*> managers; // as many as there are opening angles
	static const int xxxxxxx; // dummy 
	static const int initializeManagers(); // purpose: initialize the elements of 'managers'
	
	const int gamma;
	const double blue;
	const XFileManager::twoFilesPlusRemainder bluenessInfo;
	const XFile &X0, &X1;
	const double blueRemainder;// element [0,1]; indicates 'how far we are from X0 towards X1'
	const double sunStep, earthStep, aziStep;
	// returns the XFileManager corresponding to gamma
	inline XFileManager& manager(int gamma) const;
public:
	LagerrosLookedup (int gamma, double blueness);
	~LagerrosLookedup(){};
protected:
	double interpolate (double cosSun, double cosEarth, double cosAzi) const;
	double operator () (double cosSun, double cosEarth, double cosAzi) const
		{return interpolate(cosSun, cosEarth, cosAzi);};
public:
	double CorrectionFactor (const vector3& dA, const vector3& A2Earth, const vector3& A2Sun);
	static const gammaini& getGamma() 
		{return gammaList;};
	static const bluenessini& getBlueness() 
		{return bluenessList;};
private: // prohibitive
	LagerrosLookedup();
	LagerrosLookedup(const LagerrosLookedup&);
	LagerrosLookedup& operator= (const LagerrosLookedup&);
		
	// li'l dummy class to simplify life in LagerrosLookedup::interpolate
	class values
	{
	public:
		int iLow, iHigh;
		double remainder;
		values(double value, double first, double last, double step);
		virtual ~values(){};
	}; // class LagerrosLookedup::values

	// 1D-interpolation
	static inline double interpol(double low, double high, double remainder);
	// 2D-interpolation
	static inline double LagerrosLookedup::interpol (double lowlow, double lowhigh,
						double highlow, double highhigh,
						double remainder1, double remainder2);

}; // class LagerrosLookedup

class gammaini
{
protected:
	static const std::string directoryRoot; // must terminate with a '/'
	int n;
	arrays::array<int>    gammas;
	arrays::array<std::string> names;       // must NOT terminate with a '/'
public:
	// returns the name of the directory corresponding to gamma
	const std::string dir(int gamma) const;
	gammaini (const std::string& gammaName);
	const int getN() const
	{return n;};
private: // prohibitive
	gammaini();
	gammaini(const gammaini&);
	gammaini& operator= (const gammaini&);
	friend LagerrosLookedup;
}; // class gammaini



class bluenessini
{
protected:
	int n;
	int first, last, step;
	arrays::array<int>    blue;
	arrays::array<std::string> names;
public:
	bluenessini (const std::string& fileName);
	const std::string& filename (int blue) const;
	const int getFirst() const
	{return first;};
	const int getLast() const
	{return last;};
	const int getStep() const
	{return step;};
private:
	bluenessini();
	bluenessini(const bluenessini&);
	bluenessini& operator= (const bluenessini&);
	friend class LagerrosLookedup;
}; //class bluenessini






#endif // !defined(AFX_LAGERROSLOOKEDUP_H__21B45DA7_2FFF_4914_ACFC_044AE91C77E0__INCLUDED_)
