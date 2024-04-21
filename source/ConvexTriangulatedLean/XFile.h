// XFile:
//
// Contents of one file produced by "LookupTables":
// Correction factors for a hemispherical crater for different aspects
//
//////////////////////////////////////////////////////////////////////

// Author: michael.mueller@dlr.de
// Date:   2004 Mar 16

// Not supposed to be a base-class for anything, non-virtual dtor!!!

#if !defined(AFX_XFILE_H__1C847101_2DBF_4FC9_AE17_257CE0C07F48__INCLUDED_)
#define AFX_XFILE_H__1C847101_2DBF_4FC9_AE17_257CE0C07F48__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "arrays.h"


class XFile  
{
protected:
	int gamma, blueness, nSun, nEarth, nAzi;
	double sunStep, earthStep, aziStep;
	arrays::array3D<double> field;
public:
	XFile(const char filename[]);

	const double& operator() (int iSun, int iEarth, int iAzi) const
		{return field[iSun][iEarth][iAzi];};
	double getSunStep() const 
		{return sunStep;};
	double getEarthStep() const 
		{return earthStep;};
	double getAziStep() const 
		{return aziStep;};
	// craters' opening angle in degrees (180 deg = hemispheres)
	int getGamma() const
		{return gamma;};
	// blueness = TSS*lambda (K*micron) -- Planck curve peaks at about 3000 (Wien)
	int getBlueness () const
		{return blueness;};
protected:
	double& entry(int iSun, int iEarth, int iAzi)
		{return field[iSun][iEarth][iAzi];};
private:
	XFile();
	XFile(const XFile&);
	XFile& operator= (const XFile&);
public:
	~XFile(){};
};

#endif // !defined(AFX_XFILE_H__1C847101_2DBF_4FC9_AE17_257CE0C07F48__INCLUDED_)
