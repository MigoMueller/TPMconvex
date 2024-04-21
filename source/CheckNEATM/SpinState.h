// Spin State contains:
// eclipticVector (inherited): spin axis in ecliptic coordinates 
//							   (length = angular velocity / h^-1)
// double JD0 = JD for zero time
// double phi0= rotational phase for zero time (deg)

// This definition pretty much follows that of Mikko Kaasalainen!

// Autor: michael.mueller@dlr.de
// Date:  2004 Feb 4


#if !defined(AFX_SPINSTATE_H__677C8B19_4502_4F07_A401_C01B1BB4F5B1__INCLUDED_)
#define AFX_SPINSTATE_H__677C8B19_4502_4F07_A401_C01B1BB4F5B1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "eclipticVector.h"

class SpinState : public eclipticVector  
{
public:
	SpinState(double longitude, double latitude, double periodH, double JD0=0, double phi0=0)
		: eclipticVector(longitude, latitude, vector3::PI2/periodH),
		  JD0(JD0), phi0(phi0), periodH(periodH)
		  {};

	SpinState (const eclipticVector& ecl, double JD0=0, double phi0=0)
		: eclipticVector(ecl), JD0(JD0), phi0(phi0), periodH(vector3::PI2/ecl.getLength())
		{};

	double getJD0() const
		{return JD0;};
	double getTh() const
		{return periodH;};
	double getOmegaH() const
		{return length;};
	double getTDays() const
		{return periodH/24.;};
	double getOmegaDays() const
		{return 24*length;};

	virtual ~SpinState(){};
protected:
	const double JD0;
	const double phi0;
	const double periodH;
private:
	SpinState();
};

#endif // !defined(AFX_SPINSTATE_H__677C8B19_4502_4F07_A401_C01B1BB4F5B1__INCLUDED_)
