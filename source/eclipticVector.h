// Little utility class -- a solar system vector's ecliptic coordinates
// includes a method to output a vector (class vector3)

// Conventions:
// x=cos(longitude)*cos(latitude)
// y=sin(longitude)*cos(latitude)
// z=sin(latitude)
//
// i.e.: longitude from 0--360 (or whatever)
// latitude from -90 -- +90 deg

// units: degrees and AU
// also:  length as angular velocity (h^{-1})


#if !defined(AFX_ECLIPTICVECTOR_H__C6D1504F_67B8_46AD_8EDA_E554E2CC41AD__INCLUDED_)
#define AFX_ECLIPTICVECTOR_H__C6D1504F_67B8_46AD_8EDA_E554E2CC41AD__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<stdexcept>
#include "vector3.h"

class eclipticVector  
{
protected:
	// units: degrees and AU
	// also:  length as angular velocity (h^{-1})
	double longitude, latitude, length;
public:
	// accepts negative length, too (for spin axis vectors)
	eclipticVector(double longitude, double latitude, double length_AU=1.) throw ()
		: longitude(longitude), latitude(latitude), length(length_AU)
		{};
	vector3 ecl2vector() const throw()
		{return vector3(longitude, latitude)*=length;};
	operator vector3 () const throw()
		{return ecl2vector();};
	virtual ~eclipticVector() throw() {};

	double getLength() const throw()
		{return length;};
	double getLongitude() const throw()
		{return longitude;};
	double getLatitude() const throw()
		{return latitude;};
	// returns the cosine of the phase angle
	// input parameters are either BOTH pointing away from or to the target
	static double cosPhaseAngle (const eclipticVector& Sun, const eclipticVector& Observer)
		{return Sun.ecl2vector().dot(Observer.ecl2vector())/(Sun.length*Observer.length);};
	// returns the phase angle in degrees
	// input parameters are either BOTH pointing away from or to the target
	static double phaseDegrees (const eclipticVector& Sun, const eclipticVector& Observer)
		{return acos(cosPhaseAngle(Sun, Observer))*vector3::RAD;};
private:
	// prohibitive
	eclipticVector();
};

#endif // !defined(AFX_ECLIPTICVECTOR_H__C6D1504F_67B8_46AD_8EDA_E554E2CC41AD__INCLUDED_)
