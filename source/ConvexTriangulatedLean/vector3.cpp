// Implementation of vector3.h
// (mainly: class vector3 + corresponding functions)


#include<cmath>
#include "vector3.h"
using namespace std;

const double vector3::PI=4*atan(1.);
const double vector3::PI2=2*vector3::PI;
const double vector3::DEGREES=vector3::PI/180.;
const double vector3::RAD=180./vector3::PI;


// range of latitude: -90 (South pole) -> 0 (equator) -> 90 (North)
// longitude=0--> y=0 (x>0)
vector3::vector3(double longitude, double latitude)
throw()
{	
	longitude*=DEGREES; latitude *= DEGREES; 
	// deg -> rad
	double cosine=cos(latitude);
	x=cos(longitude)*cosine;
	y=sin(longitude)*cosine;
	z=sin(latitude);
}; /* vector3::vector3(double longitude, double latitude) */


// unusual ctor, yields cross-product of lhs and rhs (in this order)
vector3::vector3 (const vector3& lhs, const vector3& rhs) 
throw()
{
	x=lhs.y*rhs.z - lhs.z*rhs.y;
	y=lhs.z*rhs.x - lhs.x*rhs.z;
	z=lhs.x*rhs.y - lhs.y*rhs.x;
};
		



void vector3::rxcw(double rad)
{
	double cosine=cos(rad);
	double sine  =sin(rad);
	double newy= y*cosine+z*sine;
	double newz=-y*sine+z*cosine;
	y=newy; z=newz;
}; /* vector3::rxcw(double rad) */
	
void vector3::rycw(double rad)
{
	double cosine=cos(rad);
	double sine  =sin(rad);
	double newx=x*cosine-z*sine;
	double newz=x*sine+z*cosine;
	x=newx; z=newz;
}; /* vector3::rycw(double rad) */

void vector3::rzcw(double rad)
{
	double cosine=cos(rad);
	double sine  =sin(rad);
	double newx= x*cosine+y*sine;
	double newy=-x*sine+y*cosine;
	x=newx; y=newy;
}; /* vector3::rzcw(double rad) */



void vector3::rotate(int axis, double degrees)
throw (invalid_argument)
{
	double rad=degrees*DEGREES;
	if (axis<0)
	{
		axis*=-1;
		rad*=-1.;
	}
	switch(axis)
	{
	case 1:
		rxcw(rad);
		break;
	case 2:
		rycw(rad);
		break;
	case 3:
		rzcw(rad);
		break;
	default:
		throw(invalid_argument("Hey, this is 3D-space! You tried rotating about an invalid axis."));
	}
	return;
}; /* vector3::rotate(int axis, double degrees) */
		

vector3& vector3::operator=(const vector3& rhs)
 throw ()
{
	if(this==&rhs)
		return *this;
	x=rhs.x;
	y=rhs.y;
	z=rhs.z;
	return *this;
};


bool vector3::_isNull() const throw ()
	{ return x==0. && y==0. && z==0.;};


bool vector3::operator== (const int& null) const
throw (invalid_argument)
{
	if (null==0)
		return _isNull();
	else
		throw(invalid_argument("vector3==int is only valid for int=0!!!"));
};

bool vector3::operator ==(const vector3& rhs) const throw ()
{
	if(this==&rhs)
		return true;
	else
		return (x==rhs.x) && (y==rhs.y) && (z==rhs.z);
};

//vector3::~vector3()
//{
//};


void vector3::setZ(const vector3& newZAxis) throw()
{
	if (&newZAxis==this)
	{
		z=modulus();
		x=y=0;
		return;
	};
	double axisY,phi;
	if (newZAxis.x != 0)
	{
		phi=atan2(newZAxis.y, newZAxis.x);
		rzcw(phi-PI/2.);
		// rotate ccw by PI/2 - phi -> sets x to 0
		axisY=cos(phi)*newZAxis.x + sin(phi)*newZAxis.y;
	}
	else
		axisY=newZAxis.y;

	if (axisY != 0)
	{
		phi=atan2(newZAxis.z, axisY);
		rxcw(phi-PI/2.);
	}
}; // vector3::setZ

