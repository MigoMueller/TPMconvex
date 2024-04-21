// class vector3
// by michael.mueller@dlr.de
// 2003 Oct. 15
// Two ways to construct vectors:
// vector3(x,y,z) -- euclidean vector
// vector3(longitude, latitude) in degrees --> unit vector (x=cos(longitude)cos(latitude) and so on...)
// (Nearly) all reasonable vector-operations are defined,
// (nearly) all reasonable operators are adequately overloaded.

// One add'l ctor:
// vector3(vector1, vector2) is cross-product of the two
// More concise, however: vector-valued fct. cross(vec1, vec2)!

// body file: vector3.cpp

// Not supposed to be a base-class for anything, non-virtual dtor!!!

#if !defined(AFX_vector3_H__B2ED94A8_E1DA_4371_97BB_302E23983C18__INCLUDED_)
#define AFX_vector3_H__B2ED94A8_E1DA_4371_97BB_302E23983C18__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<cmath>
#include<stdexcept>


class vector3  
{
public:
	double x,y,z;
	static const double PI;
	static const double PI2; // =2*PI
	static const double DEGREES;//=PI/180.;
	static const double RAD;//=180./PI;

	// Euclidean vector
	vector3(double x, double y, double z) throw()
		: x(x), y(y), z(z) {};

	// Polar unit vector
	vector3(double longitude, // in degrees!
		   double latitude   // in degrees!
		   ) throw();
	
	// unusual ctor, yields cross-product of lhs and rhs (in this order)
	vector3 (const vector3& lhs, const vector3& rhs) throw ();

	vector3& operator+= (const vector3& rhs) throw()
	{
		x+=rhs.x; y+=rhs.y; z+=rhs.z;
		return *this;
	};

	vector3& operator-= (const vector3& rhs) throw()
	{
		x-=rhs.x; y-=rhs.y; z-=rhs.z;
		return *this;
	};

	vector3& operator*= (const double &factor) throw()
	{
		x*=factor; y*=factor; z*=factor;
		return *this;
	};

	vector3& operator/= (const double &factor) 
		throw (std::invalid_argument)
	{
		if (factor==0)
			throw(std::invalid_argument("Division by zero"));
		x/=factor; y/=factor; z/=factor;
		return *this;
	};

	double dot (const vector3& rhs) const throw()
		{return x*rhs.x + y*rhs.y + z*rhs.z;};

	double squared() const throw()
		{return (x*x+y*y+z*z);};

	double modulus() const throw()
		{return sqrt(squared());};

	// rotates vector cw around x-axis
	// argument in rad!
	void rxcw(double rad);
	// rotates vector cw around y-axis
	// argument in rad!
	void rycw(double rad);
	// rotates vector cw around z-axis
	// argument in rad!
	void rzcw(double rad);

	// axis: 1=x, 2=y, 3=z, degrees>0: cw rotation
	// (special feature: axis<0 means ccw rotation, equivalently degrees<0...)
	// throws std::invalid_argument for axis==0 or |axis|>3
	void rotate(int axis, double degrees)
		throw (std::invalid_argument);

	// rotates vector into a coordinate system in which newZAxis is z-axis
	// (order: first about z-axis, then about x-axis)
	void setZ(const vector3& newZAxis) throw();

	vector3& operator=(const vector3& rhs)
		 throw ();

	vector3 operator-() const throw ()
		{ return vector3(-x, -y, -z);};

	bool operator==(const vector3& rhs) const throw ();
	bool operator!=(const vector3& rhs) const throw ()
		{ return !(*this==rhs);};
	bool _isNull() const throw ();
	bool operator== (const int& null) const
		throw (std::invalid_argument);
	~vector3() throw() {};
	
	static double dot(const vector3& lhs, const vector3& rhs)
		{return lhs.dot(rhs);};
	// angle btw rhs and lhs in degrees
	static double angle (const vector3& rhs, const vector3& lhs)
	{ 
		double cosine=dot(rhs, lhs)/rhs.modulus()/lhs.modulus();
		if(cosine>=1.)
			return 0;
		if(cosine<=-1.)
			return 180;
		return acos(cosine)*vector3::RAD;
	};
	// cross product
	static const vector3 cross(const vector3& lhs, const vector3& rhs)
	{ return vector3(lhs, rhs); };
private:
	vector3() {}; // prohibitive
}; // class vector3


// Overload arithmetric operators for vectors
inline const vector3 operator+(const vector3& lhs, const vector3& rhs)
	{ return vector3(lhs)+=rhs;  };
inline const vector3 operator-(const vector3& lhs, const vector3& rhs)
	{ return vector3(lhs)-=rhs;  };
inline const vector3 operator*(const vector3& vec, const double &factor)
	{ return vector3(vec)*=factor; };
inline const vector3 operator* (const double factor, const vector3& vec)
	{ return vec*factor; };
// make multiplication commutative
inline const vector3 operator/ (const vector3& vec, double factor)
throw (std::invalid_argument)
	{ return vector3(vec)/=factor; };



#endif // !defined(AFX_vector3_H__B2ED94A8_E1DA_4371_97BB_302E23983C18__INCLUDED_)
