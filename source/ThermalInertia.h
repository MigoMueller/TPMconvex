// abstract base class for solving the one-dimensional problem of heat conduction

// Author: michael.mueller@dlr.de
// Date:   2004 Aug 30

// Assumptions: du/dt = d^2u / dz^2 (i.e.: natural units - cf. Spencer 1990)
// Boundary conditions: du/dz(z=infty) = 0
//                      u(z=0)^4       = f(t) + thermal parameter * du/dz(z=0)

// implicit assumptions: 
//       thermal conductivity, density, and heat capacity are constant with depth and temperature
//       skin depth is small compared to surface details like craters and the like
//       the annual heat wave can be neglected (day << year -- this would be wrong for Mercury!!!)


#if !defined(AFX_ThermalInertia_H__33B49FB6_9021_454C_A73B_1FC680D53233__INCLUDED_)
#define AFX_ThermalInertia_H__33B49FB6_9021_454C_A73B_1FC680D53233__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<vector>
#include<stdexcept>
#include "vectorN.h"

// BAUSTELLE: add member functions to alter / read the constant paramters!

class ThermalInertia  
{
protected:
	class TimeStep;       // defined below
	friend class TimeStep; 
protected:	
	const double parameter; // thermal parameter as defined by Spencer 1990
	const int nTime;        // # time steps (time = omega * physical time)
	const int nZ;           // # z steps    (z in units of skin depths)
	const double zMax;      // 'infinity' -- where constant temperature is assumed (~4 should be fine)
	const double accuracyGoal; 

	const double dT;        // 2 Pi / nTime
	const double dZ;        // zMax / (nZ-1)
	const double dZ2;       // dZ*dZ
	const double dT_dZ2;    // dT/dZ2 -- needed in solving the PDE; must be smaller than 0.5!!!
	const double pdZ;       // parameter * dZ -- needed for boundary condition

	mutable std::vector<TimeStep> cycle;       // temperature profile over one rotation period
	mutable vectorN              inRadiation; // incoming solar flux over one rotation period

public:
	virtual ~ThermalInertia(){};
	ThermalInertia(double parameter, int nTime, int nZ, double zMax, double accuracyGoal);
	int getNTime() const
		{return nTime;};
	double getDT() const throw()
		{return dT;};
	const vectorN& getInRadiation() const throw()
		{return inRadiation;};
	void setAccuracy(double accuracyGoal) const
		{const_cast<double&>(this->accuracyGoal) = accuracyGoal;}
private:
	ThermalInertia();
	ThermalInertia(const ThermalInertia&);
	//ThermalInertia operator= (const ThermalInertia&);
}; // class ThermalInertia



class ThermalInertia::TimeStep
{
protected:
	const ThermalInertia& conductivity;
	std::vector<double> profile;
public:
	// inRadiation will be clipped to be >= 0
	void calculate (const TimeStep& old, double inRadiation);
	double operator()() const
		{return profile[0];};
	double& operator()()
		{return profile[0];};
	TimeStep(const ThermalInertia& conductivity)
		: conductivity(conductivity),
		profile(conductivity.nZ, 0)
	{};
private:
	TimeStep();
public:
	TimeStep(const TimeStep& rhs)
		: conductivity(rhs.conductivity),
		profile(conductivity.nZ, 0)
	{
		for (int i=0; i<conductivity.nZ; i++)
			profile[i]=rhs.profile[i];
	};
	TimeStep& operator=(const TimeStep& rhs);
	TimeStep& operator=(double constant);
}; // class ThermalInertia::TimeStep





#endif // !defined(AFX_ThermalInertia_H__33B49FB6_9021_454C_A73B_1FC680D53233__INCLUDED_)
