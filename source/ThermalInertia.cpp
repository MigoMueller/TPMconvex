#include "ThermalInertia.h"
#include<cmath>
#include<sstream>
#include<iostream>
using namespace std;


ThermalInertia::ThermalInertia(double parameter, int nTime, int nZ, double zMax, double accuracyGoal)
	: parameter(parameter),
	  nTime(nTime),
	  nZ(nZ),
	  zMax(zMax),
	  accuracyGoal(accuracyGoal),
	  dT(8*atan(1.)/nTime),  // tan(pi/4)=1
	  dZ(zMax / (nZ-1)), 
	  dZ2(dZ*dZ),
	  dT_dZ2(dT/dZ2),
	  pdZ (parameter/dZ),
	  inRadiation(nTime)
{
	cycle.resize(nTime, TimeStep(*this));	
	if (dT_dZ2 > 0.5)
	{
		ostringstream dummy;
		dummy<<"ThermalInertia::ctor: dT / dZ^2 is "<<dT_dZ2<<", it must be <= 0.5 for numerical stability!\n";
		dummy<<"Increase the #time-steps!"<<ends;
		throw std::invalid_argument(dummy.str());
	};
}; // BAUSTELLE: add some checking of positiveness and the like!



void ThermalInertia::TimeStep::calculate(const TimeStep &old, double inRadiation)
{
	if (inRadiation<0)
		inRadiation = 0;
	const int   & nZ  = conductivity.nZ;
	//const double& dT  = conductivity.dT;
	//const double& dZ  = conductivity.dZ;
	//const double& dZ2 = conductivity.dZ2;
	//const double& parameter = conductivity.parameter;
	const double& dT_dZ2 = conductivity.dT_dZ2;
	const double& pdz    = conductivity.pdZ;

	// determine inner points of profile according to the PDE
	{
		for (int i=1; i<=nZ-2; i++)
		{
			profile[i] = old.profile[i]+ dT_dZ2*(old.profile[i-1] + old.profile[i+1] - 2*old.profile[i]);
			if (profile[i] < 0)
			{
				cerr<<"Encountered T<0 below surface; continue anyway with T=0."<<endl;
				profile[i]=0;
			};
			if (profile[i] > 1)
			{
				cerr<<"Encountered T>TSS below surface; continue anyway with T=TSS."<<endl;
				profile[i]=1;
			}
		};
	};
	// determine lowest point according to the PDE + boundary condition:
	// think of a virtual profile[nZ] == profile[nZ-1] due to BC (d/dz = 0)
	profile[nZ-1] = old.profile[nZ-1] + dT_dZ2*(old.profile[nZ-2] - old.profile[nZ-1]);
	if (profile[nZ-1] < 0)
	{
		cerr<<"Encountered T<0 below surface; continue anyway with T=0."<<endl;
		profile[nZ-1]=0;
	};
	if (profile[nZ-1] > 1)
	{
		cerr<<"Encountered T>TSS below surface; continue anyway with T=TSS."<<endl;
		profile[nZ-1]=1;
	}

	// Solve boundary condition @ surface using Newton's method:
	// f(u) = u^4 - inRadiation + parameter/dZ (u-u1) == 0	
	const double u1  = profile[1];
	//const double pdz = parameter/dZ;
	// use 'old' surface temperature as first guess:
	double u = old.profile[0];
	double deltaU, u3;
	double oldU, fu, fpu;
	do
	{
		u3 = u*u*u;
		fu  = u3*u - inRadiation + pdz*(u-u1); // f(u)
		fpu = 4*u3 + pdz;					  // df/du
		deltaU = -fu/fpu;        // f(u+deltaU) = 0 to first order
		oldU = u;
		u+=deltaU;
		if (u<0)
		{
			u = oldU/2;
			deltaU=oldU; // make sure we don't stop after this
		}
		if (u>1)
		{
			u = (1+oldU)/2;
			deltaU=oldU; // make sure we don't stop after this
		};
	}
	while ( fabs(deltaU/oldU) > conductivity.accuracyGoal);
	profile[0] = u;
}; //ThermalInertia::TimeStep::calculate



ThermalInertia::TimeStep& ThermalInertia::TimeStep::operator =(const ThermalInertia::TimeStep &rhs)
{
	if (&conductivity != &rhs.conductivity)
		throw invalid_argument("ThermalInertia::Timestep::operator=: lhs and rhs must belong to same ThermalInertia!");
	for (int i=0; i<conductivity.nZ; i++)
		profile[i] = rhs.profile[i];
	return *this;
}; //ThermalInertia::TimeStep::operator =



ThermalInertia::TimeStep& ThermalInertia::TimeStep::operator =(double constant)
{
	if (constant<0)
	{
		ostringstream dummy;
		dummy<<"ThermalInertia::Timestep::operator=: number to assign to must be non-negative; you tried "<<constant<<ends;
		throw invalid_argument(dummy.str());
	};
	for (int i=0; i<conductivity.nZ; i++)
		profile[i]=constant;
	return *this;
};
