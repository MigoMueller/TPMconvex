// hemispheric craters without thermal inertia
// Closely following Lagerros' analytical solution

// Author: michael.mueller@dlr.de
// Date:   2004 June 7


#if !defined(AFX_HEMISPHERICNOINERTIA_H__109EB529_5460_41F3_910F_B88FF3264775__INCLUDED_)
#define AFX_HEMISPHERICNOINERTIA_H__109EB529_5460_41F3_910F_B88FF3264775__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "TriangulatedConvex.h"
#include "HemisphericCrater.h"


class HemisphericNoInertia : 
	public HemisphericCrater, 
	public TriangulatedConvex::ThermalModelConvexNoInertia  
{
public:
    class HemisphericNoInertiaLagerros; // defined below

/*	double CorrectionFactor (double cosSun, double cosEarth, double cosAzi,
							 double lambdaMu, double rAU) const
	{
		double value=fluxModFactorsPerArea(cosSun, cosEarth, cosAzi, lambdaMu, rAU);
		if (value==0)
			return 0;
		return value/(cosEarth*HemisphericCrater::planck(cosSun));
	};
*/
	HemisphericNoInertia(double emissivity, double uncorrected_bond, double gammaDeg);
	HemisphericNoInertia(const TriangulatedConvex& shape, double gammaDeg);
	virtual ~HemisphericNoInertia()throw(){};
protected:
	double trapezTheta (double oldTrapez, double& start, double& stepSize, unsigned long& nPoints) const;
	double trapezPhi (unsigned long& nPoints, double oldTrapez, double& cosStep, double& sinStep,
									   double cosTheta, double sinTheta) const;
	double phiIntegral (double cosTheta) const;
	virtual double fluxModFactorsPerArea (double cosSun, double cosEarth,
										  double cosAzi, double lambdaMu, double rAU) const;
	virtual inline double planck (const vector3& unitNormal) const;

	const double TempTerm0;
	mutable double TempTerm;
	mutable double planck0;
private:
	// maximum #iterations for trapez-integrals
	static const int maxTrapez;
	inline double integrand (const vector3& unitNormal) const;
	inline void init();
};

// Lagerros is like the 'real' thing, but multipleReflection0 is re-defined in hack()
class HemisphericNoInertia::HemisphericNoInertiaLagerros : public HemisphericNoInertia
{
public:
    HemisphericNoInertiaLagerros(double emissivity, double uncorrected_bond, double gammaDeg)
        : HemisphericNoInertia(emissivity, uncorrected_bond, gammaDeg)
        {hack();};
    HemisphericNoInertiaLagerros(const TriangulatedConvex& shape, double gammaDeg)
        : HemisphericNoInertia(shape, gammaDeg)
        {hack();};
    virtual ~HemisphericNoInertiaLagerros()throw(){};
private:
    void hack()
    {
        const_cast<double&> (this->multipleReflection0) = (1-emissivity)*S*(1-S);
    }
}; // class HemisphericNoInertiaLagerros



#endif // !defined(AFX_HEMISPHERICNOINERTIA_H__109EB529_5460_41F3_910F_B88FF3264775__INCLUDED_)
