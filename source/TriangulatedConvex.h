// Lean implementation of convex asteroids' shape models and thermal modelling on them.


#if !defined(AFX_TRIANGULATEDCONVEX_H__62813EF5_A003_4F6B_8F84_D569F494BB56__INCLUDED_)
#define AFX_TRIANGULATEDCONVEX_H__62813EF5_A003_4F6B_8F84_D569F494BB56__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "asteroid.h"
#include "ConvexFile.h"
#include "vectorN.h"
//#include "ThermalInertiaOnlyConvex.h"
#include<iosfwd> // forward declaration of, e.g., ostream (used in ThermalLightCurve)



class TriangulatedConvex : public asteroid  
{
public:
	class LambertianEmitter;           // dto., inherits from above
	class ThermalModelConvex;          // defined after end of TriangulatedConvex
	class ThermalModelConvexNoInertia; // dto., inherits from above

	TriangulatedConvex (const ConvexFile& file, const SpinState& spinAxis,
		double H, double G, double pV, double emissivity=0.9);

	double ThermalFlux (const vector3& A2Sun, const vector3& A2Earth, 
		double lambdaMu, double rAU, double deltaKM, const ThermalModelConvex& model) const;
	std::vector<double> ThermalFlux (const vector3& A2Sun, const vector3& A2Earth, 
		const std::vector<double>& lambdaMu, double rAU, double deltaKM, const ThermalModelConvex& model) const;
	double ThermalFlux (const asteroid::fluxRecord& record, const ThermalModelConvex& model) const
		{return ThermalFlux (record.Sun2Aster, record.Earth2Aster, record.JD, record.lambdaMu, model);};

	// template: will work with single lambda or with vector of lambdas
	template<class lambdaType>
	lambdaType ThermalFlux (const eclipticVector& Sun2Aster, 
		const eclipticVector& Earth2Aster, 
		double JD, const lambdaType& lambdaMu, const ThermalModelConvex& model) const
	{
		vector3 A2Sun   = bodyFix(Sun2Aster, JD);
		A2Sun *= -1; // unit vector pointing away from asteroid at zero time
		vector3 A2Earth = bodyFix(Earth2Aster, JD);
		A2Earth *= -1;	
		double rAU= Sun2Aster.getLength();
		double deltaKM=Earth2Aster.getLength()*AsteroidConstants::AU;
		return ThermalFlux (A2Sun, A2Earth, lambdaMu, rAU, deltaKM, model);
	}; // TriangulatedConvex::ThermalFlux

	//void ThermalLightCurve (const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster,
	//	double lambdaMu, double JD0, const ThermalModelConvex& model)
	
	template<class modelType>
	vector3 YarkovskyForce (const vector3& A2Sun, double rAU, const modelType& model) const
	{
		vector3 result (0,0,0);
		double TSS = model.getTSS(rAU);
		double factor = scaleFactor * 1e3 * TSS*TSS; // 1e3: convert diameter(km) -> diameter (m)
		factor *= factor * 2./3.* AsteroidConstants::sigma * emissivity / AsteroidConstants::c;
		// factor = 2/3 e sigma TSS^4 *scaleFactor^2 / c --> prefactor of (Lambertian) Yarkovsky force

		std::vector<vector3>::const_iterator it, ende=dA.end();
		for (it = dA.begin(); it != ende; it++)
			result += model.Yarkovsky_u4dA (A2Sun, *it);		
		return factor*result;
	};

	template<class modelType>
        std::vector<vectorN> getULightCurves(const eclipticVector& Sun2Aster, const eclipticVector& Obs2Aster, double JD, const modelType& model,
					       std::vector<vector3>& A2Suns, std::vector<vector3>& A2Obses, std::vector<double>& JDs) const
	{
	    std::vector<vectorN> result;
	    std::vector<vector3>::const_iterator it, ende=dA.end();
	    vector3 A2Sun   = bodyFix(Sun2Aster, JD);
	    A2Sun *= -1;
	    vector3 A2Obs = bodyFix(Obs2Aster, JD);
	    A2Obs *= -1;
	    for (it=dA.begin(); it!=ende; it++)
	      result.push_back(model.uLightCurve(A2Sun, *it));
	    const unsigned int nTime = model.getNTime();
	    A2Suns = std::vector<vector3>();
	    A2Obses = std::vector<vector3>();
	    JDs = std::vector<double>();
	    const double periodH = spinAxis.getTh();
	    const double dTDay = periodH/(24*nTime); // time increment in days
	    const double angularIncrement = model.getDT(); // angle increment in rad by which to rotate A2Sun and A2Obs
	    for (int i=0; i<nTime; i++)
	    {
	      JDs.push_back(JD);
	      A2Suns.push_back(A2Sun);
	      A2Obses.push_back(A2Obs);
	      JD += dTDay;
	      A2Sun.rzcw(angularIncrement);
	      A2Obs.rzcw(angularIncrement);
	    };
	    return result;
        };

	virtual ~TriangulatedConvex()throw(){};

	virtual inline void setPV(double newAlbedo) throw (std::invalid_argument)
		{asteroid::setPV(newAlbedo); setScaleFactor();};
	virtual inline void setH(double H)
		{asteroid::setH(H); setScaleFactor();};
protected:
	const std::vector<vector3> dA;
	const double intrinsicDiameter;
	double scaleFactor;  // = asteroid::diameter/intrinsicDiameter (km/model unit length)
	void setScaleFactor()
		{scaleFactor=diameter/intrinsicDiameter;};

public:
// class definitions:
	


	class ThermalModelConvex : public asteroid::ThermalModel
	{
	public:
		// returns the flux emitted by the facet described by dA up to constant factors:
		// 2 pi emissivity h c^2 / (delta^2 lambda^5) * scaleFactor^2 
		// These must be multiplied in calling function!
		virtual std::vector<double> fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
			const std::vector<double>& lambdaMu, double rAU) const 
			= 0;
		virtual vectorN ThermalLightCurveModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
			double lambdaMu, double rAU, unsigned int nPointsEach) const
			= 0;
		inline std::vector<double> CorrectionFactor (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
			const std::vector<double>& lambdaMu, double rAU) const; 
		virtual ~ThermalModelConvex()throw(){};
	protected:
		ThermalModelConvex(double emissivity, double bond)
			: ThermalModel(emissivity, bond){};
		ThermalModelConvex (const TriangulatedConvex& shape)
			: ThermalModel(shape){};
	}; // class TriangulatedConvex::ThermalModelConvex
	
	
	
	class ThermalModelConvexNoInertia : public ThermalModelConvex
	{
	public:
		// split up calculations for vector<lambda> into one per lambda
		inline virtual std::vector<double> fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
			const std::vector<double>& lambdaMu, double rAU) const;
		// calculate cosines of Earth, Sun, and their azimuth + have fluxModFactorsPerArea do the job
		double fluxModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
			double lambdaMu, double rAU) const;
		virtual vectorN ThermalLightCurveModFactors (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
			double lambdaMu, double rAU, unsigned int nPointsEach) const;

		~ThermalModelConvexNoInertia()throw(){};
	protected:
		// this is where the deriving class takes over
		virtual double fluxModFactorsPerArea (double cosSun, double cosEarth, double cosAzi,
			double lambdaMu, double rAU) const = 0;
		ThermalModelConvexNoInertia (double emissivity, double bond) : ThermalModelConvex(emissivity, bond) {};
		ThermalModelConvexNoInertia (const TriangulatedConvex& shape) : ThermalModelConvex(shape) {};
	}; // class TriangulatedConvex:ThermalModelConvexNoInertia

	
	
	class LambertianEmitter : public ThermalModelConvexNoInertia
	{
	public:
		LambertianEmitter(const TriangulatedConvex& shape) : ThermalModelConvexNoInertia(shape) {};
		LambertianEmitter(double emissivity, double bond)  : ThermalModelConvexNoInertia(emissivity, bond) {};
		virtual ~LambertianEmitter () throw() {};
	protected:
		virtual double fluxModFactorsPerArea (double cosSun, double cosEarth, double cosAzi,
			double lambdaMu, double rAU) const 
		{return cosEarth /( exp( AsteroidConstants::hck/(lambdaMu * TSS(rAU, cosSun))) - 1);};
	}; // class TriangulatedConvex::LambertianEmitter


	class NEATM : public ThermalModelConvexNoInertia
	{
	protected:
		const double eta;
	public:
		NEATM(const TriangulatedConvex& shape, double eta)
			: ThermalModelConvexNoInertia(shape), eta(eta)
		{};
		NEATM(double emissivity, double bond, double eta)
			: ThermalModelConvexNoInertia(emissivity, bond), eta(eta)
		{};
		virtual ~NEATM() throw() {};
	protected:
		virtual double fluxModFactorsPerArea (double cosSun, double cosEarth, double cosAzi,
			double lambdaMu, double rAU) const 
		{return cosEarth /( exp( AsteroidConstants::hck/(lambdaMu * TSS(rAU, cosSun,eta))) - 1);};

	}; // class TriangulatedConvex::NEATM

	// makes two revolutions with nPointsEach points each (only one will be calculated, however)
	// (if thermal inertia is involved: nTime points each - nPointsEach will simply be ignored)
	void ThermalLightCurve(std::ostream& out, const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster, 
		double lambdaMu, double JD0, const ThermalModelConvex& model, int nPointsEach) const;
	// as above, but mixing two models
	// (think of model1 as the one with craters, model without craters, and density1 = crater density)
	void ThermalLightCurve(std::ostream& out, const eclipticVector& Sun2Aster, const eclipticVector& Earth2Aster, 
		double lambdaMu, double JD0, 
		const ThermalModelConvex& model1, const ThermalModelConvex& model2, double density1, 
		int nPointsEach) const;

}; // class TriangulatedConvex



// !!!!!!!!!!!!!!!!!!!!!!!!!
// inline implementations
// !!!!!!!!!!!!!!!!!!!!!!!!!



inline std::vector<double> TriangulatedConvex::ThermalModelConvex::CorrectionFactor 
	(const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
	 const std::vector<double>& lambdaMu, double rAU) const
{
	LambertianEmitter Lambert(emissivity, bond);
	std::vector<double> result = fluxModFactors(A2Sun, A2Earth, dA, lambdaMu, rAU);
	std::vector<double> reference = Lambert.fluxModFactors (A2Sun, A2Earth, dA, lambdaMu, rAU);
	for (unsigned int i=0; i<lambdaMu.size(); i++)
	{
		if (reference[i] == 0)
			result[i] = 0;
		else
			result[i] = result[i]/reference[i];
	};
	return result;
};  // TriangulatedConvex::ThermalModelConvex::CorrectionFactor



inline std::vector<double> TriangulatedConvex::ThermalModelConvexNoInertia::fluxModFactors 
   (const vector3& A2Sun, const vector3& A2Earth, const vector3& dA,
    const std::vector<double>& lambdaMu, double rAU) const 
{
	const std::vector<double>::size_type N=lambdaMu.size();
	std::vector<double> result(N,0);
	for (unsigned int i=0; i<N; i++)
		result[i] = fluxModFactors(A2Sun, A2Earth, dA, lambdaMu[i], rAU);
	return result;
}; // TriangulatedConvex::ThermalModelConvexNoInertia::fluxModFactors 





#endif // !defined(AFX_TRIANGULATEDCONVEX_H__62813EF5_A003_4F6B_8F84_D569F494BB56__INCLUDED_)
