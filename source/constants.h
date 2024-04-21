#if !defined(AsteroidConstants_H_INCLUDED)
#define AsteroidConstants_H_INCLUDED

namespace AsteroidConstants
{
	// some natural constants in (sort of) SI-units
	const double h=6.6260755e-34;  // Js
	const double c=2.99792458e8;   // m/s
	const double k=1.380658e-23;   // J/K
	const double solar=1367;       // W/m^2 // Thuillier et al. (2004)
	const double sigma=5.67051e-8;    // W/m^2/K^4
	const double AU=149597870.691; // AU->km
	const double degrees=3.14159265358979/180.;  //degrees->rad
	
	// derived constants in practical units
	const double solar_over_stefan_boltzmann=solar/sigma; // in K^4
	const double hck=1e6*h*c/k;    // mu*K (value: 14387.686603333910352889709109714)
	const double hc2=1e24*h*c*c;   // in W/(m^2 mu)*(mu)^5
	// note: will be divided by lambda^5 (in mu)
	
	// Asteroid-specific
	inline double q(double G){ return (.29+.684*G);}// phase integral
} // namespace AsteroidConstants

#endif // AsteroidConstants_H_INCLUDED
