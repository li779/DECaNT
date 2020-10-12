#ifndef constants_h
#define constants_h

#include <cmath> // using pow

namespace constants
{
	// constants();

	static const double pi=3.141592;
	static const double inv_pi=1./pi;

	static const double eV=1.6*std::pow(10,-19.0); //[Joules]
	static const double hb=6.5*std::pow(10,-16.0)*eV; //[Joules.s]
	static const double kb=1.3865*std::pow(10,-23.0); //[Joules/Kelvin]

	static const double eps0 = 8.85*std::pow(10,-12); // permittivity of free space
	static const double q0 = 1.6*std::pow(10,-19); // charge of electron
}

#endif // constants_h
