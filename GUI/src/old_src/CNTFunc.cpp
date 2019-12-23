/**
CNTFunc.cpp
Purpose: Provides function implementation for CNTFunc.h
@author Alex Gabourie
*/

#include <math.h>

#include "btBulletCollisionCommon.h"

#include "CNTFunc.h"


/**
Calculates the translation vector and its length
*/
void CNT::setTvecLen()
{
	double yo = floor(-1.5);
	double t1 = (2 * m + n) / float(dr);
	double t2 = -(2 * n + m) / float(dr);
	TvecLen = sqrt(pow(sqrt(3.) / 2 * A_CC*(t1 + t2), 2) + pow(A_CC / 2 * (t1 - t2), 2));

}

/**
Gets the translation vector length

@return The translation vector length
*/
double CNT::getTvecLen() const
{
	return TvecLen;
}

/**
Calculates the number of atoms in the CNT unit cell and sets the class member for numAtomsUnitCell
*/
void CNT::setNumAtomsUnitCell()
{
	setdr();
	numAtomsUnitCell = static_cast<int>( 2 * (2 * pow(diameter * SIMD_PI, 2)) / (A_CC*A_CC*dr));

}

/**
Gets the number of atoms in the CNT unit cell

@return The number of atoms in the CNT unit cell
*/
int CNT::getNumAtomsUnitCell() const
{
	return numAtomsUnitCell;
}


/**
Calculates the greatest common divisor of two integers

@param m A number
@param n Another number
@return The greatest common divisor of m and n
*/
int CNT::gcd(int m, int n)     	// function definition
{                         	// block begin
	int  r;                	// declaration of remainder

	while (n != 0) {       	// not equal
		r = m % n;          	// modulus operator
		m = n;              	// assignment
		n = r;
	}                      	// end while loop
	return m;              	// exit gcd with value m
}

/**
sets greatest common denominator of (2m+n) and (2n+m)
*/
void CNT::setdr()
{
	int d = gcd(m, n);
	if ((n - m) % (3 * d) == 0)
	{
		dr = 3 * d;
	} else
	{
		dr = d;
	}
}

/**
calculates cnt helicity [0,pi/6]
*/
void CNT::setHelicity()
{
	helicity = acos(2 * (n + m*.5) / (sqrt(3 * pow((m + n)*1., 2) + pow((m - n)*1., 2))));
}

/**
calculates the diameter of the cnt [Angstroms]
*/
void CNT::setDiameter()
{
	diameter = A_CC*sqrt(pow(n*1., 2) + pow(m*1., 2) + n*m) / SIMD_PI;
}

/**
sets maximum curvature [radians/angstroms]
*/
void CNT::setCurvature()
{
	curvature = (1.49 / pow(diameter, 2))*(1 + 9.89 / pow(diameter, 5) * 1000.0 * cos(6.0 * helicity));

	/*Curvature method may need to be re-evaluated as our diameters are in the range of 14-20 
	angstroms and the equation is fit for 10-15 angstroms*/
}

/**
cnt constructor

@param hamada_n The n value for the chiral vector
@param hamada_m The m value for the chiral vector
*/
CNT::CNT(int hamada_n, int hamada_m)
{
	m = hamada_m;
	n = hamada_n;
	setHelicity();
	setDiameter();
	setCurvature();
	setNumAtomsUnitCell();
	setTvecLen();
}
	
//! cnt destructor : Implicitly called, error if try to call it

/**
gets nanotube helicity

@return The helicity of the nanotube
*/
double CNT::getHelicity() const
{
	return helicity;
}

/**
gets nanotube diameter

@return The diameter of the nanotube
*/
double CNT::getDiameter() const
{
	return diameter;
}

/**
gets nanotube curvature

@return The maximum curvature of the nanotube
*/
double CNT::getCurvature() const
{
	return curvature;
}

/**
gets m

@return m
*/
int CNT::getm() const
{
	return m;
}

/**
gets n

@return n
*/
int CNT::getn() const
{
	return  n;
}