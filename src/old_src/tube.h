#ifndef __TUBE_H__
#define __TUBE_H__

#include <list>
#include "btBulletDynamicsCommon.h"
#include <memory>
#include "rectangle.h"
#include <vector>

using namespace std;

// Allows for multiple return values for the fixed point for tube separation calculation
struct tubeSepResult
{
	int iter;
	bool converge;
	long double result;

	tubeSepResult()
	{
		iter = 0;
		converge = false;
		result = 0; //converged value, true/real answer
	}
};

//Allows for multiple return values for fixed point for height calculation
struct heightResult
{
	int iter;
	bool converge;
	long double result;

	heightResult()
	{
		iter = 0;
		converge = false;
		result = 0;
	}
};

//relevant parameters for tube construction
struct tubeParams
{
	long double height;
	long double tubeSeparation;
	int numSections;

	tubeParams()
	{
		height = 0;
		tubeSeparation = 0;
		numSections = 0;
	}
};

class tube
{
	int tubeNum;
	double length;
	double cylHeight;
	double tubeSpacing;
	double minSpacing;
	double diameter;
	double curvature;
	int numSections;
	bool finished;
	shared_ptr<rectangle> tubeAabb;
	shared_ptr<vector<btRigidBody*>> cylList;


private:
    void getCylHeight(const double heightGuess, const int numSections, heightResult &hres);
	void getTubeSeparation(const double height, const double guess, tubeSepResult &tsres);
	void extractTubeParams(tubeParams &res);


	/*There are special ways that I want to handle this object as the btRigidBody 
	pointers in the cylList will be copies of the pointers in the dynamics world.
	I decided that it would be easiest to use all pointers and not try to use
	smart pointers as I would have to change the underlying structure of the dynamics
	engine and I do not want to do that. CNT object can use the shared pointer though
	as that is my object. 
	
	For the cylList, I need to make sure that the only remaining shared pointer is the
	one left in the tube object, so when the tube object is destroyed, I can ensure that
	all of the btRigidBody pointers are null so we don't have double deleting occuring.
	This should be the only special treatment of this object*/
public:
	tube(int num, double newDiameter, double newCurvature, double newLength, double newMinSpacing);
	~tube();
	void addCyl(btRigidBody* body);
	shared_ptr<vector<btRigidBody*>> getCylList();
	int getTubeNum();
	double getLength();
	double getCylHeight();
	double getTubeSpacing();
	double getMinSpacing();
	double getDiameter();
	double getCurvature();
	double getNumSections();
	void setFinished();
	bool isFinished();
	shared_ptr<rectangle> getRectangle();

};

#endif
