#include "tube.h"
#include <iostream>
#include "rectangle.h"

using namespace std;

/**
Creates the structure responsible for remembering all bullet physics related information
for each nanotube

@param num Represents the numth tube created in the simulation. Bookkeeping number
@param newDiameter The diameter of the tube
@param curvature The max curvature for the nanotube bundle
@param newLength The length of the nanotube
@param newMinSpacing The distance the spacing cylinder goes out from the tube diameter
*/
tube::tube(int num, double newDiameter, double newCurvature, double newLength, double newMinSpacing)
{
	tubeNum = num;
	length = newLength;
	minSpacing = newMinSpacing;
	curvature = newCurvature;
	diameter = newDiameter;
	cylList = static_cast<shared_ptr<vector<btRigidBody*>>>(new vector<btRigidBody*>());
	double a_min = 40.0; //minumum lenght of the bundle segments [nanometers]
	tubeSpacing = 0.1;
	numSections = static_cast<int>(floor(length / a_min));
	cylHeight = (length - tubeSpacing*(numSections-1)) / numSections;
	tubeAabb = make_shared<rectangle>(rectangle());
	tubeAabb->setTubeNum(num);
	finished = false;
}


/**
Deletes pointers for the cylinder list such that they are not double deleted later. (Sketchy, I know)
*/
tube::~tube()
{
	// Need to go assign all of the btRigidBody objects to null as they are already deleted 
	//  elsewhere. This may look like a memory leak because by setting the pointers to null,
	// the actual pointer cannot be deleted. However, there are two pointers pointing to the
	// memory locations in question and the second set of pointers handles the memory
	for (auto itr = cylList->begin(); itr != cylList->end(); ++itr)
	{
		*itr = nullptr;
	}

	/* cnt and tubeNum should delete fine on their own. The list can now clean itself up.*/
}

/**
Gets the number of the tube

@return The tube's creation number
*/
int tube::getTubeNum()
{
	return tubeNum;
}

/**
Adds a pointer of the cylinder to the end of the cylList

@param body The pointer to the cylinder to be added
*/
void tube::addCyl(btRigidBody* body)
{
	cylList->push_back(body);
}

/**
Gets the entire list of cylinders

@return A pointer to a list of rigid body object pointers
*/
shared_ptr<vector<btRigidBody*>> tube::getCylList()
{
	return cylList;
}

/**
Gets the length of the tube

@return the length of the nanotube
*/
double tube::getLength()
{
	return length;
}

/**
Gets the height of the cylinders in the tube (all the same)

@return the cylinder height
*/
double tube::getCylHeight()
{
	return cylHeight;
}

/**
Gets the spacing between other nanotubes

@return the minimum spacing between this tube and others
*/
double tube::getMinSpacing()
{
	return minSpacing;
}

/**
Gets the spacing between adjacent cylinders within the same nanotube

@return the spacing between adjacent cylinders within the same nanotube
*/
double tube::getTubeSpacing()
{
	return tubeSpacing;
}


double tube::getDiameter()
{
	return diameter;
}

/**
Calculates all of the necessary parameters for the simulation to build correct CNT bundles

@param tubeParams A struct that holds relevant tube building parameters. (cylHeight, cylSpacing, etc)
*/
void tube::extractTubeParams(tubeParams& res)
{
	double a_min = 20.0; //nanometers
	double a = a_min;
	double range = a_min*1.5;
	if (length < a_min)
	{
		throw new runtime_error(string("Desired length must be greater than 20 nm.\n"));
	}

	double steps = 1000;
	double stepSize = range / steps;

	double prevDec = 1;
	double prevVal = 0;
	int numSections = 0;
	//t_a increases monotonically with a, so we will not have any unexpected results

	tubeSepResult tubeSepRes = tubeSepResult();
	heightResult heightInfo = heightResult();

	//finds number of tubes at an approximate height
	while (a < a_min + range)
	{
		tubeSepRes.converge = false;
		getTubeSeparation(a, btScalar(0.), tubeSepRes);
		//checks to make sure calculation went alright
		if (!tubeSepRes.converge)
		{
			throw new runtime_error(string("Tube separation calculation did not converge.\n"));
		}
		double t_a = tubeSepRes.result;

		double val = (length + t_a) / (a + t_a);
		double currDec = val - floor(val);
		if (currDec > prevDec)
		{
			numSections = static_cast<int>(floor(prevVal));
			break;
		}
		prevDec = currDec;
		prevVal = val;
		a += stepSize;
	}

	
	try
	{
		getCylHeight(a, numSections, heightInfo);
	}
	catch (runtime_error err)
	{
		string error = err.what();
		throw new runtime_error(error);
	}
	//check if height calc went alright -> This is currently a double check
	if (!heightInfo.converge)
	{
		throw new runtime_error(string("Section height calculation did not converge.\n"));
	}
	res.height = heightInfo.result;

	tubeSepRes.converge = false;
	getTubeSeparation(res.height, btScalar(0.), tubeSepRes);
	//checks to make sure calculation went alright
	if (!tubeSepRes.converge)
	{
		throw new runtime_error(string("Tube separation calculation did not converge.\n"));
	}
	res.tubeSeparation = tubeSepRes.result;

	res.numSections = numSections;
}

/**
Calculates the required cylinder height based on multiple parameters

@param heightGuess The guess for newton's fixed point method
@param numSections The number of sections to be used to make the length of the tube correct
@param hres The object that contains the height of each cylinder and some other method information
*/
void tube::getCylHeight(const double heightGuess, const int numSections, heightResult& hres)
{
	double tol = .0000000001;
	int iterlim = 500;
	double a = heightGuess;
	double h = .001;

	tubeSepResult t = tubeSepResult();
	tubeSepResult tplus = tubeSepResult();
	tubeSepResult tminus = tubeSepResult();

	double separationGuess;
	for (int i = 0; i < iterlim && !hres.converge; i++)
	{
		t.converge = false;
		tplus.converge = false;
		tminus.converge = false;
		separationGuess = tan(curvature*a / 2.0)*diameter;
		getTubeSeparation(a, separationGuess, t);
		getTubeSeparation(a + h, separationGuess, tplus);
		getTubeSeparation(a - h, separationGuess, tminus);
		//checks to make sure calculation went alright
		if (!(t.converge && tplus.converge && tminus.converge))
		{
			throw new runtime_error(string("Tube separation calculation did not converge.\n"));
		}
		double t_a = t.result;
		double t_ap = tplus.result;
		double t_am = tminus.result;

		double faplus = (length + t_ap) / ((a + h) + t_ap) - float(numSections);
		double faminus = (length + t_am) / ((a - h) + t_am) - float(numSections);
		double fprime = (faplus - faminus) / (2 * h);
		double anew = a - ((length + t_a) / (a + t_a) - float(numSections)) / fprime;

		if (abs(anew - a) < tol)
		{
			hres.converge = true;
			hres.result = anew;
			hres.iter = i;
		}
		a = anew;
	}
}

/**
Gets the necessary separation between two cylinders to get the correct length


@param height Height of the cylinder
@param guess The guess for newton's fixed point method
@param tsres The tube separation number and some other call specific stats 
*/
void tube::getTubeSeparation(const double height, const double guess, tubeSepResult& tsres)
{
	double tol = .0000000001;
	int iterlim = 500;

	double x = guess;

	double radius = diameter / 2.0;

	for (int i = 0; i < iterlim && !tsres.converge; i++)
	{
		double xnew = x - (atan(x / (2.0*radius)) / curvature - x / 2 - height) /
			((1 / curvature)*(radius*2.0) / (pow(x, 2) + 4 * pow(radius, 2)) - .5);
		if (abs(xnew - x) < tol)
		{
			tsres.converge = true;
			tsres.iter = i;
			tsres.result = xnew;
		}
		x = xnew;
	}
}


/**
Gets the curvature max curvature of the tube
*/
double tube::getCurvature()
{
	return curvature;
}

/**
Gets the number of sections in the tube
*/
double tube::getNumSections()
{
	return numSections;
}

//checks if the tube is completed building
bool tube::isFinished()
{
	return finished;
}

//tube is finished building and aabb rectangle can be built
void tube::setFinished()
{
	if (finished == true)
	{
		return;
	}

	finished = true;

	auto firstCyl = cylList->front();
	btVector3 firstAabbMin, firstAabbMax;
	firstCyl->getAabb(firstAabbMin, firstAabbMax);
	btVector3 firstPos = firstCyl->getCenterOfMassPosition();

	auto lastCyl = cylList->back();
	btVector3 lastAabbMin, lastAabbMax;
	lastCyl->getAabb(lastAabbMin, lastAabbMax);
	btVector3 lastPos = lastCyl->getCenterOfMassPosition();

	double xMax;
	double xMin;
	double zMax;
	double zMin;

	if (firstPos.getX() > lastPos.getX())
	{
		xMax = firstAabbMax.getX();
		xMin = lastAabbMin.getX();
	}
	else
	{
		xMax = lastAabbMax.getX();
		xMin = firstAabbMin.getX();
	}

	if (firstPos.getZ() > lastPos.getZ())
	{
		zMax = firstAabbMax.getZ();
		zMin = lastAabbMin.getZ();
	}
	else
	{
		zMax = lastAabbMax.getZ();
		zMin = firstAabbMin.getZ();
	}

	tubeAabb->setHeight(zMax - zMin);
	tubeAabb->setWidth(xMax - xMin);

}

//Gets the Aabb rectangle for entire straight tube
shared_ptr<rectangle> tube::getRectangle()
{
	return tubeAabb;
}

