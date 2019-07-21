/**
MeshEnv.cpp
Purpose: Creates the world, objects, and camera manipulations for the CNT mesh demo

@author Alex Gabourie
@version 1.01 4/15/15

This is a heavily edited version of BasicDemo.cpp provided by Bullet Physics
*/

// Definitions:
// Spacing Cylinder - The large cylinder used to ensure a minimum distance between nanotubes.
//		Non-rendered.
// Spine Cylinder - The smaller cylinder that is used for curvature limitations and represents
//		the actual radius of the nanotubes. Rendered

#define SIMD_QUARTER_PI (SIMD_PI * btScalar(0.25)) 

///btBulletDynamicsCommon.h is the main Bullet include file, contains most common include files.

#include <stdio.h>
#include <limits>
#include <time.h>
#include <iostream>
#include <sstream>
#include <regex>
#include <experimental/filesystem>


// #include "GlutStuff.h"
#include "btBulletDynamicsCommon.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

#include "LinearMath/btIDebugDraw.h"
#include "LinearMath/btAabbUtil2.h"
#include "LinearMath/btMatrix3x3.h"
#include "LinearMath/btVector3.h"

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rectangle.h"
#include "MeshEnv.h"


using namespace rapidxml;
using namespace std;
// static GLDebugDrawer gDebugDraw;

// why using global variables?????
string outputPath;
bool toDebugDraw = true;
string runID;
bool manualEndSimulation = false;

bool compareRect(shared_ptr<rectangle> r1, shared_ptr<rectangle> r2);
void ClearScreen();

/**
Multiplies two btMatrix3x3 matricies together. m*v

@param m First matrix
@param v Second Matrix
@return The resulting matrix
*/
SIMD_FORCE_INLINE btVector3 
MeshEnv::multOperator(const btMatrix3x3& m, const btVector3& v)
{
	return btVector3(m[0].dot(v), m[1].dot(v), m[2].dot(v));
}

/**
Custom collision filter callback to ensure spine cylinders don't collide with each other
*/
struct AGCustomFilterCallback : public btOverlapFilterCallback
{
	//const after function name means that the function will not change any of the member variables of the
	// class/struct. In this case, there are no member variables to change, so this is true.
	//return true when pairs need collision

	/**
	Adding additional checks to see if collision detections will be added to rest of collision handling for
	spine cylinders.

	@param proxy0 First object
	@param proxy1 Second object
	@return True if need collision, false otherwise
	*/
	virtual bool needBroadphaseCollision(btBroadphaseProxy* proxy0, btBroadphaseProxy* proxy1) const override
	{
		//lines to still implement the collision mask filtering completed before
		// obj0 collides with obj1 and vice versa
		bool collides = (proxy0->m_collisionFilterGroup & proxy1->m_collisionFilterMask) != 0;
		collides = collides && (proxy0->m_collisionFilterMask & proxy1->m_collisionFilterGroup);

		//additional logic to handle spine-spine collisions

		//if they are both spine cylinders
		if ((proxy0->m_collisionFilterGroup == proxy1->m_collisionFilterGroup) 
			&& proxy0->m_collisionFilterGroup == COL_SPINE)
		{
			//check their spine index and their virtibrae index
			// HACK: CNT program only uses collisionObjects in proxies. Will not work if others used.
			btCollisionObject* obj0 = static_cast<btCollisionObject*>(proxy0->m_clientObject);
			btCollisionObject* obj1 = static_cast<btCollisionObject*>(proxy1->m_clientObject);

			// //Check if same tube
			// if ((obj0->getSpine() == obj1->getSpine()) && obj0->getSpine()!=0)
			// //if ((obj0->getSpine() == obj1->getSpine()))
			// {
			// 	int vertebraeDiff = abs(obj0->getVertebrae() - obj1->getVertebrae());
			// 	//check if adjacent vertebrae
			// 	if (vertebraeDiff == 1)
			// 	{
			// 		return false;
			// 	}
				
			// }
		}

		return collides;
	}
};

/**
Rotates the given rigid body object about the origin as the rotation matrix specifies

@param rotMatrix Matrix that specifies rotation transformation
@param obj Rigid body to be rotated
@return void
*/
void MeshEnv::rotateAboutOrigin(btMatrix3x3 &rotMatrix, btRigidBody* obj)
{
	btMatrix3x3 basis = obj->getWorldTransform().getBasis();
	btVector3 origin = obj->getWorldTransform().getOrigin();
	basis *= rotMatrix;
	origin = multOperator(rotMatrix, origin);
	obj->getWorldTransform().setBasis(basis);
	obj->getWorldTransform().setOrigin(origin);

}

/**
Rotates the given rigid body object about the origin as the rotation matrix specifies

@param rotMatrix Matrix that specifies rotation transformation
@param obj Rigid body to be rotated
@param shift The amount the object should be translated
@return void
*/
void MeshEnv::rotateAndShift(btMatrix3x3 &rotMatrix, btRigidBody* obj, btVector3 &shift)
{
	rotateAboutOrigin(rotMatrix, obj);
	btTransform* objTransform = &obj->getWorldTransform();
	objTransform->setOrigin(objTransform->getOrigin() + shift);
}

/**
Shifts a rigid body

@param obj The object to be shifted
@param shift The amount to be shifted
@return void
*/
void MeshEnv::shift(btRigidBody* obj, btVector3& shift)
{
	btTransform* objTransform = &obj->getWorldTransform();
	objTransform->setOrigin(objTransform->getOrigin() + shift);
}

/**
Converts numbers with some units to angstroms

@param unit The current unit
@param val The current value
@return the value in angstroms
*/
double MeshEnv::convertUnits(string unit, double val)
{
	if (unit.compare("mm") == 0 || unit.compare("millimeter") == 0)
	{
		return val*1000000;
	}
	else if (unit.compare("um") == 0 || unit.compare("micrometer") == 0)
	{
		return val*1000;
	}
	else if (unit.compare("nm") == 0 || unit.compare("nanometer") == 0)
	{
		return val * 1;
	}
	else if (unit.compare("pm") == 0 || unit.compare("picometer") == 0)
	{
		return val * .001;
	}
	else if (unit.compare("A") == 0 || unit.compare("angstrom") == 0)
	{
		return val * .01;
	}
	else
	{
		return std::numeric_limits<int>::min();
	}
}

//Initializes the physics world with all of the carbon nanotubes as specified by the input XML.
void	MeshEnv::initPhysics()
{
	MeshEnv::initPhysics(btScalar(100.));
}


//Initializes the physics world with all of the carbon nanotubes as specified by the input XML.
//@param camDistance The initial camera distance
void	MeshEnv::initPhysics(float camDistance)
{

	//////////////////// ENVIRONMENT PARAMETERS /////////////////////////////////////////////
	// setTexturing(true); // replace this with enableTexture(bool enable) in GL_ShapeDrawer.h
	// setShadows(true);
	// setCameraDistance(camDistance);

	//Initialize random number generation
	initRandomNumGen();

	//// XML Parameter Instantiation ////
	string outputFolderPath = "";
	numTubes = 1;
	m_xdim = 150.0;
	m_ydim = 100.0;
	m_zdim = 150.0;
	vector<int> chirality(0);
	int chirCount = 0;
	m_tubeList = make_shared<vector<shared_ptr<tube>>>(vector<shared_ptr<tube>>(0));
	m_waveList = vector<wave>(0);
	bool done = false;
	
	while (!done){
		try{
			shared_ptr<xml_document<>> doc = make_shared<xml_document<>>(); //create xml object
			shared_ptr<file<>> xmlFile = make_shared<file<>>(file<>(inputXMLPath.c_str())); //open file
			doc->parse<0>(xmlFile->data()); //parse contents of file
			xml_node<>* currNode = doc->first_node(); //gets the node "Document" or the root node
			
			//OUTPUT FOLDER
			currNode = currNode->first_node(); //Output folder validated below
			outputFolderPath = currNode->value();
			outputFolderPath = fixPath(outputFolderPath);
			if (outputFolderPath.back() != '/')
			{
				outputFolderPath += "/";
			}

			//NUMTUBES NODE
			currNode = currNode->next_sibling(); 
			numTubes = atoi(currNode->value());
			//Input validation
			if (numTubes <= 0)
			{
				printf("Configuration Error: Must enter positive integer number of tubes.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			m_tubeList = make_shared<vector<shared_ptr<tube>>>(vector<shared_ptr<tube>>(0));
			//FRICTION NODE
			currNode = currNode->next_sibling(); 
			friction = atof(currNode->value());
			//incorrect range
			if (friction <= 0)
			{
				printf("Configuration Error: Must enter positive friction value.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}

			//GRAVITY NODE
			currNode = currNode->next_sibling(); 
			gravity = atof(currNode->value());
			//incorrect range
			if (gravity >= 0)
			{
				printf("Configuration Error: For object to fall down, gravity must be negative\n");
				system("pause");
				exit(EXIT_FAILURE);
			}

			// SPACING NODE //
			currNode = currNode->next_sibling();
			minSpacing = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value())) / 2.0;
			//incorrect units
			if (minSpacing == std::numeric_limits<int>::min() / 2.0)
			{
				printf("Configuration Error: Incorrect units for spacing.\nRefer to manual"
					" for valid unit entries.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			//incorrect range
			else if (minSpacing <= 0)
			{
				printf("Configuration Error: Must enter positive spacing value.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END SPACING NODE //

			// BUNDLE LENGTHS NODE //
			currNode = currNode->next_sibling()->next_sibling();
			lmin = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value()));
			lmax = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->value()));
			//incorrect units
			if (lmin == std::numeric_limits<int>::min() || lmax == std::numeric_limits<int>::min())
			{
				printf("Configuration Error: Incorrect units for lengths.\nRefer to manual"
					" for valid unit entries.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			//incorrect pairing
			else if (lmin > lmax)
			{
				printf("Configuration Error: Lmin must be less than or equal to Lmax\n");
				system("pause");
				exit(EXIT_FAILURE);
			} 
			//incorrect range
			else if(lmin <= 0 || lmax <= 0)
			{
				printf("Configuration Error: Must enter positive length values.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END BUNDLE LENGTHS NODE //

			// BUNDLE LOOP CIRCUMFERENCE //

			currNode = currNode->next_sibling();
			double loopCircumference = convertUnits(string(currNode->first_node()->value()), 
				atof(currNode->first_node()->next_sibling()->value()));

			if (loopCircumference <= 0)
			{
				printf("Configuration Error: Circumference must be larger than 0.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			curvMax = 2 * SIMD_PI / (loopCircumference);

			// END BUNDLE LOOP CIRCUMFERENCE //

			// DEVICE DIMENSIONS NODE //
			currNode = currNode->next_sibling();
			m_xdim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value())) / 2.0;
			m_ydim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->value())) / 2.0;
			m_zdim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->next_sibling()->value())) / 2.0;
			//incorrect units
			if (m_xdim == std::numeric_limits<int>::min() / 2.0 || m_ydim == std::numeric_limits<int>::min() / 2.0 || m_zdim == std::numeric_limits<int>::min() / 2.0)
			{
				printf("Configuration Error: Incorrect units for device dimensions.\nRefer to manual"
					" for valid unit entries.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			//incorrect range
			else if (m_xdim <= 0 || m_ydim <= 0 || m_zdim <= 0)
			{
				printf("Configuration Error: Must enter positive device dimensions.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			else if (m_xdim < 1.5*lmax || m_zdim < 1.5*lmax)
			{
				printf("Device dimensions are too small for the given bundle lengths. Please increase device size.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END DEVICE DIMENSIONS NODE //

			// BUNDLE PLACEMENT NODE //
			currNode = currNode->next_sibling();
			x_placement_extent = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value())) / 2.0;
			z_placement_extent = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->value())) / 2.0;
			placement_width = vector<double>(0.);
			placement_width.push_back(x_placement_extent * 2.0);
			placement_width.push_back(z_placement_extent * 2.0);
			//incorrect units
			if (x_placement_extent == std::numeric_limits<int>::min() / 2.0 || z_placement_extent == std::numeric_limits<int>::min() / 2.0)
			{
				printf("Configuration Error: Incorrect units for bundle placement area.\nRefer to manual"
					" for valid unit entries.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			//incorrect range
			else if (x_placement_extent <= 0 || z_placement_extent <= 0)
			{
				printf("Configuration Error: Must enter positive values for bundle placement area.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			else if (x_placement_extent < lmax || z_placement_extent < lmax)
			{
				printf("Lengths in bundle placement area are too small for the given bundle lengths. Please increase device size.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			else if (x_placement_extent >= m_xdim || z_placement_extent >= m_zdim)
			{
				printf("Dimensions in bundle placement area must be smaller than the lengths in the device dimensions.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END BUNDLE PLACEMENT NODE //

			// CHIRALITY NODE //
			currNode = currNode->next_sibling();
			for (xml_node<>* child = currNode->first_node(); child; child = child->next_sibling())
			{
				chirCount++;
			}
			currNode = currNode->first_node(); //Go to first cnt
			for (int i = 0; i < chirCount; i++)
			{
				chirality.insert(chirality.begin() + 2 * i, atoi(currNode->first_node()->value()));
				chirality.insert(chirality.begin() + 2 * i + 1, atoi(currNode->first_node()->next_sibling()->value()));
				
				int temp_n = atoi(currNode->first_node()->value());
				int temp_m = atoi(currNode->first_node()->next_sibling()->value());

				//invalid m
				if (temp_n < temp_m)
				{
					printf("Configuration Error: For CNTs, n must be larger than m.\n");
					system("pause");
					exit(EXIT_FAILURE);
				} 
				//invalid m or n range
				else if (temp_n < 0 || temp_m < 0)
				{
					printf("Configuration Error: Please provide hamada parameters that are greater than or equal to 0.\n");
					system("pause");
					exit(EXIT_FAILURE);
				}

				currNode = currNode->next_sibling();
			}
			delete currNode;
			done = true;
		}
		catch (runtime_error err)
		{
			cout << err.what();
			cout << "\n";
			cout << "Continue? [y/n]: ";
			std::string temp;
			cin >> temp;
			if (temp.compare("y") != 0)
			{
				system("pause");
				exit(EXIT_FAILURE);
			}
			char *inputXMLPathArray = new char[xmlArrayLength];
			system("cls");
			cout << "Enter config xml path (Example in program files directory):\n";
			cin.ignore();
			cin.getline(inputXMLPathArray, xmlArrayLength);
			inputXMLPath = inputXMLPathArray;
			delete[] inputXMLPathArray;
		}
	}


	runID = "d";
	string response;
	{
		time_t timer;
		struct tm currTime;
		if (time(&timer) != -1)
		{
			struct tm* err = std::localtime(&timer);
			// if (err)
			// {
			// 	printf("Invalid argument to localtime.\n");
			// 	system("pause");
			// 	exit(EXIT_FAILURE);
			// }
		}

		// #pragma warning(suppress: 6001)
		// runID = runID + to_string(currTime.tm_mday) + "." + to_string(currTime.tm_mon + 1) + "."
		// 	+ to_string(currTime.tm_year % 100) + "_t" + to_string(currTime.tm_hour) + "." +
		// 	to_string(currTime.tm_min) + "." + to_string(currTime.tm_sec) + "_c" + to_string(numTubes) + 
		// 	"_x" + to_string(static_cast<int>(m_xdim / 5.0)) + "y" + to_string(static_cast<int>(m_ydim / 5.0)) + "z" +
		// 	to_string(static_cast<int>(m_zdim / 5.0));
	}

	outputPath = outputFolderPath + runID + "/";
	std::wstring wide_string(outputPath.begin(), outputPath.end());
	if (std::experimental::filesystem::create_directories(wide_string.c_str()) == 0)
	{
		printf("Invalid output folder path.\n");
		system("pause");
		exit(EXIT_FAILURE);
	}

	//Copy the input xml config file to the destination folder.
	wstring w_inputXMLPathString;
	w_inputXMLPathString.assign(inputXMLPath.begin(), inputXMLPath.end());
	const wchar_t* w_inputXMLPath = w_inputXMLPathString.c_str();

	string xmlCopy = outputPath + runID + ".xml";
	wstring w_runIDstring;
	w_runIDstring.assign(xmlCopy.begin(), xmlCopy.end());
	const wchar_t* w_runID = w_runIDstring.c_str(); //allocated on stack. Don't need to delete
	bool success = std::experimental::filesystem::copy_file(w_inputXMLPath, w_runID, std::experimental::filesystem::copy_options::overwrite_existing);
	//Check to see if file copy was a success
	if (!success)
	{
		std::cout << "Configuration file could not be copied into results directory.\n";
		exit(EXIT_FAILURE);
	}

	/////////////////////////////// END OF ENVIRONMENT PARAMETERS //////////////////////////////////

	///collision configuration contains default setup for memory, collision setup
	m_collisionConfiguration = new btDefaultCollisionConfiguration();

	///use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
	m_dispatcher = new	btCollisionDispatcher(m_collisionConfiguration);

	m_broadphase = new btDbvtBroadphase();

	///the default constraint solver. For parallel processing you can use a different solver (see Extras/BulletMultiThreaded)
	btSequentialImpulseConstraintSolver* sol = new btSequentialImpulseConstraintSolver;
	m_solver = sol;

	m_dynamicsWorld = new btDiscreteDynamicsWorld(m_dispatcher, m_broadphase, m_solver, m_collisionConfiguration);
	// m_dynamicsWorld->setDebugDrawer(&gDebugDraw);
	// gDebugDraw.setDebugMode(1024);

	//setting up the new filter callback for non-adjacent spine cylinder collisions
	btOverlapFilterCallback* filterCallback = new AGCustomFilterCallback();
	m_dynamicsWorld->getPairCache()->setOverlapFilterCallback(filterCallback);

	m_dynamicsWorld->setGravity(btVector3(0, gravity, 0));

	///create a few basic rigid bodies
	//This is the ground plane box, made the same way as in Blender: sizes are from center of object to outside
	//btScalar gndDim = btScalar(25.);
	btScalar planeHlfExtThick = btScalar(.5); //half extent thickness of plane
	btBoxShape* groundShape = new btBoxShape(btVector3(m_xdim, planeHlfExtThick , m_zdim));

	groundShape->setMargin(.0);
	//adding the ground shape to the aligned object array for easy cleanup later
	m_collisionShapes.push_back(groundShape);
	btTransform planeTransform;
	planeTransform.setIdentity();
	planeTransform.setOrigin(btVector3(0, -planeHlfExtThick, 0));
	btScalar planeMass = 0.;

	//creates the necessary parts of an object and adds the object to the dynamics world
	//  #not sure what to do with handling the memory of the rigid body, the dynamics world
	//  # might take care of it becuase it has all of the pointers in the world

	int planeCollidesWith = COL_EVERYTHING & ~COL_PLANE;

	//I created this overloaded method to allow for the collision masking
	btRigidBody* body = localCreateRigidBody(planeMass, planeTransform, groundShape, COL_PLANE, planeCollidesWith);
	body->setFriction(friction);
	
	//
	//Create Walls of Simulation
	//
	planeTransform.setOrigin(btVector3(0, 0, 0));

	//x wall plane
	btBoxShape* xWallShape = new btBoxShape(btVector3(planeHlfExtThick, m_ydim, m_zdim));
	xWallShape->setMargin(0.);
	m_collisionShapes.push_back(xWallShape);

	//+x plane
	body = localCreateRigidBody(planeMass, planeTransform, xWallShape, COL_PLANE, planeCollidesWith);
	//btMatrix3x3 planeRot = btMatrix3x3(btQuaternion(btVector3(0, 0, 1), SIMD_HALF_PI));
	btVector3 planeShift(m_xdim + planeHlfExtThick, m_ydim, 0);
	shift(body, planeShift);
	//body->setDrawable(FALSE); // feature added by AG

	//-x plane
	body = localCreateRigidBody(planeMass, planeTransform, xWallShape, COL_PLANE, planeCollidesWith);
	planeShift = btVector3(-(m_xdim + planeHlfExtThick), m_ydim, 0);
	shift(body, planeShift);
	// body->setDrawable(false); // feature added by AG

	//z wall plane
	btBoxShape* zWallShape = new btBoxShape(btVector3(m_xdim, m_ydim, planeHlfExtThick));
	zWallShape->setMargin(0.);
	m_collisionShapes.push_back(zWallShape);

	//+z plane
	body = localCreateRigidBody(planeMass, planeTransform, zWallShape, COL_PLANE, planeCollidesWith);
	//planeRot = btMatrix3x3(btQuaternion(btVector3(1, 0, 0), SIMD_HALF_PI));
	planeShift = btVector3(0, m_ydim, m_zdim + planeHlfExtThick);
	shift(body, planeShift);
	// body->setDrawable(FALSE); // feature added by AG

	//-z plane
	body = localCreateRigidBody(planeMass, planeTransform, zWallShape, COL_PLANE, planeCollidesWith);
	planeShift = btVector3(0, m_ydim, -(m_zdim + planeHlfExtThick));
	shift(body, planeShift);
	// body->setDrawable(false); // feature added by AG
	
	//////////////////////////////////////////////////////////////////////////////////////////////

	//Set up first wave
	nextWave = true;
	currWave = 0;
	placementArea = placement_width[0] * placement_width[1];
	auto firstWave = wave();
	firstWave.startIdx = 1;
	m_waveList.push_back(firstWave);
	
	//MeshEnv::keyboardCallback('D', 0, 0); //starts without rendering anything, release code
	// MeshEnv::keyboardCallback('d', 0, 0); //debug code, no ending
}



void MeshEnv::addWave()
{
	//////////////////////// Nanotube Creation //////////////////////////////////////////////////

	//Collision mask
	int spineCollidesWith = COL_SPINE | COL_PLANE | COL_BOX;

	//boolean that specifies whether or not tubes in this wave have the sum of their Aabb's 
	// below the placement area limit
	bool belowArea = true;
	double currWaveArea = 0;

	//
	// NANOTUBE CREATION LOOP
	//
	for (auto i = m_waveList[currWave].startIdx; i < numTubes && belowArea; i++)
	{

		//
		// Physical Nanotube Parameters
		//
		double lrange = lmax - lmin;
		double radius = 20; //radius of the tube
		double length = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX))*lrange + lmin; //A
		if (radius > m_waveList[currWave].maxRadius) { m_waveList[currWave].maxRadius = radius; }
		shared_ptr<tube> curr_tube;
		try{
			//
			// create tube
			//
			curr_tube = make_shared<tube>(tube(i, radius*2.0, curvMax, length, minSpacing));
			m_tubeList->push_back(curr_tube);
		}
		catch (runtime_error err)
		{
			cout << err.what();
			cout << "\n";
			system("pause");
			exitPhysics();
			exit(EXIT_FAILURE);

		}
		double height = curr_tube->getCylHeight();
		double tubeSeparation = curr_tube->getTubeSpacing();
		int numSection = curr_tube->getNumSections();

		m_waveList[currWave].cylCount += numSection; //keeps total count of cylinders


		//NOTE: USE PLANE COORDINATES (AS VIEWED IN WIREFRAME MODE) FOR ALL REFERENCES

		//
		// Pre-Build Parameter Calculations
		//
		//location of the first cylinder, and the prev cylinder in the for loop below
		double prevYPos = height / 2;
		/*Each tube will refered to as a spine made up of individual vertebrae. Each tube will have a spine
		number. Each cylinder in that spine will have a member field specifying the spine number as well.
		Each cylinder also represents a vertebrae in the spine. Each vertebrae will have an index as well.
		These numbers are necessary for broadphase collision filtering.
		*/
		double currYPos = prevYPos; // center of the next cylinder in the chain. To be calculated.
		double constraintShift = (height + tubeSeparation) / 2.0;
		//Determines actual tube length after convergence methods
		double tubeLength = (numSection - 1)*(tubeSeparation + height) + height;

		//double allowableXPos = m_xdim - sqrt(pow((tubeLength / 2.0), 2) + radius*radius);
		//double allowableZPos = m_xdim - sqrt(pow((tubeLength / 2.0), 2) + radius*radius);

		//Non-random positions
		//btScalar xpos = 0;
		//btScalar zpos = 0;
		//btScalar theta = 0;
		//
		//Random number generation
		//
		//btScalar xpos = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX))*allowableXPos*2. - allowableXPos;
		//btScalar zpos = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX))*allowableZPos*2. - allowableZPos;
		btScalar theta = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX))*SIMD_PI - SIMD_PI / 2.;


		btVector3 shiftToOrigin(0, tubeLength / 2., 0);
		btMatrix3x3 rotMatrix = btMatrix3x3(btQuaternion(btVector3(0, 1, 0), theta));
		rotMatrix*=(btMatrix3x3(btQuaternion(btVector3(0, 0, 1), btScalar(-SIMD_PI / 2.))));
		shiftToOrigin = -multOperator(rotMatrix, shiftToOrigin);
		//btVector3 shift = btVector3(xpos, radius + minSpacing + 2 * (radius + minSpacing) * (i - 1), zpos) + shiftToOrigin; //ground test tubes
		btVector3 shift = shiftToOrigin;// + btVector3(xpos, 5.0 + radius + minSpacing + 2 * (radius + minSpacing) * (currWave - 1), zpos);


		/// Create Dynamic Objects
		btTransform startTransform;
		//Set the basis to the identity matrix which is a basic Euler matrix. See ms to need to be called before
		// trying to do any sort of transforms.
		startTransform.setIdentity();
		//after basis is set, setting origin makes more sense


		//Create the first cylinder manually and then the rest programmatically
		//create a few dynamic rigidbodies
		// Re-using the same collision is better for memory usage and performance
		// NOTE: To get the correct aabb for the cylinders, you must define radius for both x 
		// and z components (as this is a cylinder aligned around the y axis)
		btCylinderShape* spineShape = new btCylinderShape(btVector3(radius, prevYPos, radius));
		//Only need to push back each collision shape not each object with a collision shape
		m_collisionShapes.push_back(spineShape);
		//startTransform will be used for all subsequent created objects
		startTransform.setOrigin(btVector3(0, prevYPos, 0));
		//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
		btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
		btScalar	spineMass(cylMass); //mass of the object
		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (spineMass != 0.f);
		btVector3 localInertiaSpine(0, 0, 0);
		if (isDynamic)
			spineShape->calculateLocalInertia(spineMass, localInertiaSpine);
		//Must not reuse this as the Motion state must change for each cylinder
		shared_ptr<btRigidBody::btRigidBodyConstructionInfo> rbInfoSpine =
			make_shared<btRigidBody::btRigidBodyConstructionInfo>(btRigidBody::btRigidBodyConstructionInfo
			(spineMass, myMotionState, spineShape, localInertiaSpine));


		btRigidBody* prevSpineCyl = new btRigidBody(*rbInfoSpine);
		//rotate the object accordingly
		rotateAndShift(rotMatrix, prevSpineCyl, shift);
		//Add cylinder to world
		m_dynamicsWorld->addRigidBody(prevSpineCyl, COL_SPINE, spineCollidesWith); //Spine cylinder1
		//set spine
		// prevSpineCyl->setSpine(i);
		// prevSpineCyl->setVertebrae(1);
		curr_tube->addCyl(prevSpineCyl);

		//set the constraint angle
		btScalar constraintAngle = curr_tube->getCurvature()*(curr_tube->getCylHeight()+tubeSeparation); //[rads]

		//loop to create the rest of the chain that is the nanotube
		for (int j = 2; j <= numSection; j++)
		{
			//
			//Common parameters
			//

			currYPos += (height + tubeSeparation); //get its height
			startTransform.setOrigin(btVector3(0, currYPos, 0)); //set its height


			//
			//ith spine cylinder
			//
			myMotionState = new btDefaultMotionState(startTransform);
			rbInfoSpine = make_shared<btRigidBody::btRigidBodyConstructionInfo>
				(btRigidBody::btRigidBodyConstructionInfo(spineMass, myMotionState, spineShape, localInertiaSpine));
			btRigidBody* currSpineCyl = new btRigidBody(*rbInfoSpine);
			//Rotations on rigid body
			rotateAndShift(rotMatrix, currSpineCyl, shift);
			//Add spine cylinder to world
			m_dynamicsWorld->addRigidBody(currSpineCyl, COL_SPINE, spineCollidesWith);
			// currSpineCyl->setSpine(i);
			// currSpineCyl->setVertebrae(j);
			curr_tube->addCyl(currSpineCyl);


			//
			//Cone-twist constraint
			//

			//add (i-1)th cone twist constraint
			btTransform transformA, transformB;
			transformA.setIdentity();
			transformB.setIdentity();
			transformA.setOrigin(btVector3(0.0, constraintShift, 0.0));
			transformB.setOrigin(btVector3(0.0, -constraintShift, 0.0));
			//set cone twist constraint to be on tube axis and halfway between the two
			btConeTwistConstraint* coneTwist = new btConeTwistConstraint(*prevSpineCyl, *currSpineCyl, transformA, transformB);
			coneTwist->setLimit(constraintAngle, constraintAngle, btScalar(0.0), btScalar(1.0), btScalar(1.0), btScalar(1.0));
			//set constraint to max value so there is no breaking
			coneTwist->setBreakingImpulseThreshold(btScalar(BT_LARGE_FLOAT));
			coneTwist->setDbgDrawSize(btScalar(0.f));
			m_dynamicsWorld->addConstraint(coneTwist);


			//
			//Update spine cylinder pointer
			//

			//set current to previous for next loop
			prevSpineCyl = currSpineCyl;

		}
		curr_tube->setFinished();
		tubesCreated++;
		currWaveArea += curr_tube->getRectangle()->getArea();
		//if area is large enough, time to tile
		if (currWaveArea > .7*placementArea || i == numTubes - 1)
		{
			belowArea = false;
			m_waveList[currWave].endIdx = i;
			if (i < numTubes - 1)
			{
				auto nextWave = wave();
				nextWave.startIdx = i + 1;
				m_waveList.push_back(nextWave);
			}
		}

	}
	//////////////////////////////////////////////////////////////////////////////////////////
	tileWave();
	//This wave has finished. 

}

//Takes the Aabb  of each bundle and makes a tiles them on the xz plane to get a good spread
// of nanotubes for the wave
void MeshEnv::tileWave()
{
	double levelWidth = placement_width[0x01 & rotation];
	auto rectVec = vector<shared_ptr<rectangle>>(m_waveList[currWave].endIdx - m_waveList[currWave].startIdx + 1);
	//build rectVec
	auto range = m_waveList[currWave].endIdx - m_waveList[currWave].startIdx;
	for (auto i = 0; i <= range; i++)
	{
		rectVec[i] = (*m_tubeList)[m_waveList[currWave].startIdx+i-1]->getRectangle();
	}
	sort(rectVec.begin(), rectVec.end(), compareRect);
	auto rectList = list<shared_ptr<rectangle>>(0);
	//move vector to list
	for (auto i = 0; i < rectVec.size(); i++)
	{
		rectList.push_back(rectVec[i]);
	}

	//First-Fit Decreasing Height (FFDH)
	auto done = false;
	auto levelVector = vector<level>(0);
	auto currLevel = 0;
	auto currHeight = 0.;
	auto nextHeight = 0.;
	while (!done)
	{
		auto currLvlWidth = 0.;
		levelVector.push_back(level());
		levelVector.back().levelNum = currLevel;
		levelVector[currLevel].lvlHeight = rectList.front()->getHeight();
		nextHeight += levelVector[currLevel].lvlHeight;
		auto erasePrev = false;
		for (auto itr = rectList.begin(); itr != rectList.end(); ++itr)
		{
			if (erasePrev)
			{
				auto prevItr = prev(itr);
				rectList.erase(prevItr);
			}
			currLvlWidth += (*itr)->getWidth();
			if (currLvlWidth < levelWidth)
			{
				erasePrev = true;
				auto rectElem = (*itr);
				auto currPos = position();
				currPos.x = currLvlWidth - rectElem->getWidth();
				currPos.z = currHeight;
				rectElem->setPosition(currPos);
				levelVector[currLevel].rectsOnLevel.push_back(rectElem);
			}
			else
			{
				currLvlWidth -= (*itr)->getWidth();
				erasePrev = false;
			}
		}
		levelVector[currLevel].usedWidth = currLvlWidth;
		//Takes care of last element in the case that the last element was added to level
		if (erasePrev)
		{
			rectList.erase(prev(rectList.end()));
		}

		if (rectList.size() == 0)
		{
			done = true;
		}

		currLevel++;
		currHeight = nextHeight;
	}
	//Move through all of the rectangles and move their respective cylinders to the appropriate location
	double angle = SIMD_HALF_PI*rotation;
	btMatrix3x3 rotMatrix = btMatrix3x3(btQuaternion(btVector3(1, 0, 0), angle));
	btVector3 cylShift;
	btVector3 zero_vec = btVector3(0, 0, 0);
	double x_sign, z_sign, x_shift, z_shift;
	auto rotState = rotation % 4;
	switch (rotState)
	{
		case 0:
			x_sign = z_sign = -1;
			break;
		case 1:
			x_sign = 1;
			z_sign = -1;
			break;
		case 2:
			x_sign = z_sign = 1;
			break;
		case 3:
			x_sign = -1;
			z_sign = 1;
			break;
		default:
			x_sign = z_sign = 0;
	}
	for (auto itr = rectVec.begin(); itr != rectVec.end(); ++itr)
	{
		double width_hlf = 0;
		double height_hlf = 0;
		auto curr_rect = (*itr);
		auto rect_pos = curr_rect->getPosition();
		switch (rotState)
		{
			case 0:
				x_shift = rect_pos.x;
				z_shift = rect_pos.z;
				width_hlf = curr_rect->getWidth() / 2.0;
				height_hlf = curr_rect->getHeight() / 2.0;
				break;
			case 1:
				x_shift = -1 * rect_pos.z;
				z_shift = rect_pos.x;
				height_hlf = curr_rect->getWidth() / 2.0;
				width_hlf = curr_rect->getHeight() / 2.0;
				break;
			case 2:
				x_shift = -1 * rect_pos.x;
				z_shift = -1 * rect_pos.z;
				width_hlf = curr_rect->getWidth() / 2.0;
				height_hlf = curr_rect->getHeight() / 2.0;
				break;
			case 3:
				x_shift = rect_pos.z;
				z_shift = -1 * rect_pos.x;
				height_hlf = curr_rect->getWidth() / 2.0;
				width_hlf = curr_rect->getHeight() / 2.0;
				break;
			default:
				x_shift = z_shift = 0;
		}
		auto currCylList = (*m_tubeList)[curr_rect->getTubeNum()-1]->getCylList();
		cylShift = btVector3(x_sign*(x_placement_extent - width_hlf) + x_shift,
			meshYPeak + m_waveList[currWave].maxRadius + 15.0, z_sign*(z_placement_extent - height_hlf) + z_shift);
		for (auto itrCyl = currCylList->begin(); itrCyl != currCylList->end(); ++itrCyl)
		{
			auto cyl = (*itrCyl);
			btVector3 point = cyl->getCenterOfMassPosition();
			btTransform transform = cyl->getCenterOfMassTransform();
			transform.setOrigin(zero_vec);
			cyl->setCenterOfMassTransform(transform);
			rotateAndShift(rotMatrix, cyl, zero_vec);
			point = point.rotate(btVector3(0., 1., 0.), -angle);
			transform = cyl->getCenterOfMassTransform();
			transform.setOrigin(point);
			cyl->setCenterOfMassTransform(transform);
			shift(cyl, cylShift);
		}
	}

	rotation++; //next wave gets rotated differently
}

/**
steps the simulation, then it renders the dynamics world. Also responsible for checking if the 
simulation is done and saving the data.
*/
void MeshEnv::clientMoveAndDisplay()
{
	//simple dynamics world doesn't handle fixed-time-stepping
	btClock m_clock;
	btScalar ms = (btScalar)m_clock.getTimeMicroseconds();

	///step the simulation
	if (m_dynamicsWorld)
	{
		if (nextWave){ addWave(); nextWave = false; }
		m_dynamicsWorld->stepSimulation(ms / 1000000.f); ///////////////STEP SIMULATION////////////
		if ((checkCntr&bitMask) == 0 || manualEndSimulation){
			/*each 8th step, I would like to check if simulation is done*/
			bool waveComplete = true;
			//keeps the total number of resting cylinders
			auto sleepTotCntr = 0;
			//iterate over the list of tubes while we think that the simulation might be done
			for (auto i = m_waveList[currWave].startIdx-1; i <= m_waveList[currWave].endIdx-1 && !manualEndSimulation; i++)
			{
				auto tempCylList = (*m_tubeList)[i]->getCylList();
				for (auto itrCyl = tempCylList->begin(); itrCyl < tempCylList->end(); ++itrCyl)
				{
					if ((*itrCyl)->getActivationState() == ISLAND_SLEEPING ||
						(*itrCyl)->getActivationState() == WANTS_DEACTIVATION)
					{
						sleepTotCntr++;
					}
				}
			}
			//if 10% of the cylinders are ready to or are sleeping, then the tube is done.
			if (float(sleepTotCntr) / float(m_waveList[currWave].cylCount) < .1)
			{
				waveComplete = false;
			}

			if ((m_waveList[currWave].endIdx == numTubes - 1 && waveComplete) || manualEndSimulation)
			{
				//outputResults();
				exitPhysics();
				exit(EXIT_SUCCESS);
			}

			if (waveComplete)
			{
				nextWave = true;
				makeTubesStatic(m_waveList[currWave].startIdx, m_waveList[currWave].endIdx);
				currWave++;
			}

		}
		checkCntr++;
		//optional but useful: debug drawing
		if (toDebugDraw){
			m_dynamicsWorld->debugDrawWorld();
		}
	}
	
	// AHD: commented not sure what is does these functions do.
	// //if (toDebugDraw)
	// //{
	// 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// 	//Rendering function. Uses
	// 	renderme();
	// 	//APIENTRY functions -> cannot edit
	// 	glFlush();
	// 	swapBuffers();
	// //}

}

//Takes tubes from tube list and makes them static as they have already contributed to the mesh
void MeshEnv::makeTubesStatic(int start_tube, int end_tube)
{
	btVector3 max, min;
	for (auto i = m_waveList[currWave].startIdx-1; i <= m_waveList[currWave].endIdx-1; i++)
	{
		auto tempCylList = (*m_tubeList)[i]->getCylList();
		for (auto itrCyl = tempCylList->begin(); itrCyl < tempCylList->end(); ++itrCyl)
		{
			(*itrCyl)->getAabb(min, max);
			if (max.getY() > meshYPeak){ meshYPeak = max.getY(); }
			(*itrCyl)->setActivationState(DISABLE_SIMULATION);
		}
	}
}

// AHD: commented not sure what is does these functions do.
// /**
// renders the dynamics world as it stands currently
// */
// void MeshEnv::displayCallback(void) {


// 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

// 	renderme();

// 	//optional but useful: debug drawing to detect problems
// 	if (m_dynamicsWorld)
// 		m_dynamicsWorld->debugDrawWorld();

// 	glFlush();
// 	swapBuffers();
// }

/**
Exits the current simulation and starts over
*/
void	MeshEnv::clientResetScene()
{
	initPhysics(exitPhysics());
}

/**
Exits the current simulation, cleaning up all data
*/
float	MeshEnv::exitPhysics()
{

	//cleanup in the reverse order of creation/initialization

	//remove the rigidbodies from the dynamics world and delete them
	int i;
	for (i = m_dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; i--)
	{
		btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
		{
			delete body->getMotionState();
		}
		m_dynamicsWorld->removeCollisionObject(obj);
		delete obj;
	}

	//delete collision shapes
	for (int j = 0; j<m_collisionShapes.size(); j++)
	{
		btCollisionShape* shape = m_collisionShapes[j];
		delete shape;
	}
	m_collisionShapes.clear();

	delete m_dynamicsWorld;

	delete m_solver;

	delete m_broadphase;

	delete m_dispatcher;

	delete m_collisionConfiguration;

	return 100.0;
}

///The MyOverlapCallback is used to show how to collect object that overlap with a given bounding box defined by aabbMin and aabbMax. 
///See m_dynamicsWorld->getBroadphase()->aabbTest.
struct	MyOverlapCallback : public btBroadphaseAabbCallback
{
	btVector3 m_queryAabbMin;
	btVector3 m_queryAabbMax;

	int m_numOverlap;
	MyOverlapCallback(const btVector3& aabbMin, const btVector3& aabbMax) : m_queryAabbMin(aabbMin), m_queryAabbMax(aabbMax), m_numOverlap(0)	{}
	virtual bool	process(const btBroadphaseProxy* proxy) override
	{
		btVector3 proxyAabbMin, proxyAabbMax;
		btCollisionObject* colObj0 = static_cast<btCollisionObject*>(proxy->m_clientObject);
		colObj0->getCollisionShape()->getAabb(colObj0->getWorldTransform(), proxyAabbMin, proxyAabbMax);
		if (TestAabbAgainstAabb2(proxyAabbMin, proxyAabbMax, m_queryAabbMin, m_queryAabbMax))
		{
			m_numOverlap++;
		}
		return true;
	}
};

// /**
// Takes a screenshot of the mesh
// */
// int MeshEnv::takeScreenshot()
// {
// 	int code = 0; //return code

// 	//window handle needs wide char array
// 	wstring w_runIDstring;
// 	w_runIDstring.assign(runID.begin(), runID.end());
// 	const wchar_t* w_runID = w_runIDstring.c_str(); //allocated on stack no delete
// 	HWND hWnd = FindWindow(L"GLUT", w_runID); //get window handle
// 	//check to see if failed
// 	if (!hWnd)
// 	{
// 		return -1;
// 	}
// 	SetForegroundWindow(hWnd);


// 	RECT rect; //gets window boundaries
// 	BOOL success = GetWindowRect(hWnd, &rect);
// 	//check to see if failed
// 	if (!success)
// 	{
// 		return -1;
// 	}
	
// 	LONG width = rect.right - rect.left;
// 	LONG height = rect.bottom - rect.top;

// 	BYTE* pixels = new BYTE[3 * width*height];

// 	//Build the pixel array to be saved to PNG
// 	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	
// 	//Set up file name for screenshot
// 	string screenshotFileName = outputPath + runID + ".bmp";
// 	wstring w_screenshotFileNamestring;
// 	w_screenshotFileNamestring.assign(screenshotFileName.begin(), screenshotFileName.end());
// 	const WCHAR* w_screenshotFileName = w_screenshotFileNamestring.c_str();
	
// 	//Convert to bitmap and save
// 	long size;
// 	BYTE* buffer = ConvertRGBToBMPBuffer(pixels, width, height, &size);
// 	success = SaveBMP(buffer, width, height, size, w_screenshotFileName);
// 	if (!success)
// 	{
// 		code = -1;
// 		goto cleanup;
// 	}

// cleanup:
// 	delete[] buffer;
// 	delete[] pixels;

// 	return code;


// }

//http://tipsandtricks.runicsoft.com/Cpp/BitmapTutorial.html#chapter5
/*
BYTE* ConvertRGBToBMPBuffer ( BYTE* Buffer, int width,
int height, long* newsize )


This function takes as input an array of RGB values, it's width
and height.
The buffer gets then transformed to an array that can be used
to write to a windows bitmap file. The size of the array
is returned in newsize, the array itself is the
return value of the function.
Both input and output buffers must be deleted by the
calling function.

The input buffer is expected to consist of width * height
RGB triplets. Thus the total size of the buffer is taken as
width * height * 3.

The function then transforms this buffer so that it can be written
to a windows bitmap file:
First the RGB triplets are converted to BGR.
Then the buffer is swapped around since .bmps store
images uside-down.
Finally the buffer gets DWORD ( 32bit ) aligned,
meaning that each scanline ( 3 * width bytes ) gets
padded with 0x00 bytes up to the next DWORD boundary

*/

// BYTE* MeshEnv::ConvertRGBToBMPBuffer(BYTE* Buffer, int width, int height, long* newsize)
// {

// 	// first make sure the parameters are valid
// 	if ((nullptr == Buffer) || (width == 0) || (height == 0))
// 		return nullptr;

// 	// now we have to find with how many bytes
// 	// we have to pad for the next DWORD boundary	

// 	int padding = 0;
// 	int scanlinebytes = width * 3;
// 	while ((scanlinebytes + padding) % 4 != 0)     // DWORD = 4 bytes
// 		padding++;
// 	// get the padded scanline width
// 	int psw = scanlinebytes + padding;

// 	// we can already store the size of the new padded buffer
// 	*newsize = height * psw;

// 	// and create new buffer
// 	BYTE* newbuf = new BYTE[*newsize];

// 	// fill the buffer with zero bytes then we dont have to add
// 	// extra padding zero bytes later on
// 	memset(newbuf, 0, *newsize);

// 	// now we loop trough all bytes of the original buffer, 
// 	// swap the R and B bytes and the scanlines
// 	long bufpos;
// 	for (int y = 0; y < height; y++)
// 		for (int x = 0; x < 3 * width; x += 3)
// 		{
// 			bufpos = y * 3 * width + x;     // position in original buffer

// 			newbuf[bufpos] = Buffer[bufpos + 2];       // swap r and b
// 			newbuf[bufpos + 1] = Buffer[bufpos + 1]; // g stays
// 			newbuf[bufpos + 2] = Buffer[bufpos];     // swap b and r
// 		}

// 	return newbuf;
// }

// /*
// bool SaveBMP ( BYTE* Buffer, int width, int height,
// long paddedsize, LPCTSTR bmpfile )

// Function takes a buffer of size <paddedsize>
// and saves it as a <width> * <height> sized bitmap
// under the supplied filename.
// On error the return value is false.

// */

// bool MeshEnv::SaveBMP(BYTE* Buffer, int width, int height, long paddedsize, LPCTSTR bmpfile)
// {
// 	// declare bmp structures 
// 	BITMAPFILEHEADER bmfh;
// 	BITMAPINFOHEADER info;

// 	// andinitialize them to zero
// 	memset(&bmfh, 0, sizeof(BITMAPFILEHEADER));
// 	memset(&info, 0, sizeof(BITMAPINFOHEADER));

// 	// fill the fileheader with data
// 	bmfh.bfType = 0x4d42;       // 0x4d42 = 'BM'
// 	bmfh.bfReserved1 = 0;
// 	bmfh.bfReserved2 = 0;
// 	bmfh.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + paddedsize;
// 	bmfh.bfOffBits = 0x36;		// number of bytes to start of bitmap bits

// 	// fill the infoheader

// 	info.biSize = sizeof(BITMAPINFOHEADER);
// 	info.biWidth = width;
// 	info.biHeight = height;
// 	info.biPlanes = 1;			// we only have one bitplane
// 	info.biBitCount = 24;		// RGB mode is 24 bits
// 	info.biCompression = BI_RGB;
// 	info.biSizeImage = 0;		// can be 0 for 24 bit images
// 	info.biXPelsPerMeter = 0x0ec4;     // paint and PSP use this values
// 	info.biYPelsPerMeter = 0x0ec4;
// 	info.biClrUsed = 0;			// we are in RGB mode and have no palette
// 	info.biClrImportant = 0;    // all colors are important

// 	// now we open the file to write to
// 	HANDLE file = CreateFile(bmpfile, GENERIC_WRITE, FILE_SHARE_READ,
// 		nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
// 	if (file == nullptr)
// 	{
// 		CloseHandle(file);
// 		return false;
// 	}

// 	// write file header
// 	unsigned long bwritten;
// 	if (WriteFile(file, &bmfh, sizeof(BITMAPFILEHEADER), &bwritten, nullptr) == false)
// 	{
// 		CloseHandle(file);
// 		return false;
// 	}
// 	// write infoheader
// 	if (WriteFile(file, &info, sizeof(BITMAPINFOHEADER), &bwritten, nullptr) == false)
// 	{
// 		CloseHandle(file);
// 		return false;
// 	}
// 	// write image data
// 	if (WriteFile(file, Buffer, paddedsize, &bwritten, nullptr) == false)
// 	{
// 		CloseHandle(file);
// 		return false;
// 	}

// 	// and close file
// 	CloseHandle(file);

// 	return true;
// }

/**
Changes \ to / in strings. This is to fix file paths

@param path The path to be fixed
*/
string MeshEnv::fixPath(string &path)
{
	regex rgx("\\\\");
	return regex_replace(path, rgx, "/");
}

/**
Initialized the random number generator by providing a seed value from
the time.
*/
void MeshEnv::initRandomNumGen()
{
	//Initialize random number generation
	time_t seconds;
	time(&seconds); //assign time from clock
	//Seed the random number generator
	srand(static_cast<int>(seconds));
}

void MeshEnv::outputResults()
{

	// m_debugMode = 0x1800; //render objects
	// m_debugMode |= btIDebugDraw::DBG_NoDeactivation;
	// gDisableDeactivation = true;
	
	// //turn off context
	// m_debugMode |= btIDebugDraw::DBG_NoHelpText;

	// AHD: commented out, but it could be useful for drawing.
	// getDynamicsWorld()->getDebugDrawer()->setDebugMode(m_debugMode);

	// //set camera distance
	// setCameraDistance(sqrt(2 * pow(m_xdim, 2) + 2 * pow(m_zdim, 2) + pow(m_ydim*1.5, 2)));
	// m_cameraTargetPosition[0] = m_cameraTargetPosition[1] = m_cameraTargetPosition[2] = 0;
	// m_ele = 30;
	// m_azi = 40;

	// //Rendering function. Uses
	// renderme();
	// //APIENTRY functions -> cannot edit
	// glFlush();
	// swapBuffers();

	// takeScreenshot();
	// glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// screenshotHasBeenTaken = true;
	// m_debugMode = 1;
	// getDynamicsWorld()->getDebugDrawer()->setDebugMode(m_debugMode);
	// gDisableDeactivation = false; //allow deactivation again

	//iterate for file output
	for (auto itrTube = m_tubeList->begin(); itrTube != m_tubeList->end(); ++itrTube)
	{
		//create file with some filename
		int tubeNum = (*itrTube)->getTubeNum();
		stringstream filePathMaker;
		filePathMaker << outputPath << "CNT_Num_" << tubeNum << ".csv";
		string fileName = filePathMaker.str();
		ofstream file;
		file.open(fileName);
		// add initial info to file
		file << "Diameter:," << (*itrTube)->getDiameter() << "\nLength:," << (*itrTube)->getLength()
			<< "\nCylinder Height:," << (*itrTube)->getCylHeight() << "\nIntertube Spacing:," <<
			(*itrTube)->getMinSpacing()*2.0 << "\nIntercylinder Spacing:," << (*itrTube)->getTubeSpacing()
			<< "\nx" << tubeNum << ",y" << tubeNum << ",z" << tubeNum << endl;

		//The amount the constraint is shifted from the center of mass position from each cylinder
		btVector3 constraintShift = btVector3(0, -((*itrTube)->getCylHeight() + (*itrTube)->getTubeSpacing()) / 2.0, 0);

		//iterates over all of the cylinders of the current CNT. Outputs the positions to the file.
		auto tempCylList = (*itrTube)->getCylList();
		auto itrCyl = tempCylList->begin();

		//grab first cylinder information first and add to file
		btVector3 pos = (*itrCyl)->getCenterOfMassPosition();
		file << pos[0] << "," << pos[1] << "," << pos[2] << endl;
		++itrCyl;

		/*From the next cylinder, we can get both the center of mass positions and
		the position of the point to point constraint that connected this cylinder
		with the previous. Add all to file in order.*/
		for (itrCyl; itrCyl != tempCylList->end(); ++itrCyl)
		{
			btTransform currCyl = (*itrCyl)->getWorldTransform();
			btVector3 conPos = currCyl.operator*(constraintShift);
			file << conPos[0] << "," << conPos[1] << "," << conPos[2] << endl;
			pos = (*itrCyl)->getCenterOfMassPosition();
			file << pos[0] << "," << pos[1] << "," << pos[2] << endl;
		}
		file.close();
	}
	exitPhysics();
	std::cout << "Mesh created successfully!\n";
	exit(EXIT_SUCCESS);
}

bool compareRect(shared_ptr<rectangle> r1, shared_ptr<rectangle> r2)
{
	return r1->getHeight() > r2->getHeight();
}


// copied from DemoApplication.h
btRigidBody*	MeshEnv::localCreateRigidBody(float mass, const btTransform& startTransform,btCollisionShape* shape)
{
	btAssert((!shape || shape->getShapeType() != INVALID_SHAPE_PROXYTYPE));

	//rigidbody is dynamic if and only if mass is non zero, otherwise static
	bool isDynamic = (mass != 0.f);

	btVector3 localInertia(0,0,0);
	if (isDynamic)
		shape->calculateLocalInertia(mass,localInertia);

	//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects

// #define USE_MOTIONSTATE 1
#ifdef USE_MOTIONSTATE
	btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);

	btRigidBody::btRigidBodyConstructionInfo cInfo(mass,myMotionState,shape,localInertia);

	btRigidBody* body = new btRigidBody(cInfo);
	body->setContactProcessingThreshold(m_defaultContactProcessingThreshold);

#else
	btRigidBody* body = new btRigidBody(mass,0,shape,localInertia);	
	body->setWorldTransform(startTransform);
#endif//

	m_dynamicsWorld->addRigidBody(body);

	return body;
}

// copied from DemoApplication.h
btRigidBody*	MeshEnv::localCreateRigidBody(float mass, const btTransform& startTransform, btCollisionShape* shape, short group, short mask)
{
	btAssert((!shape || shape->getShapeType() != INVALID_SHAPE_PROXYTYPE));

	//rigidbody is dynamic if and only if mass is non zero, otherwise static
	bool isDynamic = (mass != 0.f);

	btVector3 localInertia(0, 0, 0);
	if (isDynamic)
		shape->calculateLocalInertia(mass, localInertia);

	//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects

// #define USE_MOTIONSTATE 1
#ifdef USE_MOTIONSTATE
	btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);

	btRigidBody::btRigidBodyConstructionInfo cInfo(mass, myMotionState, shape, localInertia);

	btRigidBody* body = new btRigidBody(cInfo);
	body->setContactProcessingThreshold(m_defaultContactProcessingThreshold);

#else
	btRigidBody* body = new btRigidBody(mass, 0, shape, localInertia);
	body->setWorldTransform(startTransform);
#endif//

	m_dynamicsWorld->addRigidBody(body,group,mask);

	return body;
}