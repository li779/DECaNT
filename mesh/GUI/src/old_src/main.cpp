#include <iostream>

#include <btBulletDynamicsCommon.h>

// #include "MeshEnv.h"


int main(int argc, char *argv[])
{

	std::cout << "Hello world!" << std::endl;

	// prepare the world here
	//***************************************************************************

	// Build the broadphase
	btBroadphaseInterface* broadphase = new btDbvtBroadphase();

	// Set up the collision configuration and dispatcher
	btDefaultCollisionConfiguration* collisionConfiguration = new btDefaultCollisionConfiguration();
	btCollisionDispatcher* dispatcher = new btCollisionDispatcher(collisionConfiguration);

	// The actual physics solver
	btSequentialImpulseConstraintSolver* solver = new btSequentialImpulseConstraintSolver;

	// The world.
	btDiscreteDynamicsWorld* dynamicsWorld = new btDiscreteDynamicsWorld(dispatcher, broadphase, solver, collisionConfiguration);
	dynamicsWorld->setGravity(btVector3(0, -10, 0));

	// Create the dynamic and static rigid bodies here
	//***************************************************************************

	// create the collision shapes
	btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(0, 1, 0), 1); // plane collision shape with an offset of 1 unit from the origin
	btCollisionShape* fallShape = new btSphereShape(1); // sphere collision shape with radius 1 meter

	// add the collision ground shape by creating a rigid body instance
	btDefaultMotionState* groundMotionState = new btDefaultMotionState(btTransform(btQuaternion(0, 0, 0, 1), btVector3(0, -1, 0))); // create the motion state for the ground plane (contains orientation and postion)
	btRigidBody::btRigidBodyConstructionInfo groundRigidBodyCI(0, groundMotionState, groundShape, btVector3(0, 0, 0)); // get information for constructing the rigid body
	btRigidBody* groundRigidBody = new btRigidBody(groundRigidBodyCI); // create the rigid body using the construction info
	dynamicsWorld->addRigidBody(groundRigidBody); // add the created rigid body to the world that we created before

	// add the collision sphere shape by creating a rigid body instance
	btDefaultMotionState* fallMotionState = new btDefaultMotionState(btTransform(btQuaternion(0, 0, 0, 1), btVector3(0, 50, 0)));
	btScalar mass = 1;
	btVector3 fallInertia(0, 0, 0);
	fallShape->calculateLocalInertia(mass, fallInertia); // calculate the inertia using a helper function
	btRigidBody::btRigidBodyConstructionInfo fallRigidBodyCI(mass, fallMotionState, fallShape, fallInertia);
	btRigidBody* fallRigidBody = new btRigidBody(fallRigidBodyCI);
	dynamicsWorld->addRigidBody(fallRigidBody);


	// Test debug drawing
	myDebugDrawer debugDrawer;
	// dynamicsWorld->getDebugDrawer()->setDebugMode(btIDebugDraw::DBG_DrawWireframe);


	// Stepping through the simulation here
	//***************************************************************************
	for (int i = 0 ; i < 300 ; i++)
	{

		dynamicsWorld->stepSimulation(1 / 60.f, 10);

		btTransform trans;
		fallRigidBody->getMotionState()->getWorldTransform(trans);

		std::cout << "sphere height: " << trans.getOrigin().getY() << std::endl;
	}

	// Clean up behind ourselves like good little programmers
	//***************************************************************************

	// delete the falling body
	dynamicsWorld->removeRigidBody(fallRigidBody);
	delete fallRigidBody->getMotionState();
	delete fallRigidBody;

	// delete the ground plane
	dynamicsWorld->removeRigidBody(groundRigidBody);
	delete groundRigidBody->getMotionState();
	delete groundRigidBody;

	// remember to delete collision shapes
	delete fallShape;
	delete groundShape;

	// delete the world, the solver, the dispatcher, collision configuration, and broadphase
	delete dynamicsWorld;
	delete solver;
	delete dispatcher;
	delete collisionConfiguration;
	delete broadphase;

	return 0;

	// if (argc == 1)
	// {
	//	int xmlArrayLength = 300;
	// 	char *inputXMLPathArray = new char[xmlArrayLength];
	// 	cout << "Enter config xml path (Example in program files directory):\n";
	// 	cin.getline(inputXMLPathArray, xmlArrayLength);
	// 	inputXMLPath = inputXMLPathArray;
	// 	delete[] inputXMLPathArray;
	// }
	// else if (argc > 2)
	// {
	// 	cout << "Incorrect parameters. Only enter file path of config xml\n";
	// 	system("pause");
	// 	exit(EXIT_FAILURE);
	// }
	// else
	// {
	// 	inputXMLPath = argv[1];
	// }
	
	// std::string inputXMLPath = "./CNT_Mesh_Config.xml";

	// inputXMLPath = MeshEnv::fixPath(inputXMLPath);

	// MeshEnv ccdDemo;
	// ccdDemo.inputXMLPath = inputXMLPath;
	// ccdDemo.initPhysics(btScalar(300));


// #ifdef CHECK_MEMORY_LEAKS
// 	ccdDemo.exitPhysics();
// #else
// 	// int vert = 0;
// 	// int horiz = 0;
// 	// GetDesktopResolution(vert, horiz);
// 	// return glutmain(argc, argv, vert, horiz, runID.c_str(), &ccdDemo);
// #endif

	//default glut doesn't return from mainloop
	return 0;
}