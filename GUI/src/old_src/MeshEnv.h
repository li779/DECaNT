/**
MeshEnv.h
Purpose: Header for MeshEnv.h

@author Alex Gabourie
@version 1.01 4/15/15

This is a heavily edited version of BasicDemo.h provided by Bullet Physics
*/

#ifndef MESH_ENV_H
#define MESH_ENV_H

// #include "GlutDemoApplication.h"
// #define PlatformDemoApplication GlutDemoApplication

#include <string>
#include <cstdint>


#include "BulletDynamics/Dynamics/btDynamicsWorld.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btMatrix3x3.h"
#include "LinearMath/btVector3.h"

#include "tube.h"

//preprocessor function define
#define BIT(x)	(1<<(x))

//use this to create collision maps such that the spacing cylinders do not collide with
// the spine cylinders at all. In the end, we only want spacing cylinders to collide with
// their non-adjacent couterparts.
enum collisionTypes
{
	COL_NOTHING = 0, //collide with nothing
	COL_SPINE = BIT(0), //collide with spine cylinders
	COL_SPACE = BIT(1), //collide with spacing cylinders
	COL_PLANE = BIT(2), //collides with plane
	COL_BOX = BIT(3), //collides with box
	COL_EVERYTHING = 0xFFFFFFFF
};


class btBroadphaseInterface;
class btCollisionShape;
class btOverlappingPairCache;
class btCollisionDispatcher;
class btConstraintSolver;
struct btCollisionAlgorithmCreateFunc;
class btDefaultCollisionConfiguration;

///MeshEnv takes the basic demo and works off of it.

/**
ALL FUNCTION COMMENTS ARE IN THE CPP FILE
*/

using namespace std;

// class MeshEnv : public PlatformDemoApplication
class MeshEnv
{

	//keep the collision shapes, for deletion/cleanup
	btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;
	btBroadphaseInterface*	m_broadphase;
	btCollisionDispatcher*	m_dispatcher;
	btConstraintSolver*	m_solver;
	btDefaultCollisionConfiguration* m_collisionConfiguration;
	btDynamicsWorld* m_dynamicsWorld;

	struct wave
	{
		int startIdx;
		int endIdx;
		double maxRadius;
		std::uint32_t cylCount;
	};


	shared_ptr<vector<shared_ptr<tube>>>	m_tubeList; 
	vector<wave> m_waveList;
	bool nextWave;
	int currWave;

	//Half-extent dimensions for the mesh area
	btScalar m_xdim;
	btScalar m_ydim;
	btScalar m_zdim;
	btScalar x_placement_extent;
	btScalar z_placement_extent;
	vector<double> placement_width;
	btScalar placementArea;
	int numTubes; //total number of tubes to be in simulation
	int tubesCreated = 0; //number of tubes created
	btScalar friction = 1.0;
	btScalar gravity = -9.81;
	btScalar minSpacing = 1.5;
	btScalar lmin = 100;
	btScalar lmax = 200;
	double curvMax = 2 * SIMD_PI / (500.0); // rads/angs
	uint8_t rotation = 0; // number from 0 to 3 for number of time 90 degrees is multiplied
	double meshYPeak = 0;
	double cylMass = 1; //arbitrary mass assignment

	uint8_t checkCntr = 0x00;
	uint8_t bitMask = 0x07;
	bool screenShotPrepped = false;
	bool screenshotHasBeenTaken = false;

public:

	MeshEnv()
	{
	}
	virtual ~MeshEnv()
	{
		exitPhysics();
	}
	void	initPhysics();
	void	initPhysics(float camDistance);

	float	exitPhysics();

	virtual void clientMoveAndDisplay();

	// virtual void displayCallback();
	virtual void	clientResetScene();

	// static DemoApplication* Create()
	// {
	// 	MeshEnv* demo = new MeshEnv;
	// 	demo->myinit();
	// 	demo->initPhysics();
	// 	return demo;
	// }

	///*Standard 3x3 matrix multiplication with a 3x1 vector. Needed as we have
	//types that are specific to bullet physics*/
	////@param a - 3x3 matrix
	////@param x - 3x1 vector
	////@return - a 3x1 vector that is a*x
	SIMD_FORCE_INLINE btVector3 multOperator(const btMatrix3x3& m, const btVector3& v);

	//rotates the object about the origin
	void rotateAboutOrigin(btMatrix3x3 &rotMatrix, btRigidBody* obj);

	//FIRST rotates the object about the origin  THEN shifts the object
	void rotateAndShift(btMatrix3x3 &rotMatrix, btRigidBody* obj, btVector3 &shift);

	//Shifts the object
	void shift(btRigidBody* obj, btVector3 &shift);

	//Newtons method to find tube separation
	//MeshEnv::tubeSepResult* 
	//	getTubeSeparation(btScalar height, btScalar diameter, btScalar curvature, btScalar guess);

	//COMMENT BLOCK NOT RELEVANT FOR BUNDLED TUBES
	/*
	A couple notes: I know the curvature of the CNT will give me angle/distance. The distance is 2*height+tubeSeparation.
	I also know that the tubeSeparation is radius*tan(angle/2). Solving these two equations for height, you get some
	expression like height = atan(tubeSeparation/radius)/curvature - tubeSeparation/2. If I'm inputting a tubeSeparation,
	everything is determined easily. However, I feel like it would be more useful to fix the length of the tube sections.
	From my calculations, I see that, for the relevant chiralities, the tube separations will be about 13-15 % of the tube
	height and not a big deal. A majority of what we see will be the tube sections. The diameters of the tubes also range from
	8.1 to 10.9 angstroms. From my MATLAB plots, I will be able to get a good initial guess for the tubeSeparation needed for
	the height I want. This is somewhat manual, but it was always going to be manual. So I choose height, get tube separation

	It also appears that guessing between -2 and 2 should allow for convergence for all of the relevant nanotubes
	*/

	std::string inputXMLPath;
	int xmlArrayLength;

	//converts the units of the xml doc to the units of the simulation
	double convertUnits(std::string unit, double val);

	int takeScreenshot();

	// bool SaveBMP(BYTE* Buffer, int width, int height, long paddedsize, LPCTSTR bmpfile);

	// BYTE* ConvertRGBToBMPBuffer(BYTE* Buffer, int width, int height, long* newsize);

	string static fixPath(string &path);

	void addWave();

	void initRandomNumGen();

	void outputResults();

	void tileWave();

	void makeTubesStatic(int start_tube, int end_tube);

	btRigidBody* localCreateRigidBody(float mass, const btTransform& startTransform,btCollisionShape* shape);
	btRigidBody* localCreateRigidBody(float mass, const btTransform& startTransform, btCollisionShape* shape, short group, short mask);

};

#endif //BASIC_DEMO_H