#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <experimental/filesystem>
#include <fstream>
#include <cstddef>

#include "../lib/json.hpp"
#include "./helper/prepare_directory.hpp"
#include "cnt_mesh.h"

#include "btBulletDynamicsCommon.h"
#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h" 
#include "../misc_files/CommonInterfaces/CommonRigidBodyBase.h"

void cnt_mesh::initPhysics() {

	createEmptyDynamicsWorld();


	if (m_dynamicsWorld->getDebugDrawer())
		m_dynamicsWorld->getDebugDrawer()->setDebugMode(btIDebugDraw::DBG_DrawWireframe + btIDebugDraw::DBG_DrawContactPoints);
}

// create a rectangular container using half planes
void cnt_mesh::create_ground_plane() {

	Ly = 0.;

	// create a few basic rigid bodies
	//**********************************************************************************************

	// create the bottom ground plane normal to the y axis
	{
		btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(0, 1, 0), 0); // plane collision shape with an offset of 0 unit from the origin
		m_collisionShapes.push_back(groundShape);

		btTransform groundTransform;
		groundTransform.setIdentity();
		groundTransform.setOrigin(btVector3(0, 0, 0));

		btScalar mass(0.);
		createRigidBody(mass, groundTransform, groundShape, btVector4(0, 0, 1, 1)); // I think the last input is not used for anything. On paper it is supposed to be the collor
	}
}

void cnt_mesh::create_z_plane() {

	 // create the z direction side wall planes
	 {
	 	btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(0, 0, 1), 0); // plane collision shape with an offset of 0 unit from the origin
	 	m_collisionShapes.push_back(groundShape);

	 	btScalar mass(0.);

	 	btTransform groundTransform;
	 	groundTransform.setIdentity();
	 	groundTransform.setOrigin(btVector3(0,0,-_half_Lz));
	 	createRigidBody(mass,groundTransform,groundShape, btVector4(0,0,1,1)); // I think the last input is not used for anything. On paper it is supposed to be the collor
	 }

	 {
	 	btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(0, 0, -1), 0); // plane collision shape with an offset of 0 unit from the origin
	 	m_collisionShapes.push_back(groundShape);

	 	btScalar mass(0.);

	 	btTransform groundTransform;
	 	groundTransform.setIdentity();
	 	groundTransform.setOrigin(btVector3(0,0,_half_Lz));
	 	createRigidBody(mass,groundTransform,groundShape, btVector4(0,0,1,1)); // I think the last input is not used for anything. On paper it is supposed to be the collor
	 }


//	 // create the x direction side wall planes
//	 {
//	 	btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(1, 0, 0), 0); // plane collision shape with an offset of 0 unit from the origin
//	 	m_collisionShapes.push_back(groundShape);
//
//	 	btScalar mass(0.);
//
//	 	btTransform groundTransform;
//	 	groundTransform.setIdentity();
//	 	groundTransform.setOrigin(btVector3(-_half_Lx,0,0));
//	 	createRigidBody(mass,groundTransform,groundShape, btVector4(0,0,1,1)); // I think the last input is not used for anything. On paper it is supposed to be the collor
//	 }
//
//	 {
//	 	btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(-1, 0, 0), 0); // plane collision shape with an offset of 0 unit from the origin
//	 	m_collisionShapes.push_back(groundShape);
//
//	 	btScalar mass(0.);
//
//	 	btTransform groundTransform;
//	 	groundTransform.setIdentity();
//	 	groundTransform.setOrigin(btVector3(_half_Lx,0,0));
//	 	createRigidBody(mass,groundTransform,groundShape, btVector4(0,0,1,1)); // I think the last input is not used for anything. On paper it is supposed to be the collor
//	 }

	
}

// make tubes static in the simulation and only leave number_of_active_tubes as dynamic in the simulation.
void cnt_mesh::freeze_tubes(unsigned number_of_active_tubes) {
	if (tubes.size() <= number_of_active_tubes)
		return;

	auto it = std::prev(tubes.end(), number_of_active_tubes + 1);

	while (it->isDynamic) {
		// make the sections static by putting their mass equal to zero
		for (auto& b : it->bodies) {
			b->setMassProps(0, btVector3(0, 0, 0));
		}

		//delete constraints between the tube sections
		for (auto& c : it->constraints) {
			m_dynamicsWorld->removeConstraint(c);
			delete c;
			c = nullptr;
		}

		it->constraints.clear();
		it->isDynamic = false;

		if (it == tubes.begin())
			break;
		else
			--it;
	}
	
}

void cnt_mesh::freeze_bundles(unsigned number_of_active_bundles) {
	if (bundles.size() <= number_of_active_bundles)
		return;

	auto it = std::prev(bundles.end(), number_of_active_bundles + 1);
	while (it->isDynamic) {
		auto my_tube = std::prev(it->subtubes.end());
		while (my_tube->isDynamic) {
			// make the sections static by putting their mass equal to zero
			for (auto& b : my_tube->bodies) {
				b->setMassProps(0, btVector3(0, 0, 0));
			}

			//delete constraints between the tube sections
			for (auto& c : my_tube->constraints) {
				m_dynamicsWorld->removeConstraint(c);
				delete c;
				c = nullptr;
			}

			my_tube->constraints.clear();
			my_tube->isDynamic = false;

			if (my_tube == it->subtubes.begin())
				break;
			else
				--my_tube;
		}
		if (it == bundles.begin())
			break;
		else
			--it;
	}
	bundle reduce = bundles.front();
	while(!reduce.isDynamic) {
		bundles.pop_front();
		reduce = bundles.front();
	}
	return;
}


// remove the tubes from the simulation and only leave max_number_of_tubes in the simulation
void cnt_mesh::remove_tubes(unsigned max_number_of_tubes) {
	if (tubes.size() <= max_number_of_tubes)
		return;

	int n = tubes.size() - max_number_of_tubes;

	while (tubes.size() > max_number_of_tubes) {
		tube& my_tube = tubes.front();
		if (!my_tube.isDynamic) {
			for (auto& b : my_tube.bodies) {
				deleteRigidBody(b);
				b = nullptr;
			}
			my_tube.bodies.clear();
		}
		tubes.pop_front();
	}
}

// this method gives the appropriate coordinate for releasing the next tube
btVector3 cnt_mesh::drop_coordinate() {
	return btVector3(_half_Lx *((2.0 * float(std::rand()) / float(RAND_MAX)) - 1.0),
		drop_height, //drop_height + Ly,
		_half_Lz * ((2.0 * float(std::rand()) / float(RAND_MAX)) - 1.0)
	);
}

// this method gives the appropriate coordinate for releasing the next tube
btVector3 cnt_mesh::drop_para_coordinate() {
	return btVector3(_half_Lx,
					 drop_height, //drop_height + Ly,
					 _half_Lz * ((2.0 * float(std::rand()) / float(RAND_MAX)) - 1.0)
	);
}

void cnt_mesh::renderScene() {
	CommonRigidBodyBase::renderScene();
}

void cnt_mesh::resetCamera() {
	float dist = 41;
	float pitch = -35;
	float yaw = 52;
	float targetPos[3] = { 0,0.46,0 };
}

// save properties of the input tube
void cnt_mesh::save_one_tube(tube& t) {

	if (number_of_saved_tubes % 10000 == 0) {
		number_of_cnt_output_files++;

		position_file.close();
		std::string filename = "tube" + std::to_string(number_of_cnt_output_files) + ".pos.dat";
		output_file_path = _output_directory / filename;
		position_file.open(output_file_path, std::ios::out);
		position_file << std::showpos << std::scientific;

//		orientation_file.close();
//		filename = "tube" + std::to_string(number_of_cnt_output_files) + ".orient.dat";
//		output_file_path = _output_directory / filename;
//		orientation_file.open(output_file_path, std::ios::out);
//		orientation_file << std::showpos << std::scientific;

//		length_file.close();
//		filename = "tube" + std::to_string(number_of_cnt_output_files) + ".len.dat";
//		output_file_path = _output_directory / filename;
//		length_file.open(output_file_path, std::ios::out);
//		length_file << std::showpos << std::scientific;

		chirality_file.close();
		filename = "tube" + std::to_string(number_of_cnt_output_files) + ".chiral.dat";
		output_file_path = _output_directory / filename;
		chirality_file.open(output_file_path, std::ios::out);
		chirality_file << std::showpos << std::scientific;
	}
	number_of_saved_tubes++;
	position_file << "tube number: " << number_of_saved_tubes << " ; ";
//	orientation_file << "tube number: " << number_of_saved_tubes << " ; ";
	chirality_file << "tube number: " << number_of_saved_tubes << " ; ";

	int i = 0;
	btTransform trans;
	for (const auto& b : t.bodies) {
		b->getMotionState()->getWorldTransform(trans);
		position_file << trans.getOrigin().x() << " , " << trans.getOrigin().y() << " , " << trans.getOrigin().z() << " ; ";

//		btQuaternion qt = trans.getRotation();
//		btVector3 ax(0, 1, 0); // initial axis of the cylinder
//		ax = ax.rotate(qt.getAxis(), qt.getAngle());
//		orientation_file << ax.x() << " , " << ax.y() << " , " << ax.z() << " ; ";

//		length_file << t.body_length[i++] << ";";

		chirality_file << t.chirality[0] << " , "<< t.chirality[1]  << " ; ";
	}

	position_file << std::endl;
//	orientation_file << std::endl;
//	length_file << std::endl;
	chirality_file << std::endl;
}

void cnt_mesh::get_Ly() { //TODO the top is set to be the top of settled pile
	btTransform trans;
	float avgY = 0;
	int count = 0;
	for (const auto& t : tubes) {
		if (t.isDynamic) {
			for (const auto& b : t.bodies) {
				b->getMotionState()->getWorldTransform(trans);
				avgY += trans.getOrigin().getY();
				count++;
			}
		}
	}

	if (count == 0) {
		Ly = 0;
	}
	else {
		Ly = avgY / float(count);
	}

}

// set and save the json properties that is read and parsed from the input_json file.
void cnt_mesh::save_json_properties(nlohmann::json j) {
	std::ofstream json_file;
	json_file.open(_output_directory.path() / "input.json", std::ios::out);
	json_file << std::setw(4) << j << std::endl;
	json_file.close();
};

// add a tube with a diameter randomly picked from prebuilt colShapes for tube sections
void cnt_mesh::add_tube() {

	tubes.push_back(tube());
	tube& my_tube = tubes.back();

	int d = std::rand() % _tube_section_collision_shapes.size(); // index related to the diameter of the tube
	my_tube.diameter = _tube_diameter[d];

	int l = std::rand() % _tube_length.size(); // index related to the length of the tube
	float length = _tube_length[l];

	int c = std::rand() % _tube_chirality.size(); // index related to the chirality of the tube
	my_tube.chirality= _tube_chirality[c];

	btVector3 _drop_coordinate = drop_coordinate();

	// create a few dynamic rigidbodies
	//*********************************************************************************************


	btScalar mass(1.f);

	//rigidbody is dynamic if and only if mass is non zero, otherwise static
	bool isDynamic = (mass != 0.f);

	btVector3 localInertia(0, 0, 0);

	// Re-using the same collision is better for memory usage and performance
	btCollisionShape* colShape = nullptr;

	float c_length = 0;

	while (c_length < length) {
		int sl = std::rand() % _section_length.size();
		colShape = _tube_section_collision_shapes[d][sl];


		// Create Dynamic Objects
		btTransform startTransform;
		startTransform.setIdentity();

		// btScalar x_loc = _drop_coordinate[0];
		// btScalar y_loc = _drop_coordinate[1] + c_length + _section_length[sl]/2. + my_tube.diameter;
		// btScalar z_loc = _drop_coordinate[2];

		btScalar x_loc = 0;
		btScalar y_loc = 10;
		btScalar z_loc = c_length + _section_length[sl] / 2. + my_tube.diameter;

		btVector3 origin(x_loc, y_loc, z_loc);
		startTransform.setOrigin(origin);
		// startTransform.setRotation(btQuaternion(1, 1, 0, 0)); // set cylinder axis along x-direction
		startTransform.setRotation(btQuaternion(0, 0, 0, 1)); // set cylinder axis along y-direction
		// startTransform.setRotation(btQuaternion(0, 1, 1, 0)); // set cylinder axis along z-direction


		my_tube.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
		my_tube.bodies.back()->setMassProps(1.0, btVector3(1, 0, 1)); // turn off rotation along the y-axis of the cylinder shapes
		my_tube.body_length.push_back(_section_length[sl]);

		c_length += _section_length[sl] + my_tube.diameter;
	}

	my_tube.length = c_length;

	// //add N-1 spring constraints
	// for(int i=0;i<my_tube.bodies.size()-1;++i) {
	//   btRigidBody* b1 = my_tube.bodies[i];
	//   btRigidBody* b2 = my_tube.bodies[i+1];

	//   // spring constraint
	//   // btPoint2PointConstraint* centerSpring = new btPoint2PointConstraint(*b1, *b2, btVector3(0,(my_tube.body_length[i]+my_tube.diameter)/2,0), btVector3(0,-(my_tube.body_length[i+1]+my_tube.diameter)/2,0));
	//   // centerSpring->m_setting.m_damping = 1.5; //the damping value for the constraint controls how stiff the constraint is. The default value is 1.0
	//   // centerSpring->m_setting.m_impulseClamp = 0; //The m_impulseClamp value controls how quickly the dynamic rigid body comes to rest. The defual value is 0.0

	//   // const float pi = 3.14159265358979323846;

	//   // btTransform frameInA, frameInB;
	//   // frameInA = btTransform::getIdentity();
	//   // frameInA.getBasis().setEulerZYX(1, 0, 1);
	//   // frameInA.setOrigin(btVector3(0,1.1*(my_tube.body_length[i])/2,0));
	//   // frameInB = btTransform::getIdentity();
	//   // frameInB.getBasis().setEulerZYX(1,0, 1);
	//   // frameInB.setOrigin(btVector3(0,-1.1*(my_tube.body_length[i+1])/2,0));





	//   // btConeTwistConstraint* centerSpring = new btConeTwistConstraint(*b1, *b2, frameInA, frameInB);

	//   // centerSpring->setLimit(
	//   //                         pi/20, // _swingSpan1
	//   //                         pi/20, // _swingSpan2
	//   //                         pi, // _twistSpan
	//   //                         2, // _softness
	//   //                         0.3, // _biasFactor
	//   //                         1.0F // _relaxationFactor
	//   //                       );


	//   // m_dynamicsWorld->addConstraint(centerSpring);
	//   // my_tube.constraints.push_back(centerSpring);
	// }


	
}


// this method adds bundle in the xz plane
void cnt_mesh::add_bundle_in_xz(bool parallel) {

	const float pi = 3.14159265358979323846;

	tubes.push_back(tube());
	tube& my_tube1 = tubes.back();

	tubes.push_back(tube());
	tube& my_tube2 = tubes.back();

	tubes.push_back(tube());
	tube& my_tube3 = tubes.back();

	tubes.push_back(tube());
	tube& my_tube4 = tubes.back();

	tubes.push_back(tube());
	tube& my_tube5 = tubes.back();

	tubes.push_back(tube());
	tube& my_tube6 = tubes.back();

	tubes.push_back(tube());
	tube& my_tube7 = tubes.back();

	bundles.push_back(bundle());
	bundle& my_bundle = bundles.back();
	my_bundle.subtubes.push_back(my_tube1);
	my_bundle.subtubes.push_back(my_tube2);
	my_bundle.subtubes.push_back(my_tube3);
	my_bundle.subtubes.push_back(my_tube4);
	my_bundle.subtubes.push_back(my_tube5);
	my_bundle.subtubes.push_back(my_tube6);
	my_bundle.subtubes.push_back(my_tube7);

	int d = std::rand() % _tube_section_collision_shapes.size(); // index related to the diameter of the tube
  // We will add 7 tubes at a time to adhere to the hexagonally packed bundle morphology
  // Tube 1 will be the center tube, and tubes 2-7 will utilize a hinge constraint to attach to tube 1
	my_tube1.diameter = _tube_diameter[d];
	my_tube2.diameter = _tube_diameter[d];
	my_tube3.diameter = _tube_diameter[d];
	my_tube4.diameter = _tube_diameter[d];
	my_tube5.diameter = _tube_diameter[d];
	my_tube6.diameter = _tube_diameter[d];
	my_tube7.diameter = _tube_diameter[d];


	int l = std::rand() % _tube_length.size(); // index related to the length of the tube
	float length = _tube_length[l];

	//int c = std::rand() % _tube_chirality.size(); // index related to the chirality of the tube
	my_tube1.chirality= _tube_chirality[d];
	my_tube2.chirality= _tube_chirality[d];
	my_tube3.chirality= _tube_chirality[d];
	my_tube4.chirality= _tube_chirality[d];
	my_tube5.chirality= _tube_chirality[d];
	my_tube6.chirality= _tube_chirality[d];
	my_tube7.chirality= _tube_chirality[d];

	// set drop orientation of the tube
	float angle = parallel ? 0 : (float(std::rand() % 1000) / 1000. * pi);
	btVector3 ax(std::cos(angle), 0, std::sin(angle)); // axis vector for the tube sections

	// set a quaternion to determine the orientation of tube sections, note that the initial orientation of the tube sections are along the y-axis
	btQuaternion qt;
	btVector3 q_axis = ax.rotate(btVector3(0, 1, 0), pi / 2); // axis vector for the quaternion describing orientation of tube sections
	qt.setRotation(q_axis, pi / 2);

	btVector3 drop_coor = parallel ? drop_para_coordinate() : drop_coordinate();
	// btVector3 drop_coor(0,Ly,0);

	// set the density of the material making the tubes
	btScalar density = 1;

	// Re-using the same collision is better for memory usage and performance
	btCollisionShape* colShape = nullptr;

	float c_length = 0;

	while (c_length < length) {
		int sl = std::rand() % _section_length.size();
		btScalar sec_length_plus_distances = 1. * _section_length[sl];

		colShape = _tube_section_collision_shapes[d][sl];

		btScalar mass = density * _section_length[sl];

		btScalar ax_loc = c_length + sec_length_plus_distances / 2. - length / 2;

		btVector3 origin(ax_loc * ax + drop_coor);

		// create a btTransform that discribes the orientation and location of the rigidBody
		btTransform startTransform;
		startTransform.setOrigin(origin);
		startTransform.setRotation(qt);


		// create rigid bodies
		my_tube1.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
		my_tube2.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
		my_tube3.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
		my_tube4.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
		my_tube5.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
		my_tube6.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
		my_tube7.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
			// my_tube.bodies.back()->setMassProps(mass,btVector3(1,0,1)); // turn off rotation along the y-axis of the cylinder shapes
		my_tube1.body_length.push_back(sec_length_plus_distances);
		my_tube2.body_length.push_back(sec_length_plus_distances);
		my_tube3.body_length.push_back(sec_length_plus_distances);
		my_tube4.body_length.push_back(sec_length_plus_distances);
		my_tube5.body_length.push_back(sec_length_plus_distances);
		my_tube6.body_length.push_back(sec_length_plus_distances);
		my_tube7.body_length.push_back(sec_length_plus_distances);

		c_length += my_tube1.body_length.back();
	}

	my_tube1.length = c_length;
	my_tube2.length = c_length;
	my_tube3.length = c_length;

	//add N-1 constraints between the rigid bodies
	for (int i = 0; i < my_tube1.bodies.size() - 1; ++i) {
		btRigidBody* b1 = my_tube1.bodies[i];
		btRigidBody* b2 = my_tube1.bodies[i + 1];
		btRigidBody* b3 = my_tube2.bodies[i];
		btRigidBody* b4 = my_tube2.bodies[i + 1];
		btRigidBody* b5 = my_tube3.bodies[i];
		btRigidBody* b6 = my_tube3.bodies[i + 1];
		btRigidBody* b7 = my_tube4.bodies[i];
		btRigidBody* b8 = my_tube4.bodies[i + 1];
		btRigidBody* b9 = my_tube5.bodies[i];
		btRigidBody* b10 = my_tube5.bodies[i + 1];
		btRigidBody* b11 = my_tube6.bodies[i];
		btRigidBody* b12 = my_tube6.bodies[i + 1];
		btRigidBody* b13 = my_tube7.bodies[i];
		btRigidBody* b14 = my_tube7.bodies[i + 1];


		// // spring constraint
		// btPoint2PointConstraint* centerSpring = new btPoint2PointConstraint(*b1, *b2, btVector3(0,(my_tube.body_length[i])/2,0), btVector3(0,-(my_tube.body_length[i+1])/2,0));
		// centerSpring->m_setting.m_damping = 1.5; //the damping value for the constraint controls how stiff the constraint is. The default value is 1.0
		// centerSpring->m_setting.m_impulseClamp = 0; //The m_impulseClamp value controls how quickly the dynamic rigid body comes to rest. The defualt value is 0.0


		// cone constarint
		btTransform frameInA, frameInB;
		frameInA = btTransform::getIdentity();
		frameInA.getBasis().setEulerZYX(1, 0, 1);
		frameInA.setOrigin(btVector3(0, my_tube1.body_length[i] / 2, 0));
		frameInB = btTransform::getIdentity();
		frameInB.getBasis().setEulerZYX(1, 0, 1);
		frameInB.setOrigin(btVector3(0, -my_tube1.body_length[i + 1] / 2, 0));

		btConeTwistConstraint* centerSpring = new btConeTwistConstraint(*b1, *b2, frameInA, frameInB);
		btConeTwistConstraint* centerSpring2 = new btConeTwistConstraint(*b3, *b4, frameInA, frameInB);
		btConeTwistConstraint* centerSpring3 = new btConeTwistConstraint(*b5, *b6, frameInA, frameInB);
		btConeTwistConstraint* centerSpring4 = new btConeTwistConstraint(*b7, *b8, frameInA, frameInB);
		btConeTwistConstraint* centerSpring5 = new btConeTwistConstraint(*b9, *b10, frameInA, frameInB);
		btConeTwistConstraint* centerSpring6 = new btConeTwistConstraint(*b11, *b12, frameInA, frameInB);
		btConeTwistConstraint* centerSpring7 = new btConeTwistConstraint(*b13, *b14, frameInA, frameInB);
		centerSpring->setLimit(
			0, // _swingSpan1
			0, // _swingSpan2
			pi / 2, // _twistSpan
			1, // _softness
			0.3000000119F, // _biasFactor
			1.0F // _relaxationFactor
		);

		// hinge constarint
		frameInA = btTransform::getIdentity();
		frameInA.setOrigin(btVector3(my_tube1.diameter, 0, 0));
		frameInB = btTransform::getIdentity();
		frameInB.setOrigin(btVector3(0, 0, 0));

		//The hinge constraints connect each of the segments in tubes 2-7 to the corresponding
		//segment on tube 1

		btHingeConstraint* circleSpring = new btHingeConstraint(*b3, *b1, frameInA, frameInB);
		btHingeConstraint* circleSpring2 = new btHingeConstraint(*b5, *b1, frameInA, frameInB);
		btHingeConstraint* circleSpring3 = new btHingeConstraint(*b7, *b1, frameInA, frameInB);
		btHingeConstraint* circleSpring4 = new btHingeConstraint(*b9, *b1, frameInA, frameInB);
		btHingeConstraint* circleSpring5 = new btHingeConstraint(*b11, *b1, frameInA, frameInB);
		btHingeConstraint* circleSpring6 = new btHingeConstraint(*b13, *b1, frameInA, frameInB);
		btHingeConstraint* circleSpring7 = new btHingeConstraint(*b4, *b2, frameInA, frameInB);
		btHingeConstraint* circleSpring8 = new btHingeConstraint(*b6, *b2, frameInA, frameInB);
		btHingeConstraint* circleSpring9 = new btHingeConstraint(*b8, *b2, frameInA, frameInB);
		btHingeConstraint* circleSpring10 = new btHingeConstraint(*b10, *b2, frameInA, frameInB);
		btHingeConstraint* circleSpring11 = new btHingeConstraint(*b12, *b2, frameInA, frameInB);
		btHingeConstraint* circleSpring12 = new btHingeConstraint(*b14, *b2, frameInA, frameInB);


		m_dynamicsWorld->addConstraint(centerSpring, false);
		m_dynamicsWorld->addConstraint(centerSpring2, false);
		m_dynamicsWorld->addConstraint(centerSpring3, false);
		m_dynamicsWorld->addConstraint(centerSpring4, false);
		m_dynamicsWorld->addConstraint(centerSpring5, false);
		m_dynamicsWorld->addConstraint(centerSpring6, false);
		m_dynamicsWorld->addConstraint(centerSpring7, false);
		my_tube1.constraints.push_back(centerSpring);
		my_tube2.constraints.push_back(centerSpring2);
		my_tube3.constraints.push_back(centerSpring3);
		my_tube4.constraints.push_back(centerSpring4);
		my_tube5.constraints.push_back(centerSpring5);
		my_tube6.constraints.push_back(centerSpring6);
		my_tube7.constraints.push_back(centerSpring7);

		m_dynamicsWorld->addConstraint(circleSpring, false);
		m_dynamicsWorld->addConstraint(circleSpring2, false);
		m_dynamicsWorld->addConstraint(circleSpring3, false);
		m_dynamicsWorld->addConstraint(circleSpring4, false);
		m_dynamicsWorld->addConstraint(circleSpring5, false);
		m_dynamicsWorld->addConstraint(circleSpring6, false);
		m_dynamicsWorld->addConstraint(circleSpring7, false);
		m_dynamicsWorld->addConstraint(circleSpring8, false);
		m_dynamicsWorld->addConstraint(circleSpring9, false);
		m_dynamicsWorld->addConstraint(circleSpring10, false);
		m_dynamicsWorld->addConstraint(circleSpring11, false);
		m_dynamicsWorld->addConstraint(circleSpring12, false);

		//tube 1 has 12 hinge constraints because it is the central CNT
		my_tube1.constraints.push_back(circleSpring);
		my_tube1.constraints.push_back(circleSpring2);
		my_tube1.constraints.push_back(circleSpring3);
		my_tube1.constraints.push_back(circleSpring4);
		my_tube1.constraints.push_back(circleSpring5);
		my_tube1.constraints.push_back(circleSpring6);
		my_tube1.constraints.push_back(circleSpring7);
		my_tube1.constraints.push_back(circleSpring8);
		my_tube1.constraints.push_back(circleSpring9);
		my_tube1.constraints.push_back(circleSpring10);
		my_tube1.constraints.push_back(circleSpring11);
		my_tube1.constraints.push_back(circleSpring12);

		//All other tubes have just two hinge constraints

		my_tube2.constraints.push_back(circleSpring);
		my_tube2.constraints.push_back(circleSpring7);

		my_tube3.constraints.push_back(circleSpring2);
		my_tube3.constraints.push_back(circleSpring8);

		my_tube4.constraints.push_back(circleSpring3);
		my_tube4.constraints.push_back(circleSpring9);

		my_tube5.constraints.push_back(circleSpring4);
		my_tube5.constraints.push_back(circleSpring10);

		my_tube6.constraints.push_back(circleSpring5);
		my_tube6.constraints.push_back(circleSpring11);

		my_tube7.constraints.push_back(circleSpring6);
		my_tube7.constraints.push_back(circleSpring12);

	}




}


// TODO chirality is not added yet
// this method adds parallel tube like blinds in the xz plane
void cnt_mesh::add_single_tube_in_xz(bool parallel) {
	const float pi = 3.14159265358979323846;

	tubes.push_back(tube());
	tube& my_tube = tubes.back();

	int d = std::rand() % _tube_section_collision_shapes.size(); // index related to the diameter of the tube
	my_tube.diameter = _tube_diameter[d];
	my_tube.chirality = _tube_chirality[d];

	int l = std::rand() % _tube_length.size(); // index related to the length of the tube
	float length = _tube_length[l];

	// set drop orientation of the tube
	float angle = parallel ? 0 : (float(std::rand() % 1000) / 1000. * pi);
	btVector3 ax(std::cos(angle), 0, std::sin(angle)); // axis vector for the tube sections

	// set a quaternion to determine the orientation of tube sections, note that the initial orientation of the tube sections are along the y-axis
	btQuaternion qt;
	btVector3 q_axis = ax.rotate(btVector3(0, 1, 0), pi / 2); // axis vector for the quaternion describing orientation of tube sections
	qt.setRotation(q_axis, pi / 2);

	btVector3 drop_coor = parallel ? drop_para_coordinate() : drop_coordinate();
	// btVector3 drop_coor(0,Ly,0);

	// set the density of the material making the tubes
	btScalar density = 1;

	// Re-using the same collision is better for memory usage and performance
	btCollisionShape* colShape = nullptr;

	float c_length = 0;

	while (c_length < length) {
		int sl = std::rand() % _section_length.size();
		btScalar sec_length_plus_distances = 1. * _section_length[sl];

		colShape = _tube_section_collision_shapes[d][sl];

		btScalar mass = density * _section_length[sl];

		btScalar ax_loc = c_length + sec_length_plus_distances / 2. - length / 2;

		btVector3 origin(ax_loc * ax + drop_coor);

		// create a btTransform that discribes the orientation and location of the rigidBody
		btTransform startTransform;
		startTransform.setOrigin(origin);
		startTransform.setRotation(qt);


		// create rigid bodies
		my_tube.bodies.push_back(createRigidBody(mass, startTransform, colShape));	// no static object
		// my_tube.bodies.back()->setMassProps(mass,btVector3(1,0,1)); // turn off rotation along the y-axis of the cylinder shapes
		my_tube.body_length.push_back(sec_length_plus_distances);

		c_length += my_tube.body_length.back();
	}

	my_tube.length = c_length;

	//add N-1 constraints between the rigid bodies
	for (int i = 0; i < my_tube.bodies.size() - 1; ++i) {
		btRigidBody* b1 = my_tube.bodies[i];
		btRigidBody* b2 = my_tube.bodies[i + 1];

		// // spring constraint
		// btPoint2PointConstraint* centerSpring = new btPoint2PointConstraint(*b1, *b2, btVector3(0,(my_tube.body_length[i])/2,0), btVector3(0,-(my_tube.body_length[i+1])/2,0));
		// centerSpring->m_setting.m_damping = 1.5; //the damping value for the constraint controls how stiff the constraint is. The default value is 1.0
		// centerSpring->m_setting.m_impulseClamp = 0; //The m_impulseClamp value controls how quickly the dynamic rigid body comes to rest. The defualt value is 0.0


		// cone constarint
		btTransform frameInA, frameInB;
		frameInA = btTransform::getIdentity();
		frameInA.getBasis().setEulerZYX(1, 0, 1);
		frameInA.setOrigin(btVector3(0, my_tube.body_length[i] / 2, 0));
		frameInB = btTransform::getIdentity();
		frameInB.getBasis().setEulerZYX(1, 0, 1);
		frameInB.setOrigin(btVector3(0, -my_tube.body_length[i + 1] / 2, 0));

		btConeTwistConstraint* centerSpring = new btConeTwistConstraint(*b1, *b2, frameInA, frameInB);
		centerSpring->setLimit(
			0, // _swingSpan1
			0, // _swingSpan2
			pi / 2, // _twistSpan
			1, // _softness
			0.3000000119F, // _biasFactor
			1.0F // _relaxationFactor
		);


		m_dynamicsWorld->addConstraint(centerSpring, false);
		my_tube.constraints.push_back(centerSpring);
	}



}

// make tubes static in the simulation and only leave number_of_active_tubes as dynamic in the simulation.
void cnt_mesh::save_tubes(int number_of_unsaved_tubes) {
	if (tubes.size() <= number_of_unsaved_tubes)
		return;

	auto it = prev(tubes.end(), number_of_unsaved_tubes);

	int n = 0;
	while (!it->isSaved) {
		save_one_tube(*it);
		it->isSaved = true;
		it--;
		if (it == tubes.begin())
			break;
	}
}
