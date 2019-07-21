#ifndef cnt_mesh_h
#define cnt_mesh_h

#include <cstdlib>
#include <ctime>
#include <vector>
#include <array>
#include <list>
#include <experimental/filesystem>
#include <fstream>

#include "btBulletDynamicsCommon.h"
#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "../lib/json.hpp"
#include "./helper/prepare_directory.hpp"

#include "../misc_files/CommonInterfaces/CommonRigidBodyBase.h"

struct cnt_mesh : public CommonRigidBodyBase
{

	private:

	// infomation storing the simulation information
	std::experimental::filesystem::directory_entry _output_directory; // this is the address of the output directory
	std::experimental::filesystem::path output_file_path; // this is the full address of the output file.
	std::fstream position_file; // this is the output file that the coordinate of cnt sections are written into.
	std::fstream orientation_file; // this is the output file that the orientation of cnt sections are written into.
	std::fstream length_file; // this is the output file that the length of cnt sections are written into.
	int number_of_saved_tubes; // this is the total number of cnts whos coordinates are saved into output file.
	int number_of_cnt_output_files; // this is the number of output files that the cnt coordinates has been written into.

	nlohmann::json _json_prop; // json object containing simulation input properties

	// container properties
	float _half_Lx,_half_Lz; // container size. y-axis is the direction in which the top of the container is open (vertical direction and the direction in which gravity is applied), Lx and Lz are the direction in which the container is enclosed
	float Ly; // this is an average height for the stack of cnt mesh

	float drop_height=0;

	std::vector<float> _tube_diameter;
	std::vector<float> _section_length;
	std::vector<float> _tube_length;
	
	std::vector<std::vector<btCollisionShape*>> _tube_section_collision_shapes; // first index determines the diameter, the second index determines the length of the section

	// class to store information and the rigid bodies of each separate cnt.
	struct tube {
		int number_of_sections;
		float diameter=0; // diameter of the tube which is the same for all body objects
		float length=0;
		bool isDynamic=true;
		bool isSaved=false;
		std::vector<btRigidBody*> bodies; // btRigidBody objects that make the tube
		std::vector<float> body_length;
		std::vector<btTypedConstraint*> constraints; // movement constraints that connect the bodies
	};
	// list to store all the tubes that we will in the simulation
	std::list<tube> tubes;

	btVector3 drop_coordinate(); // this method gives the appropriate coordinate for releasing the next tube

  public:
	// constructor
	cnt_mesh(struct GUIHelperInterface* helper, nlohmann::json j): CommonRigidBodyBase(helper) {
		std::srand(std::time(0)); // use current time as seed for random generator
		std::cout << "seeded the random number generator!!!" << std::endl;

		_json_prop = j;

		// initialize the output parameters
		number_of_saved_tubes = 0;
		number_of_cnt_output_files = 0;
	}

	// set the simulation properties according to _json_prop object which is constructed from input.json
	void parse_json_prop(){
		
		std::string output_path = _json_prop["output directory"];
		bool keep_old_files=_json_prop["keep old files"];
		_output_directory = prepare_directory(output_path, keep_old_files);

		float container_half_width = float(_json_prop["container width [nm]"])/2.;
		_half_Lx = container_half_width;
		_half_Lz = container_half_width;
		
		_tube_diameter.push_back(float(_json_prop["cnt diameter [nm]"]));

		int min_tube_length = _json_prop["cnt total length [nm]"][0];
		int max_tube_length = _json_prop["cnt total length [nm]"][1];
		for (int i=min_tube_length; i<=max_tube_length; ++i){
			_tube_length.push_back(float(i));
		}

		int min_section_length = _json_prop["cnt section length [nm]"][0];
		int max_section_length = _json_prop["cnt section length [nm]"][1];
		for (int i=min_section_length; i<=max_section_length; ++i){
			_section_length.push_back(float(i));
		}

		drop_height = float(_json_prop["drop height [nm]"]);
	}

	// create all the btCollisionShape that are used to make tubes
	void create_tube_colShapes(){
		btCollisionShape* colShape=nullptr;
		for (float d: _tube_diameter){
			_tube_section_collision_shapes.push_back(std::vector<btCollisionShape*>());
			for (float l: _section_length){
				colShape = new btCylinderShape(btVector3(d/2.0 ,l/2.0, d/2.0));
				m_collisionShapes.push_back(colShape);
				_tube_section_collision_shapes.back().push_back(colShape);
			}
		}

		std::cout << "new collision shape created!" << std::endl;
	}

	void initPhysics();
	void renderScene();
	
	void resetCamera();

	inline void stepSimulation(float deltaTime) {
		if (m_dynamicsWorld)
		{
			m_dynamicsWorld->stepSimulation(deltaTime,10,deltaTime);
		}
	}

	// set and save the json properties that is read and parsed from the input_json file.
	void save_json_properties(nlohmann::json j);

	// this method adds a tube to the system.
	void add_tube();

	// this method adds a tube in the xz plane
	void add_tube_in_xz();

	// this method creates an open top container for the cnts
	void create_container();

	// gets the number of tubes in the simulation
	inline int num_tubes() {
		return tubes.size();
	}

	// make tubes static in the simulation and only leave _number_of_active_tubes as dynamic in the simulation.
	void freeze_tubes(unsigned number_of_active_tubes);

	// remove the tubes from the simulation and only leave _max_number_of_tubes in the simulation
	void remove_tubes(unsigned max_number_of_tubes);

	// save the coordinates of the tube to an output file.
	void save_one_tube(tube &t);

	// update Ly, which is roughly the height of the filled container
	void get_Ly();

	// read Ly, which is roughly the height of the filled container
	inline const float& read_Ly() {
		return Ly;
	};
	
	// get number of saved tubes
	inline const int& no_of_saved_tubes() {
		return number_of_saved_tubes;
	};

	void save_tubes(int number_of_unsaved_tubes);

};


#endif //cnt_mesh_h
