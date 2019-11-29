#include <stdio.h>
#include <iostream>
#include <ctime>
#include <array>

#include "../misc_files/CommonInterfaces/CommonExampleInterface.h"
#include "../misc_files/CommonInterfaces/CommonGUIHelperInterface.h"
#include "../misc_files/Utils/b3Clock.h"

#include "../misc_files/ExampleBrowser/OpenGLGuiHelper.h"

#include "../lib/json.hpp"

#include "cnt_mesh.h"


// this block of code and the global variable and function is used for handling mouse input and
// moving objects via mouse. you can comment it if this capability is not needed any more.
//*************************************************************************************************
// CommonExampleInterface*    example;
cnt_mesh*    example;
int gSharedMemoryKey=-1;

//*************************************************************************************************

int main(int argc, char* argv[]) {

	// print the start time and start recording the run time
	std::clock_t start = std::clock();
	std::time_t start_time = std::time(nullptr);
	std::cout << std::endl << "start time:" << std::endl << std::asctime(std::localtime(&start_time)) << std::endl;


	// get the input JSON filename
	std::string filename;
	if (argc <= 1){
		filename = "input.json";
	} else {
		filename = argv[1];
	}

	// read the input JSON file
	std::ifstream input_file(filename.c_str());
	nlohmann::json j;
	input_file >> j;	
	
	int number_of_tubes_added_together = j["number of tubes added together"];
	int number_of_active_bundles = j["number of active bundles"];
	int number_of_tubes_before_deletion = j["number of tubes before deletion"];
	int number_of_unsaved_tubes = j["number of unsaved tubes"];
	int number_of_bundles = j["number of bundles"];
        int number_of_steps = j["number_of_steps"];
        btScalar time_step = j["time_step"];


	// flag to let the graphic visualization happen
	bool visualize = j["visualize"];

	// CommonExampleInterface* example;
	example = new cnt_mesh(NULL, j);	

	example->parse_json_prop();
	example->save_json_properties(j);

	
	example->initPhysics();
	example->create_container(); //container size is set in input.json
	example->create_tube_colShapes();


	// example->get_Ly();
	// example->add_tube_in_xz();

	int step_number = 0;
	while(true)
	{
		step_number++;
		btScalar dtSec = time_step;
		// btScalar dtSec = 0.01;
		example->stepSimulation(dtSec);

		if (step_number % number_of_steps == 0) // add new tubes every couple of steps.
		{	
			example->get_Ly();

			// add this many cnt's at a time
			for (int i=0; i<number_of_tubes_added_together; i++)
			{
				example->add_bundle_in_xz();
			}
			example->save_tubes(number_of_unsaved_tubes);
			example->freeze_bundles(number_of_active_bundles); // keep only this many of tubes active (for example 100) and freeze the rest of the tubes
			example->remove_tubes(number_of_tubes_before_deletion); // keep only this many of tubes in the simulation (for example 400) and delete the rest of objects
			
			std::cout << "number of saved tubes: " << example->no_of_saved_tubes() << ",  height [nm]:" << example->read_Ly() << "      \r" << std::flush;
			
			
		}
	if(example->no_of_saved_tubes()/7 > number_of_bundles)
		break;
	}

	
	// print the end time and the runtime
	std::clock_t end = std::clock();
	std::time_t end_time = std::time(nullptr);
	std::cout << std::endl << "end time:" << std::endl << std::asctime(std::localtime(&end_time));
	std::cout << "runtime: " << std::difftime(end_time,start_time) << " seconds" << std::endl << std::endl;
	
	

	example->exitPhysics();
	delete example;
	return 0;
}

