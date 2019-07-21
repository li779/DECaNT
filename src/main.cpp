#include <stdio.h>
#include <iostream>
#include <ctime>
#include <array>

#include "../misc_files/CommonInterfaces/CommonExampleInterface.h"
#include "../misc_files/CommonInterfaces/CommonGUIHelperInterface.h"
#include "../misc_files/Utils/b3Clock.h"

#include "../misc_files/OpenGLWindow/SimpleOpenGL3App.h"
#include "../misc_files/ExampleBrowser/OpenGLGuiHelper.h"

#include "../lib/json.hpp"

#include "cnt_mesh.h"


// this block of code and the global variable and function is used for handling mouse input and
// moving objects via mouse. you can comment it if this capability is not needed any more.
//*************************************************************************************************
// CommonExampleInterface*    example;
cnt_mesh*    example;
int gSharedMemoryKey=-1;

b3MouseMoveCallback prevMouseMoveCallback = 0;
static void OnMouseMove( float x, float y)
{
	bool handled = false; 
	handled = example->mouseMoveCallback(x,y); 	 
	if (!handled)
	{
		if (prevMouseMoveCallback)
			prevMouseMoveCallback (x,y);
	}
}

b3MouseButtonCallback prevMouseButtonCallback  = 0;
static void OnMouseDown(int button, int state, float x, float y) {
	bool handled = false;

	handled = example->mouseButtonCallback(button, state, x,y); 
	if (!handled)
	{
		if (prevMouseButtonCallback )
			prevMouseButtonCallback (button,state,x,y);
	}
}
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
	int number_of_active_tubes = j["number of active tubes"];
	int number_of_tubes_before_deletion = j["number of tubes before deletion"];
	int number_of_unsaved_tubes = j["number of unsaved tubes"];



	SimpleOpenGL3App* app;
	GUIHelperInterface* gui;

	// flag to let the graphic visualization happen
	bool visualize = j["visualize"];
	
	// SimpleOpenGL3App is a child of CommonGraphicsApp virtual class.
	app = new SimpleOpenGL3App("carbon nanotube mesh",1024,768,true);

	prevMouseButtonCallback = app->m_window->getMouseButtonCallback();
	prevMouseMoveCallback = app->m_window->getMouseMoveCallback();

	app->m_window->setMouseButtonCallback((b3MouseButtonCallback)OnMouseDown);
	app->m_window->setMouseMoveCallback((b3MouseMoveCallback)OnMouseMove);
	
	gui = new OpenGLGuiHelper(app,false); // the second argument is a dummy one
	// gui = new DummyGUIHelper();

	CommonExampleOptions options(gui);

	// CommonExampleInterface* example;
	example = new cnt_mesh(options.m_guiHelper, j);
	
	example->parse_json_prop();
	example->save_json_properties(j);


	example->initPhysics();
	example->create_container(); //container size is set in input.json
	example->create_tube_colShapes();

	if (visualize) {
		example->resetCamera();
	}
	
	int step_number = 0;


	// example->get_Ly();
	// example->add_tube_in_xz();


	while(true)
	{
		step_number ++;
	
		btScalar dtSec = 0.05;
		// btScalar dtSec = 0.01;
		example->stepSimulation(dtSec);

		if (step_number % 50 == 0) // add new tubes every couple of steps.
		{	
			example->get_Ly();

			// add this many cnt's at a time
			for (int i=0; i<number_of_tubes_added_together; i++)
			{
				example->add_tube_in_xz();
			}
			example->save_tubes(number_of_unsaved_tubes);
			example->freeze_tubes(number_of_active_tubes); // keep only this many of tubes active (for example 100) and freeze the rest of the tubes
			// example->remove_tubes(number_of_tubes_before_deletion); // keep only this many of tubes in the simulation (for example 400) and delete the rest of objects
			
			std::cout << "number of saved tubes: " << example->no_of_saved_tubes() << ",  height [nm]:" << example->read_Ly() << "      \r" << std::flush;
			
			if (visualize)
			{
				app->m_instancingRenderer->init();
				app->m_instancingRenderer->updateCamera(app->getUpAxis());
				example->renderScene();
				
				// draw some grids in the space
				DrawGridData dg;
				dg.upAxis = app->getUpAxis();
				app->drawGrid(dg);
				
				app->swapBuffer();

			}
		}

	}


	// if we did not visualize the simulation all along now visualize it one last time.
	if (not visualize)
	{
		example->resetCamera();
		app->m_instancingRenderer->init();
		app->m_instancingRenderer->updateCamera(app->getUpAxis());
		example->renderScene();
		
		// draw some grids in the space
		DrawGridData dg;
		dg.upAxis = app->getUpAxis();
		app->drawGrid(dg);
		
		app->swapBuffer();

	}
	
	// print the end time and the runtime
	std::clock_t end = std::clock();
	std::time_t end_time = std::time(nullptr);
	std::cout << std::endl << "end time:" << std::endl << std::asctime(std::localtime(&end_time));
	std::cout << "runtime: " << std::difftime(end_time,start_time) << " seconds" << std::endl << std::endl;
	
	std::cin.ignore();

	example->exitPhysics();
	delete example;
	delete app;


	return 0;
}

