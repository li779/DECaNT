// main file defining the monte carlo simulation

#include <stdio.h>
#include <chrono>
#include <ctime>
#include <omp.h>
#include <thread>

#include <experimental/filesystem>

#include "../lib/json.hpp"
#include "./monte_carlo/monte_carlo.h"
#include "./helper/utility.h"

#include "./exciton_transfer/cnt.h"
#include "./helper/prepare_directory.hpp"

int main(int argc, char *argv[]) {
  // set the number of threads for the parallel regions
  int n_threads = omp_get_max_threads();
  // int n_threads = std::min(omp_get_max_threads(), 15);
  omp_set_num_threads(n_threads);


  // print the start time and start recording the run time
  std::time_t start_time = std::time(nullptr);
  std::cout << "\n***\nstart time:\n" << std::asctime(std::localtime(&start_time)) << "***\n\n";

  // std::srand(std::time(0));
  std::srand(100);

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

	// get the json part related to exciton mc simulation
	if (j.count("exciton monte carlo")==0){
		throw std::invalid_argument("json input file does not contain \"exciton monte carlo\"");
	}
	nlohmann::json json_mc = j["exciton monte carlo"];

	// if exciton transfer type is davoody get cnt json information and add it to json_mc
	if (j["exciton monte carlo"]["rate type"].get<std::string>() == "davoody"){
		json_mc["cnts"] = j["cnts"];
	}

	//***********************************************************************************************
	// create monte carlo object and run the MC simulation
	//***********************************************************************************************

  assert(j.count("exciton monte carlo")>0);

  double time_step = json_mc["monte carlo time step"];
  int num_particle = json_mc["Number of particle in the simulation"];

  mc::monte_carlo sim(json_mc);

  /*sim.kubo_init();
  sim.save_json_properties();
  sim.kubo_create_particles();

  while (sim.time() < sim.kubo_max_time()) {
    sim.kubo_step(time_step);
    sim.kubo_save_avg_dispalcement_squared();

    std::cout << "kubo simulation: current time [seconds]: " << std::scientific << sim.time() << " .... "
              << "max time [seconds]: " << sim.kubo_max_time() << "\r" << std::flush;
  }

  std::cout << std::endl;
  std::cout << "Green-Kubo simulation finished!" << std::endl;
*/
 // std::exit(0);



	// initialize and run simulation for the exciton hopping
	sim.init();
	sim.save_json_properties();


  std::cout << "saving particle trajectories ..."<<std::endl;
  sim.calc_diffusion(time_step,num_particle);
  std::cout << "done!" << std::endl;


 /* std::cout << "\nrunning Monte Carlo:" << std::endl;

  while (true) {
    sim.step(time_step);
    sim.save_metrics(time_step);

    sim.repopulate_contacts();

    std::cout << "simulation time [seconds]: " << std::scientific << sim.time() << " .... "
              << "number of particles: " << sim.number_of_particles() << "\r" << std::flush;
  }

  // print the end time and the runtime
  std::time_t end_time = std::time(nullptr);
  std::cout << "\nruntime: " << std::difftime(end_time, start_time) << " seconds" << std::endl;
  std::cout << "\n***\nend time:\n" << std::asctime(std::localtime(&end_time)) << "***\n\n";
*/

  return 0;
}
