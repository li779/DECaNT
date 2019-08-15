#include <iostream>
#include <ctime>
#include <armadillo>

#include "cnt.h"
#include "exciton_transfer.h"
#include "constants.h"
#include "../../lib/json.hpp"

int main_exciton_transfer(int argc, char *argv[])
{
	std::time_t start_time = std::time(nullptr);
	std::cout << "\nstart time:" << std::endl << std::asctime(std::localtime(&start_time)) << std::endl;

	using json = nlohmann::json;

	std::string filename;
	if (argc <= 1){
		filename = "input.json";
	} else {
		filename = argv[1];
	}

	// read a JSON file
	std::ifstream input_file(filename.c_str());
	json j;
	input_file >> j;

	// get the parent directory for cnts
	std::string parent_directory = j["cnts"]["directory"];
	j["cnts"].erase("directory");

	// create excitons and calculate exciton dispersions
	std::vector<cnt> cnts;
	cnts.reserve(j["cnts"].size()); // this is reservation of space is crucial to ensure we do not move cnts, since the move constructor is not implemented yet
	for (const auto& j_cnt: j["cnts"])
	{
		cnts.emplace_back(cnt(j_cnt,parent_directory));
		cnts.back().calculate_exciton_dispersion();
	};

	// get the parent directory for cnts
	parent_directory = j["exciton transfer"]["directory"];
	j["exciton transfer"].erase("directory");

	for (const auto& j_ex_transfer:j["exciton transfer"])
	{
		exciton_transfer ex_transfer(j_ex_transfer, cnts, parent_directory);
		ex_transfer.run();
	}


	std::time_t end_time = std::time(nullptr);
	std::cout << std::endl << "end time:" << std::endl << std::asctime(std::localtime(&end_time));
	std::cout << "runtime: " << std::difftime(end_time,start_time) << " seconds" << std::endl << std::endl;

	return 0;
}
