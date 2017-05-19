#include <iostream>
#include <exception>
#include "System.h"
#include <boost/program_options.hpp>

int main(int argc, char** argv) {

	System system;

	namespace po = boost::program_options;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "Produce help message")
		("config-file,c", po::value< std::vector<std::string>>(), "config file");

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

	try {
		system.init(vm);
		system.start();
		system.terminate();
	}
	catch (std::exception e) {
		std::cerr << "Runtime Error" << std::endl;
		std::cerr << e.what() << std::endl;
		getchar();
	}

	return 0;
}
