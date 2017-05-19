#include <iostream>
#include <exception>
#include "System.h"

int main(int argc, char** argv) {

	System system;

	try {
		system.init();
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
