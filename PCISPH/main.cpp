#include <gl/glew.h>

#include "System.h"

#include <iostream>
#include <stdexcept>
#include <glm/glm.hpp>

int main() {
	
	System system;

	try {
		system.init("./config.json");
		system.start();
		system.terminate();
	}
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
	}

	return 0;
}