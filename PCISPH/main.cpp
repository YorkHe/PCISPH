#include <gl/glew.h>

#include "System.h"
#include <windows.h>
#include <iostream>
#include <stdexcept>
#include <glm/glm.hpp>

extern "C" {
	__declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
}

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