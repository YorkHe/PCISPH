
#include <windows.h>
#include <iostream>
#include <stdexcept>


#include "System.h"

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