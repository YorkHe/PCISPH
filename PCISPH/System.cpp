#include "System.h"
#include <iostream>


System::System()
{
	mWindow = new Window;
	mEngine = new Engine;
	mScene = new Scene;
}


System::~System()
{
}

void System::init(boost::program_options::variables_map vm)
{

	//std::cout << vm["config-file"].as<std::vector<std::string>>()[0];
	mScene->parseFromFile(vm["config-file"].as<std::vector<std::string>>()[0]);

	// Init the window and engine;
	mWindow->init();
	mEngine->init(mScene);

}

void System::start()
{
	while (!glfwWindowShouldClose(mWindow->mGLWindow))
	{
		glfwPollEvents();



		glfwSwapBuffers(mWindow->mGLWindow);
	}
}

void System::terminate()
{
	glfwTerminate();
}
