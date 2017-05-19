#include "System.h"



System::System()
{
	mWindow = new Window;
	mEngine = new Engine;
}


System::~System()
{
}

void System::init()
{

	// Init the window and engine;
	mWindow->init();
	mEngine->init();

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
