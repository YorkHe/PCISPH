#include "System.h"



System::System(): scene(nullptr), simulator(new Simulator), window(new Window), renderer(new Renderer)
{
}


System::~System()
{
	delete this->simulator;
	//delete this->window;  // Delete window will cause error, because of the ErrorCallback of GLFW
	delete this->renderer;
}

void System::init(const std::string &configPath) {
	try {
		// Initialize the new scene
		delete this->scene;
		this->scene = new Scene(configPath);

		// Refresh simulator and renderer
		this->simulator->init(this->scene);
		this->window->init();
		this->renderer->init(this->window, this->scene);
	}
	catch (std::exception e) {
		throw e;
	}
}

void System::start() {
	// Set OpenGL States
	glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_POINT_SMOOTH);

	// Begin drawing
	while (!glfwWindowShouldClose(window->mGLWindow))
	{
		glfwPollEvents();

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		this->renderer->draw(this->simulator->getParticleSet());
		this->simulator->update();

		glfwSwapBuffers(window->mGLWindow);
	}
}

void System::terminate() {
	glfwTerminate();
}