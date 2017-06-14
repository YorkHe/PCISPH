#include "Window.h"
#include <stdexcept>


Window::Window():mGLWindow(nullptr)
{
}

Window::~Window()
{
	if (mGLWindow != nullptr) {
		glfwDestroyWindow(mGLWindow);
	}
}

void Window::init() {

	if (this->mGLWindow != nullptr) {
		return;
	}

	glfwSetErrorCallback([](int error, const char* desc)
	{
		throw std::runtime_error(std::string("glfwError ") + desc);
	});

	if (glfwInit() != GL_TRUE)
	{
		throw std::runtime_error("GLFW initialization error");
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);


	mGLWindow = glfwCreateWindow(WIDTH, HEIGHT, "PCISPH", nullptr, nullptr);

	if (!mGLWindow)
	{
		throw std::runtime_error("GLFW creating window error");
	}

	glfwMakeContextCurrent(mGLWindow);

	glewExperimental = GL_TRUE;

	if (glewInit() != GLEW_OK)
	{
		throw std::runtime_error("GLEW initialization error");
	}

	glViewport(0, 0, WIDTH, HEIGHT);
}

