#pragma once
#include <gl/glew.h>
#include <GLFW/glfw3.h>
class Window
{
public:
	Window();
	~Window();

	void init();

	static const int WIDTH = 1024;
	static const int HEIGHT = 768;

	GLFWwindow* mGLWindow;	
private:
};

