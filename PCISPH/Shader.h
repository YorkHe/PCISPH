#pragma once
#include <gl/glew.h>

class Shader
{
public:
	Shader();
	~Shader();

	GLuint program;
	Shader(const GLchar* vertexPath, const GLchar* fragmentPath);
	void use();
};

