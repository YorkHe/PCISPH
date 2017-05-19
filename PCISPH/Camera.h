#pragma once
#include <gl/glew.h>
#include <glm/gtc/matrix_transform.hpp>
#include "Engine.h"

class Camera
{
public:

	GLfloat Yaw;
	GLfloat Pitch;

	GLfloat lastX;
	GLfloat lastY;

	PCISPH::Vec3 Position;
	PCISPH::Vec3 Front;
	PCISPH::Vec3 Up;

	Camera(){}
	~Camera(){}

	glm::mat4 getViewMatrix()
	{
		return glm::lookAt(Position, Position + Front, Up);
	}
};

