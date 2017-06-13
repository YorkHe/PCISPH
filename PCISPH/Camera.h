#pragma once

#include "Types.h"
#include <glm/gtc/matrix_transform.hpp>

class Camera
{
public:
	Camera();
	~Camera();

	PCISPH::Vec3 Position = PCISPH::Vec3(0.0f, 0.0f, 0.8f); // test value
	PCISPH::Vec3 Front = -Position;
	PCISPH::Vec3 Up = PCISPH::Vec3(0.0f, 1.0f, 0.0f);

	glm::mat4 getProjViewMatrix() {
		return glm::perspective(45.0f, (float)4 / 3, 0.1f, 100.0f) * glm::lookAt(Position, Position + Front, Up);
	}
};

