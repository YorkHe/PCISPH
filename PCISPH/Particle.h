#pragma once
#include <vector>
#include "Scene.h"

class Particle
{
public:
	Particle();
	Particle(float x, float y, float z);
	~Particle();

	PCISPH::Vec3 position;
	PCISPH::Vec3 velocity;
	PCISPH::Vec3 density;
	PCISPH::Vec3 pressure;
	PCISPH::Vec3 extForce;
	PCISPH::Vec3 pressureForce;
};

