#pragma once

#include <gl/glew.h>
#include <vector>

#include "Scene.h"
#include "Types.h"

class ParticleSet
{
public:
	ParticleSet();
	~ParticleSet();

	void init(size_t cnt);

	std::vector<PCISPH::Vec3> position;
	std::vector<PCISPH::Vec3> velosity;
	std::vector<PCISPH::Vec3> forces;
	std::vector<PCISPH::Vec3> pressureForce;
	std::vector<float> density;
	std::vector<float> pressure;

	static constexpr float PARTICLE_MASS = 1.0f;

	size_t count;

private:
};

