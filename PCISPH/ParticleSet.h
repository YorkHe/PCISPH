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

	float particleMass;

	void init(const float particleMass);

	const std::vector<PCISPH::Vec3>& getParticlePositions() const { return this->position; }

	void addParticle(
		const PCISPH::Vec3 &p,
		const PCISPH::Vec3 &v
	);

	std::vector<PCISPH::Vec3> position;
	std::vector<PCISPH::Vec3> predictPosition;

	std::vector<PCISPH::Vec3> velocity;
	std::vector<PCISPH::Vec3> predictVelocity;

	std::vector<PCISPH::Vec3> normal;
	std::vector<PCISPH::Vec3> forces;
	std::vector<PCISPH::Vec3> pressureForce;

	std::vector<float> density;
	std::vector<float> pressure;

	float maxDensityErr;

	size_t count;

private:

};

