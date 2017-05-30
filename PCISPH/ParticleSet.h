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

	typedef struct _Particle {
		PCISPH::Vec3 position;
		PCISPH::Vec3 velosity{ 0.0f, 0.0f, 0.0f };
		PCISPH::Vec3 forces{ 0.0f, 0.0f, 0.0f };
		PCISPH::Vec3 pressureForce{ 0.0f, 0.0f, 0.0f };
		float density = 0.0f;
		float pressure = 0.0f;

		_Particle(PCISPH::Vec3 pos) : position(pos) {}
	} Particle;

	float particleMass;

	void init(const float particleMass);

	const std::vector<PCISPH::Vec3>& getParticlePositions() const { return this->position; }
	void addParticle(const Particle &particle);

	std::vector<PCISPH::Vec3> position;
	std::vector<PCISPH::Vec3> velosity;
	std::vector<PCISPH::Vec3> forces;
	std::vector<PCISPH::Vec3> pressureForce;
	std::vector<float> density;
	std::vector<float> pressure;

	size_t count;

private:

};

