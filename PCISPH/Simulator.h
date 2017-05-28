#pragma once

#include "ParticleSet.h"
#include "Scene.h"
#include "Types.h"
#include <vector>

class Simulator
{
public:
	Simulator();
	~Simulator();

	void init(const Scene *scene);

	void update();

	const std::vector<PCISPH::Vec3>& getParticlePositions() const;

private:
	ParticleSet particleSet;
	const Scene *scene;
};

