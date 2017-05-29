#pragma once

#include "ParticleSet.h"
#include "Scene.h"
#include "Types.h"
#include "Kernel.h"
#include "Grid.h"
#include <vector>

class Simulator
{
public:
	Simulator();
	~Simulator();

	void init(const Scene *scene);

	void update();

	const ParticleSet& getParticleSet() const { return this->particleSet; }

private:
	//typedef std::vector<std::vector<size_t>> Grid;

	ParticleSet particleSet;
	Grid grid;
	const Scene *scene;
	Kernel kernel;

	void computeDensity();
	void computeDensity_test();
	void computeForces();
};

