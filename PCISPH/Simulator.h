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

	void update(const size_t maxIterations = 100);

	const ParticleSet& getParticleSet() const { return this->particleSet; }

private:

	ParticleSet particleSet;
	Grid fluidGrid;
	Grid boundaryGrid;
	const Scene *scene;
	Kernel kernel;

	// Boundary particles
	std::vector<float> boundaryMass;
	std::vector<float> boundaryDensity;
	std::vector<PCISPH::Vec3> boundaryPosition;

	float densityVarianceScale;


	void initParticlesPositions();
	void initDensityVarianceScale();
	void initBoundaryMass();
	void relax();

	void updateFluidGrid();
	void updateBoundaryGrid();

	void computeDensity();
	void computeNormal();
	void computeForces();
	void clearPressureAndPressureForce();
	void predictVelocityAndPosition();
	void updatePressure();
	void updatePressureForce();
	void updateVelocityAndPosition();

	/*update scenes which are variable. For example, add particles to the system in a FLOW scene*/
	void updateScene();
	void handleCollision();
};

