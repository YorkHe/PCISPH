#pragma once

#include "DeviceParticleSet.h"
#include "Scene.h"
#include "Kernel.h"
#include "Grid.h"
#include "Types.h"
#include "DeviceParticleSetPointer.h"

#include <vector>

class Simulator
{
public:
	Simulator();
	~Simulator();

	void init(const Scene *scene);

	void update(const size_t maxIterations = 100);

	DeviceParticleSet getParticleSet() const { return this->particleSet; }

private:

	DeviceParticleSet particleSet;
	thrust::device_ptr<DeviceParticleSetPointer> particleSetPointer;
	Grid fluidGrid;
	thrust::device_ptr<Grid> d_fluidGrid;
	Grid boundaryGrid;
	const Scene *scene;
	Kernel kernel;
	thrust::device_ptr<Kernel> d_kernel;

	thrust::device_vector<size_t> offset;
	thrust::host_vector<size_t> h_offset;

	// Boundary particles
	std::vector<float> boundaryMass;
	std::vector<float> boundaryDensity;
	std::vector<PCISPH::Vec3> boundaryPosition;

	float densityVarianceScale;


	void initParticlesPositions();
	void initDensityVarianceScale();
	void relax();

	void updateFluidGrid();

	void clearPressureAndPressureForce();

	/*update scenes which are variable. For example, add particles to the system in a FLOW scene*/
	void updateScene();
};

