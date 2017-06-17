#pragma once

#include <gl/glew.h>
#include <vector>
#include <thrust/device_vector.h>

#include "Scene.h"
#include "Types.h"

class DeviceParticleSet
{
public:
	DeviceParticleSet();
	~DeviceParticleSet();

	float particleMass;

	void init(const float particleMass);

	void addParticle(
		const PCISPH::Vec3 &p,
		const PCISPH::Vec3 &v
	);

	thrust::device_vector<PCISPH::Vec3> position;
	thrust::host_vector<PCISPH::Vec3> h_position;
	thrust::device_vector<PCISPH::Vec3> predictPosition;

	thrust::device_vector<PCISPH::Vec3> velocity;
	thrust::host_vector<PCISPH::Vec3> h_velocity;

	thrust::device_vector<PCISPH::Vec3> predictVelocity;

	thrust::device_vector<PCISPH::Vec3> normal;
	thrust::device_vector<PCISPH::Vec3> forces;
	thrust::device_vector<PCISPH::Vec3> pressureForce;

	thrust::device_vector<float> density;
	thrust::device_vector<float> pressure;

	float maxDensityErr;

	size_t count;

private:

};

