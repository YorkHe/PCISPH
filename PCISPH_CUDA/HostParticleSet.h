#pragma once

#include <thrust/host_vector.h>

#include "DeviceParticleSet.h"
#include "Types.h"


class HostParticleSet
{
public:
	HostParticleSet();
	~HostParticleSet();

	void update(DeviceParticleSet*);

	float particleMass;
	size_t count;

	thrust::host_vector<PCISPH::Vec3> position;
	thrust::host_vector<float> density;
};

