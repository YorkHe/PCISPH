#pragma once


#include "Types.h"

class DeviceParticleSetPointer
{
public:


	PCISPH::Vec3* position;
	PCISPH::Vec3* predictPosition;

	PCISPH::Vec3* velocity;
	PCISPH::Vec3* predictVelocity;

	PCISPH::Vec3* normal;
	PCISPH::Vec3* forces;
	PCISPH::Vec3* pressureForce;

	float* density;
	float* pressure;

	float particleMass;
	float maxDensityErr;
	size_t count;
};