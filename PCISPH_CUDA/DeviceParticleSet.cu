#include "DeviceParticleSet.h"



DeviceParticleSet::DeviceParticleSet()
{
}


DeviceParticleSet::~DeviceParticleSet()
{
}

void DeviceParticleSet::init(const float particleMass) {
	this->count = 0;

	this->particleMass = particleMass;

	this->maxDensityErr = 0.f;
}

void DeviceParticleSet::addParticle(
	const PCISPH::Vec3 &p,
	const PCISPH::Vec3 &v
) {

	this->h_position.push_back(p);

	this->h_velocity.push_back(v);

	this->count += 1;
}