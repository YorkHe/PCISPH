#include "HostParticleSet.h"



HostParticleSet::HostParticleSet()
{
}


HostParticleSet::~HostParticleSet()
{
}

void HostParticleSet::update(DeviceParticleSet* deviceParticleSet)
{
	particleMass = deviceParticleSet->particleMass;
	count = deviceParticleSet->count;
	position = deviceParticleSet->position;
	density = deviceParticleSet->density;
}
