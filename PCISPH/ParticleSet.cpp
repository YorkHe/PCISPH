#include "ParticleSet.h"



ParticleSet::ParticleSet()
{
}


ParticleSet::~ParticleSet()
{
}

void ParticleSet::init(const float particleMass) {
	this->count = 0;

	this->particleMass = particleMass;

	this->position.clear();
	this->predictPosition.clear();

	this->velocity.clear();
	this->predictVelocity.clear();

	this->normal.clear();
	this->forces.clear();
	this->pressureForce.clear();

	this->density.clear();
	this->pressure.clear();

	this->maxDensityErr = 0.f;
}

void ParticleSet::addParticle(
	const PCISPH::Vec3 &p,
	const PCISPH::Vec3 &v
) {

	this->position.push_back(p);
	this->predictPosition.push_back(PCISPH::Vec3(0.f));

	this->velocity.push_back(v);
	this->predictVelocity.push_back(PCISPH::Vec3(0.f));

	this->normal.push_back(PCISPH::Vec3(0.f));
	this->forces.push_back(PCISPH::Vec3(0.f));
	this->pressureForce.push_back(PCISPH::Vec3(0.f));

	this->density.push_back(0.f);
	this->pressure.push_back(0.f);

	this->count += 1;
}