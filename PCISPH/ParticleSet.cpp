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
	this->velosity.clear();
	this->forces.clear();
	this->pressureForce.clear();
	this->density.clear();
	this->pressure.clear();

	this->position.resize(0);
	this->velosity.resize(0);
	this->forces.resize(0);
	this->pressureForce.resize(0);
	this->density.resize(0);
	this->pressure.resize(0);
}

void ParticleSet::addParticle(const Particle &particle) {
	this->position.push_back(particle.position);
	this->velosity.push_back(particle. velosity);
	this->forces.push_back(particle.forces);
	this->pressureForce.push_back(particle.pressureForce);
	this->density.push_back(particle.density);
	this->pressure.push_back(particle.pressure);
	this->count += 1;
}
