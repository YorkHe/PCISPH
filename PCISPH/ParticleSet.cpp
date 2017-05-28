#include "ParticleSet.h"



ParticleSet::ParticleSet()
{
}


ParticleSet::~ParticleSet()
{
}

void ParticleSet::init(size_t cnt) {
	this->count = cnt;

	this->position.clear();
	this->velosity.clear();
	this->forces.clear();
	this->pressureForce.clear();
	this->density.clear();
	this->pressure.clear();

	this->position.resize(cnt);
	this->velosity.resize(cnt);
	this->forces.resize(cnt);
	this->pressureForce.resize(cnt);
	this->density.resize(cnt);
	this->pressure.resize(cnt);
}
