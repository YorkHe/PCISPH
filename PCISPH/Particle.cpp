#include "Particle.h"



Particle::Particle()
{
	position = PCISPH::Vec3(0, 0, 0);
	velocity = PCISPH::Vec3(0);
	density = PCISPH::Vec3(0);
	pressure = PCISPH::Vec3(0);
	extForce = PCISPH::Vec3(0);
	pressureForce = PCISPH::Vec3(0);
}


Particle::~Particle()
{
}

Particle::Particle(float x, float y, float z)
{
	Particle::Particle();
	position = PCISPH::Vec3(x, y, z);
}
