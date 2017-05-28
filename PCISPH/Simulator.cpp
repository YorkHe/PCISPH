#include "Simulator.h"



Simulator::Simulator():particleSet(), scene(nullptr)
{
}


Simulator::~Simulator()
{
}

void Simulator::init(const Scene *scene) {
	this->scene = scene;
	this->particleSet.init(pow(scene->particleNumber, 3));

	//int len = pow(scene->particleNumber, 1.0 / 3);
	int len = scene->particleNumber;
	int count = 0;
	const float step = 0.02;

	// Initialize the particle set.
	for (int x = 0; x < len; x++) {
		for (int y = 0; y < len; y++) {
			for (int z = 0; z < len; z++) {
				this->particleSet.position[count] = scene->fluidBounds[0] + step * PCISPH::Vec3(x, y, z);
				count++;
			}
		}
	}
}

void Simulator::update() {

}

const std::vector<PCISPH::Vec3>& Simulator::getParticlePositions() const {
	return this->particleSet.position;
}
