#include "Simulator.h"

#include <iostream>
#include <algorithm>

Simulator::Simulator() :particleSet(), scene(nullptr), kernel()
{
}


Simulator::~Simulator()
{
}

int cmp(const void *a, const void *b) {
	return *(float*)a - *(float*)b;
}

void Simulator::init(const Scene *scene) {
	// Initialize scene
	this->scene = scene;

	// Initialize particleSet
	float particleMass = scene->referenceDensity * pow(2 * scene->particleRadius, 3);
	this->particleSet.init(particleMass);

	// Initialize kernel
	const int KERNEL_SCALE = 8;
	kernel.init(KERNEL_SCALE * scene->particleRadius);

	// Initialize grid
	this->grid.init(scene->boxSize, kernel.H);

	// Initialize particles' positions.
	const float step = 2 * scene->particleRadius;
	for (float x = scene->particleRadius; x < scene->fluidSize.x; x += step) {
		for (float y = scene->particleRadius; y < scene->fluidSize.y; y += step) {
			for (float z = scene->particleRadius; z < scene->fluidSize.z; z += step) {
				PCISPH::Vec3 pos = scene->fluidPosition + PCISPH::Vec3(x, y, z);

				this->grid.insert(pos, this->particleSet.count);

				particleSet.addParticle(ParticleSet::Particle(pos));
			}
		}
	}

	// Initialize particles' density
	this->computeDensity();

	/*std::vector<float> d = this->particleSet.density;
	std::qsort(&d[0], d.size(), sizeof(float), cmp);
	for (int i = 0; i < d.size(); i++) {
		std::cout << d[i] << std::endl;
	}*/
}

void Simulator::update() {
	computeForces();
}

/*A naive method to calculate density*/
void Simulator::computeDensity_test() {
	for (int i = 0; i < this->particleSet.count; i++) {
		float density = 0;
		PCISPH::Vec3 posi = this->particleSet.position[i];
		for (int j = 0; j < this->particleSet.count; j++) {
			if (i == j) continue;
			PCISPH::Vec3 r = posi - this->particleSet.position[j];
			density += this->particleSet.particleMass * kernel.poly6Kernel(r);
		}
		this->particleSet.density[i] = density;
		if (i % 1000 == 0) {
			std::cout << i << std::endl;
		}
	}
}

void Simulator::computeDensity() {
	PCISPH::iVec3 gridSize = this->grid.getSize();
	for (size_t x = 0; x < gridSize.x; x++) {
		for (size_t y = 0; y < gridSize.y; y++) {
			for (size_t z = 0; z < gridSize.z; z++) {
				PCISPH::iVec3 gridPos(x, y, z);
				Grid::GridUnit unit = this->grid.getGridUnit(gridPos);
				Grid::GridUnitSet neighbors;
				this->grid.getNeighborGridUnits(gridPos, neighbors);

				Grid::GridUnit::const_iterator unitIter;
				for (unitIter = unit.begin(); unitIter != unit.end(); unitIter++) {
					size_t particleIndex = *unitIter;
					float density = 0;
					Grid::GridUnitSet::const_iterator gridIter;
					for (gridIter = neighbors.begin(); gridIter != neighbors.end(); gridIter++) {
						Grid::GridUnit::const_iterator neighborIter;
						for (neighborIter = (*gridIter).begin(); neighborIter != (*gridIter).end(); neighborIter++) {
							size_t neighborIndex = *neighborIter;
							if (neighborIndex == particleIndex) {
								continue;
							}
							PCISPH::Vec3 r = this->particleSet.position[particleIndex] - this->particleSet.position[neighborIndex];
							density += this->particleSet.particleMass * kernel.poly6Kernel(r);
						}
					}
					this->particleSet.density[particleIndex] = density;
				}
			}
		}
	}
}

void Simulator::computeForces() {

}
