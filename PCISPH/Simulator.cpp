#include "Simulator.h"
#include "utils.h"

#include <iostream>
#include <algorithm>

#define KERNEL_SCALE 4.0f
#define EPSILON 1e-30f

const bool DEBUG = false;

Simulator::Simulator() :particleSet(), scene(nullptr), kernel()
{
}


Simulator::~Simulator()
{
}

void Simulator::init(const Scene *scene) {
	// Initialize scene
	this->scene = scene;

	// Initialize particleSet
	float particleMass = scene->referenceDensity * pow(2 * scene->particleRadius, 3);
	this->particleSet.init(particleMass);
	std::cout << "paricleMass = " << particleMass << std::endl;

	// Initialize fluid particles' and boundary particles' positions.
	this->initParticlesPositions();

	// Initialize kernel
	float kernelRadius = KERNEL_SCALE * scene->particleRadius;
	kernel.init(kernelRadius);
	std::cout << "Kernel.H = " << kernel.H << std::endl;

	// Initialize grid
	this->fluidGrid.init(scene->boxSize, kernel.H);
	this->boundaryGrid.init(scene->boxSize, kernel.H);
	updateFluidGrid();
	updateBoundaryGrid();

	// Initialize density variance scale
	this->initDensityVarianceScale();
	std::cout << "densityVarianceScale = " << this->densityVarianceScale << std::endl;
	
	this->relax();

	std::cout << std::endl;
}

void Simulator::update(const size_t maxIterations) {

	updateFluidGrid();
	computeDensity();
	computeForces();
	clearPressureAndPressureForce();
	size_t iterations = 0;
	while (iterations < maxIterations) {
		particleSet.maxDensityErr = 0.f;
		predictVelocityAndPosition();
		updatePressure();
		updatePressureForce();
		static const int MIN_ITERATIONS = 3;
		static const float yita = 0.01;
		static const float TOL = yita * this->scene->referenceDensity;
		if (++iterations >= MIN_ITERATIONS && this->particleSet.maxDensityErr < TOL) {
			break;
		}
	}

	updateVelocityAndPosition();
	handleCollision();

	static size_t count = 0;
	std::cout << "maxDensityErr = " << this->particleSet.maxDensityErr << std::endl;
	std::cout << "iterations = " << iterations << std::endl;
	std::cout << "count = " << count << std::endl;
	std::cout << std::endl;
	count++;
}

void Simulator::updateFluidGrid() {
	this->fluidGrid.update(particleSet.position, [this](size_t i, size_t j) {
		std::swap(particleSet.position[i], particleSet.position[j]);
		std::swap(particleSet.velocity[i], particleSet.velocity[j]);
	});
}

void Simulator::updateBoundaryGrid() {
	this->boundaryGrid.update(boundaryPosition, [this](size_t i, size_t j) {
		std::swap(boundaryPosition[i], boundaryPosition[j]);
	});
}

void Simulator::computeDensity() {

	// compute boundary particles' densities
	for (size_t i = 0; i < boundaryDensity.size(); i++) {
		float fluidTerm = 0.0f;
		float boundaryTerm = 0.0f;

		fluidGrid.query(boundaryPosition[i], [this, i, &fluidTerm](size_t j) {
			PCISPH::Vec3 r = boundaryPosition[i] - particleSet.position[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < EPSILON) return;
			fluidTerm += kernel.poly6Kernel(r);
		});
		/*boundaryGrid.query(boundaryPosition[i], [this, i, &boundaryTerm](size_t j) {
			PCISPH::Vec3 r = boundaryPosition[i] - boundaryPosition[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < EPSILON) return;
			boundaryTerm += boundaryMass[j] * kernel.poly6Kernel(r);
		});*/

		boundaryDensity[i] = particleSet.particleMass * fluidTerm + boundaryTerm;
	}

	// compute fluid particles' densities
	for (size_t i = 0; i < particleSet.count; i++) {
		float fluidTerm = 0.0f;
		float boundaryTerm = 0.0f;

		PCISPH::Vec3 pos = particleSet.position[i];
		fluidGrid.query(pos, [this, &pos, &fluidTerm](size_t j) {
			PCISPH::Vec3 r = pos - particleSet.position[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < EPSILON) return;
			fluidTerm += kernel.poly6Kernel(r);
		});
		/*boundaryGrid.query(pos, [this, &pos, &boundaryTerm](size_t j) {
			PCISPH::Vec3 r = pos - boundaryPosition[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < EPSILON) return;
			boundaryTerm += boundaryMass[j] * kernel.poly6Kernel(r);
		});*/

		particleSet.density[i] = particleSet.particleMass * fluidTerm + boundaryTerm;
	}
}

void Simulator::computeNormal() {
	for (size_t i = 0; i < particleSet.count; i++) {
		PCISPH::Vec3 normal(0.0f);
		fluidGrid.query(particleSet.position[i], [this, i, &normal](size_t j) {
			PCISPH::Vec3 r = particleSet.position[i] - particleSet.position[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < EPSILON) return;
			normal += kernel.poly6KernelGradient(r) / particleSet.density[j];
		});
		normal *= kernel.H * particleSet.particleMass;
		particleSet.normal[i] = normal;
	}
}

void Simulator::computeForces() {
	computeNormal();
	float squaredMass = particleSet.particleMass * particleSet.particleMass;
	for (size_t i = 0; i < particleSet.count; i++) {
		PCISPH::Vec3 viscosity(0.f);
		PCISPH::Vec3 cohesion(0.f);
		PCISPH::Vec3 curvature(0.f);

		fluidGrid.query(particleSet.position[i], [this, i, &viscosity, &cohesion, &curvature](size_t j) {
			PCISPH::Vec3 r = particleSet.position[i] - particleSet.position[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < EPSILON) return;
			
			PCISPH::Vec3 vDiff = particleSet.velocity[i] - particleSet.velocity[j];
			viscosity -= vDiff * kernel.viscosityKernelLaplacian(r) / particleSet.density[j];
			
			float Kij = 2.0f * scene->referenceDensity / (particleSet.density[i] + particleSet.density[j]);
			cohesion += Kij * (r / rLen) * kernel.cohesionKernel(rLen);
			curvature += Kij * (particleSet.normal[i] - particleSet.normal[j]);
		});

		viscosity *= scene->viscosityCoefficient * squaredMass / particleSet.density[i];
		cohesion *= -scene->surfaceTensionCoefficient * squaredMass;
		curvature *= -scene->surfaceTensionCoefficient * particleSet.particleMass;

		particleSet.forces[i] = viscosity + cohesion + curvature + particleSet.particleMass * scene->gravity;
	}
}

void Simulator::clearPressureAndPressureForce() {
	memset(&(this->particleSet.pressure[0]), 0, sizeof(float) * particleSet.count);
	memset(&(this->particleSet.pressureForce[0]), 0, sizeof(PCISPH::Vec3) * particleSet.count);
}

void Simulator::predictVelocityAndPosition() {
	for (size_t i = 0; i < particleSet.count; i++) {
		PCISPH::Vec3 acceleration = (particleSet.forces[i] + particleSet.pressureForce[i]) / particleSet.particleMass;
		particleSet.predictVelocity[i] = particleSet.velocity[i] + scene->timeStep * acceleration;
		particleSet.predictPosition[i] = particleSet.position[i] + particleSet.predictVelocity[i] * scene->timeStep;
	}
}

void Simulator::updatePressure() {

	for (size_t i = 0; i < particleSet.count; i++) {
		float fDensity = 0.f;
		fluidGrid.query(particleSet.position[i], [this, i, &fDensity](size_t j) {
			PCISPH::Vec3 r = particleSet.predictPosition[i] - particleSet.predictPosition[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < EPSILON) return;
			fDensity += kernel.poly6Kernel(r);
		});
		float bDensity = 0.f;
		/*boundaryGrid.query(particleSet.position[i], [this, i, &bDensity](size_t j) {
			PCISPH::Vec3 r = particleSet.predictPosition[i] - boundaryPosition[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < EPSILON) return;
			bDensity += boundaryMass[j] * kernel.poly6Kernel(r);
		});*/
		float density = fDensity * particleSet.particleMass + bDensity;
		float densityVariation = std::max(0.f, density - scene->referenceDensity);
		particleSet.maxDensityErr = std::max(particleSet.maxDensityErr, densityVariation);

		particleSet.pressure[i] += this->densityVarianceScale * densityVariation;
	}
}

void Simulator::updatePressureForce() {
	for (size_t i = 0; i < particleSet.count; i++) {
		PCISPH::Vec3 pressureForce(0.f);
		fluidGrid.query(particleSet.position[i], [this, i, &pressureForce](size_t j) {
			PCISPH::Vec3 r = particleSet.predictPosition[i] - particleSet.predictPosition[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < 1e-5) return;
			float term1 = particleSet.pressure[i] / particleSet.density[i] / particleSet.density[i];
			float term2 = particleSet.pressure[j] / particleSet.density[j] / particleSet.density[j];
			pressureForce -= particleSet.particleMass * (term1 + term2) * kernel.spikyKernelGradient(r);
		});
		/*boundaryGrid.query(particleSet.position[i], [this, i, &pressureForce](size_t j) {
			PCISPH::Vec3 r = particleSet.predictPosition[i] - boundaryPosition[j];
			float rLen = PCISPH::length(r);
			if (rLen > kernel.H || rLen < 1e-5) return;
			float term1 = particleSet.pressure[i] / particleSet.density[i] / particleSet.density[i];
			float term2 = particleSet.pressure[i] / boundaryDensity[j] / boundaryDensity[j];
			pressureForce -= boundaryMass[j] * (term1 + term2) * kernel.spikyKernelGradient(r);
		});*/
		pressureForce *= particleSet.particleMass;
		particleSet.pressureForce[i] = pressureForce;
	}
}

void Simulator::updateVelocityAndPosition() {
	for (size_t i = 0; i < particleSet.count; i++) {
		PCISPH::Vec3 acceleration = (particleSet.forces[i] + particleSet.pressureForce[i]) / particleSet.particleMass;
		particleSet.velocity[i] += acceleration * scene->timeStep;
		particleSet.position[i] += particleSet.velocity[i] * scene->timeStep;
	}
}

void Simulator::updateScene() {
	// TODO
}

void Simulator::relax() {
	update(10000);
	for (PCISPH::Vec3 &v : particleSet.velocity) {
		v = scene->initVelocity;
	}
}

void Simulator::handleCollision() {
	static auto collisionFunc = [&](size_t i, const PCISPH::Vec3 &n, const float d) {
		particleSet.position[i] += d * n;
		particleSet.velocity[i] -= (1 + scene->restitution) * PCISPH::dot(particleSet.velocity[i], n) * n;
	};

	PCISPH::Vec3 box = scene->boxSize;
	for (size_t i = 0; i < particleSet.count; i++) {
		const auto &p = particleSet.position[i];
		if (p.x < 0) {
			collisionFunc(i, PCISPH::Vec3(1.f, 0.f, 0.f), -p.x);
		}
		if (p.x > box.x) {
			collisionFunc(i, PCISPH::Vec3(-1.f, 0.f, 0.f), p.x - box.x);
		}
		if (p.y < 0) {
			collisionFunc(i, PCISPH::Vec3(0.f, 1.f, 0.f), -p.y);
		}
		if (p.y > box.y) {
			collisionFunc(i, PCISPH::Vec3(0.f, -1.f, 0.f), p.y - box.y);
		}
		if (p.z < 0) {
			collisionFunc(i, PCISPH::Vec3(0.f, 0.f, 1.f), -p.z);
		}
		if (p.z > box.z) {
			collisionFunc(i, PCISPH::Vec3(0.f, 0.f, -1.f), p.z - box.z);
		}
	}
}

void Simulator::initParticlesPositions() {
	float r = scene->particleRadius;
	float d = 2 * r;

	// initialize boundary particles
	PCISPH::Vec3 box = scene->boxSize;

	// x-y plane
	for (float x = d; x <= box.x - d; x += d) {
		for (float y = d; y <= box.y - d; y += d) {
			boundaryPosition.emplace_back(x, y, 0.0f);
			boundaryPosition.emplace_back(x, y, box.z);
		}
	}
	// x-z plane
	for (float x = d; x <= box.x - d; x += d) {
		for (float z = d; z <= box.z - d; z += d) {
			boundaryPosition.emplace_back(x, 0.0f, z);
			boundaryPosition.emplace_back(x, box.y, z);
		}
	}
	// y-z plane
	for (float y = d; y <= box.y - d; y += d) {
		for (float z = d; z <= box.z - d; z += d) {
			boundaryPosition.emplace_back(0.0f, y, z);
			boundaryPosition.emplace_back(box.x, y, z);
		}
	}
	// x axis
	for (float x = d; x <= box.x - d; x += d) {
		boundaryPosition.emplace_back(x, 0.0f, 0.0f);
		boundaryPosition.emplace_back(x, 0.0f, box.z);
		boundaryPosition.emplace_back(x, box.y, 0.0f);
		boundaryPosition.emplace_back(x, box.y, box.z);
	}
	// y axis
	for (float y = d; y <= box.y - d; y += d) {
		boundaryPosition.emplace_back(0.0f, y, 0.0f);
		boundaryPosition.emplace_back(box.x, y, 0.0f);
		boundaryPosition.emplace_back(0.0f, y, box.z);
		boundaryPosition.emplace_back(box.x, y, box.z);
	}
	// z axis
	for (float z = d; z <= box.z - d; z += d) {
		boundaryPosition.emplace_back(0.0f, 0.0f, z);
		boundaryPosition.emplace_back(box.x, 0.0f, z);
		boundaryPosition.emplace_back(0.0f, box.y, z);
		boundaryPosition.emplace_back(box.x, box.y, z);
	}
	// corners
	boundaryPosition.emplace_back(0.0f, 0.0f, 0.0f);
	boundaryPosition.emplace_back(box.x, 0.0f, 0.0f);
	boundaryPosition.emplace_back(0.0f, box.y, 0.0f);
	boundaryPosition.emplace_back(0.0f, 0.0f, box.z);
	boundaryPosition.emplace_back(box.x, box.y, 0.0f);
	boundaryPosition.emplace_back(0.0f, box.y, box.z);
	boundaryPosition.emplace_back(box.x, 0.0f, box.z);
	boundaryPosition.emplace_back(box.x, box.y, box.z);

	boundaryMass.resize(boundaryPosition.size());
	boundaryDensity.resize(boundaryPosition.size());

	// initialize fluid particles
	PCISPH::Vec3 fBox = scene->fluidSize;
	PCISPH::Vec3 fPos = scene->fluidPosition;

	for (float x = r; x < fBox.x; x += d) {
		for (float y = r; y < fBox.y; y += d) {
			for (float z = r; z < fBox.z; z += d) {
				particleSet.addParticle(
					fPos + PCISPH::Vec3(x, y, z),				// position
					scene->initVelocity							// init velocity
				);
			}
		}
	}
}

void Simulator::initDensityVarianceScale() {
	float temp = scene->timeStep * particleSet.particleMass / scene->referenceDensity;
	float beta = 2.0f * temp * temp;
	PCISPH::Vec3 tmp1(0.0f);
	PCISPH::Vec3 tmp2(0.0f);
	float tmp3 = 0.0f;
	float r = scene->particleRadius;
	float d = 2.f * r;
	for (float x = -kernel.H - r; x <= kernel.H + r; x += d) {
		for (float y = -kernel.H - r; y <= kernel.H + r; y += d) {
			for (float z = -kernel.H - r; z <= kernel.H + r; z += d) {
				PCISPH::Vec3 pos(x, y, z);
				if (pos == PCISPH::Vec3(0.0f)) continue;

				PCISPH::Vec3 kernelValue1 = kernel.spikyKernelGradient(-pos);
				//PCISPH::Vec3 kernelValue1 = kernel.poly6KernelGradient(-pos);
				PCISPH::Vec3 kernelValue2 = kernel.poly6KernelGradient(-pos);
				tmp1 += kernelValue1;
				tmp2 += kernelValue2;
				tmp3 += PCISPH::dot(kernelValue1, kernelValue2);
			}
		}
	}
	this->densityVarianceScale = 1.0f / beta / (PCISPH::dot(tmp1, tmp2) + tmp3);
}

void Simulator::initBoundaryMass() {
	for (size_t i = 0; i < boundaryMass.size(); i++) {
		float weight = 0.0f;
		boundaryGrid.query(boundaryPosition[i], [this, i, &weight](size_t j) {
			weight += kernel.poly6Kernel(boundaryPosition[i] - boundaryPosition[j]);
		});
		boundaryMass[i] = scene->referenceDensity / weight / 1.17f;
	}
}