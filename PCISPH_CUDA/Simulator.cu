#include "Simulator.h"
#include "utils.h"

#include <iostream>
#include <cuda_runtime.h>

#include <thrust/device_new.h>

#define KERNEL_SCALE 4.0f
#define EPSILON 1e-30f

const bool DEBUG = false;

Simulator::Simulator() :particleSet(), scene(nullptr), kernel()
{
}


Simulator::~Simulator()
{
}

__global__ void initParticleSetPointer(
	DeviceParticleSetPointer* p,
	PCISPH::Vec3* position,
	PCISPH::Vec3* predictPosition,
	PCISPH::Vec3* velocity,
	PCISPH::Vec3* predictVelocity,
	PCISPH::Vec3* normal,
	PCISPH::Vec3* forces,
	PCISPH::Vec3* pressureForce,
	float* density,
	float* pressure,
	float particleMass,
	float maxDensityErr,
	size_t count
)
{
	p->position = position;
	p->predictPosition = predictPosition;
	p->velocity = velocity;
	p->predictVelocity = predictVelocity;
	p->normal = normal;
	p->forces = forces;
	p->pressureForce = pressureForce;
	p->density = density;
	p->pressure = pressure;
	p->particleMass = particleMass;
	p->maxDensityErr = maxDensityErr;
	p->count = count;
}

__global__ void initDeviceGrid(
	Grid* grid,
	float cellSize,
	PCISPH::Vec3 boxSize,
	PCISPH::uVec3 gridSize,
	size_t cellNumber
)
{
	grid->cellSize = cellSize;
	grid->boxSize = boxSize;
	grid->gridSize = gridSize;
	grid->cellNumber = cellNumber;
}

__global__
void initDeviceKernel(Kernel* k, float h)
{
	k->init(h);
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

	particleSetPointer = thrust::device_new<DeviceParticleSetPointer>();

	initParticleSetPointer << <1, 1 >> >(
		thrust::raw_pointer_cast(particleSetPointer),
		thrust::raw_pointer_cast(particleSet.position.data()),
		thrust::raw_pointer_cast(particleSet.predictPosition.data()),
		thrust::raw_pointer_cast(particleSet.velocity.data()),
		thrust::raw_pointer_cast(particleSet.predictVelocity.data()),
		thrust::raw_pointer_cast(particleSet.normal.data()),
		thrust::raw_pointer_cast(particleSet.forces.data()),
		thrust::raw_pointer_cast(particleSet.pressureForce.data()),
		thrust::raw_pointer_cast(particleSet.density.data()),
		thrust::raw_pointer_cast(particleSet.pressure.data()),
		particleSet.particleMass,
		particleSet.maxDensityErr,
		particleSet.count
		);

	// Initialize kernel
	float kernelRadius = KERNEL_SCALE * scene->particleRadius;
	kernel.init(kernelRadius);
	std::cout << "Kernel.H = " << kernel.H << std::endl;

	// Initialize grid
	this->fluidGrid.init(scene->boxSize, kernel.H);
	this->boundaryGrid.init(scene->boxSize, kernel.H);
	updateFluidGrid();


	// Initialize density variance scale
	this->initDensityVarianceScale();
	std::cout << "densityVarianceScale = " << this->densityVarianceScale << std::endl;

	offset.resize(fluidGrid.cellNumber + 1);


	d_fluidGrid = thrust::device_new<Grid>();
	initDeviceGrid << <1, 1 >> > (thrust::raw_pointer_cast(d_fluidGrid), fluidGrid.cellSize, fluidGrid.boxSize, fluidGrid.gridSize, fluidGrid.cellNumber);
	d_kernel = thrust::device_new<Kernel>();
	initDeviceKernel << <1, 1 >> > (thrust::raw_pointer_cast(d_kernel), kernelRadius);

	this->relax();

	std::cout << std::endl;
}
void Simulator::updateFluidGrid() {

	particleSet.h_velocity = particleSet.velocity;
	particleSet.h_position = particleSet.position;

	h_offset = offset;

	this->fluidGrid.update(particleSet.h_position, h_offset, [this](size_t i, size_t j) {
		std::swap(particleSet.h_position[i], particleSet.h_position[j]);
		std::swap(particleSet.h_velocity[i], particleSet.h_velocity[j]);
	});

	offset = h_offset;

	particleSet.position = particleSet.h_position;
	particleSet.velocity = particleSet.h_velocity;
}

__global__
void computeDensity(Grid* grid, DeviceParticleSetPointer* particles, size_t* offset, Kernel* kernel) {

	// compute fluid particles' densities
	for (size_t i = 0; i < particles->count; i++) {
		float fluidTerm = 0.0f;
		PCISPH::Vec3 pos = particles->position[i];

		glm::u64vec3 boundBoxMin = grid->getGridPos(pos - PCISPH::Vec3(grid->cellSize));
		glm::u64vec3 boundBoxMax = grid->getGridPos(pos + PCISPH::Vec3(grid->cellSize));

		for (auto z = boundBoxMin.z; z <= boundBoxMax.z; z++) {
			for (auto y = boundBoxMin.y; y <= boundBoxMax.y; y++) {
				for (auto x = boundBoxMin.x; x <= boundBoxMax.x; x++) {
					size_t cellIndex = grid->linearIndex(PCISPH::uVec3(x, y, z));
					for (size_t neighborIndex = offset[cellIndex]; neighborIndex < offset[cellIndex + 1]; neighborIndex++) {
						PCISPH::Vec3 r = pos - particles->position[neighborIndex];
						float rLen = PCISPH::length(r);
						if (rLen > kernel->H || rLen < EPSILON) return;
						fluidTerm += kernel->poly6Kernel(r);
					}
				}
			}
		}

		particles->density[i] = particles->particleMass * fluidTerm;
	}
}

__global__
void computeNormal(Grid* grid, DeviceParticleSetPointer* particles, size_t* offset, Kernel* kernel) {
	for (auto i = 0; i < particles->count; i++) {
		PCISPH::Vec3 n(0.0f);

		glm::u64vec3 boundBoxMin = grid->getGridPos(particles->position[i] - PCISPH::Vec3(grid->cellSize));
		glm::u64vec3 boundBoxMax = grid->getGridPos(particles->position[i] + PCISPH::Vec3(grid->cellSize));

		for (auto z = boundBoxMin.z; z <= boundBoxMax.z; z++) {
			for (auto y = boundBoxMin.y; y <= boundBoxMax.y; y++) {
				for (auto x = boundBoxMin.x; x <= boundBoxMax.x; x++) {
					size_t cellIndex = grid->linearIndex(PCISPH::uVec3(x, y, z));
					for (size_t neighborIndex = offset[cellIndex]; neighborIndex < offset[cellIndex + 1]; neighborIndex++) {
						PCISPH::Vec3 r = particles->position[i] - particles->position[neighborIndex];
						float rLen = PCISPH::length(r);
						if (rLen > kernel->H || rLen < EPSILON) return;
						n += kernel->poly6KernelGradient(r) / particles->density[neighborIndex];
					}
				}
			}
		}

		n *= kernel->H * particles->particleMass;
		particles->normal[i] = n;
	}
}

__global__
void computeForces(Grid* grid, DeviceParticleSetPointer* particleSet, size_t* offset, Kernel* kernel,
	float referenceDensity, float viscosityCoefficient, float surfaceTensionCoefficient, PCISPH::Vec3 gravity) {
	float squaredMass = particleSet->particleMass * particleSet->particleMass;

	for (auto i = 0; i < particleSet->count; i++) {
		PCISPH::Vec3 viscosity(0.f);
		PCISPH::Vec3 cohesion(0.f);
		PCISPH::Vec3 curvature(0.f);

		glm::u64vec3 boundBoxMin = grid->getGridPos(particleSet->position[i] - PCISPH::Vec3(grid->cellSize));
		glm::u64vec3 boundBoxMax = grid->getGridPos(particleSet->position[i] + PCISPH::Vec3(grid->cellSize));

		for (auto z = boundBoxMin.z; z <= boundBoxMax.z; z++) {
			for (auto y = boundBoxMin.y; y <= boundBoxMax.y; y++) {
				for (auto x = boundBoxMin.x; x <= boundBoxMax.x; x++) {
					size_t cellIndex = grid->linearIndex(PCISPH::uVec3(x, y, z));
					for (size_t j = offset[cellIndex]; j < offset[cellIndex + 1]; j++) {
						PCISPH::Vec3 r = particleSet->position[i] - particleSet->position[j];
						float rLen = PCISPH::length(r);
						if (rLen > kernel->H || rLen < EPSILON) return;

						PCISPH::Vec3 vDiff = particleSet->velocity[i] - particleSet->velocity[j];
						viscosity -= vDiff * kernel->viscosityKernelLaplacian(r) / particleSet->density[j];

						float Kij = 2.0f * referenceDensity / (particleSet->density[i] + particleSet->density[j]);
						cohesion += Kij * (r / rLen) * kernel->cohesionKernel(rLen);
						curvature += Kij * (particleSet->normal[i] - particleSet->normal[j]);
					}
				}
			}
		}



		viscosity *= viscosityCoefficient * squaredMass / particleSet->density[i];
		cohesion *= -surfaceTensionCoefficient * squaredMass;
		curvature *= -surfaceTensionCoefficient * particleSet->particleMass;

		particleSet->forces[i] = viscosity + cohesion + curvature + particleSet->particleMass * gravity;
	}
}

void Simulator::clearPressureAndPressureForce() {
	thrust::fill(particleSet.pressure.begin(), particleSet.pressure.end(), 0);
	thrust::fill(particleSet.pressureForce.begin(), particleSet.pressureForce.end(), PCISPH::Vec3(0, 0, 0));
}

__global__
void predictVelocityAndPosition(DeviceParticleSetPointer* particleSet, float timeStep) {
	for (int i = 0; i < particleSet->count; i++) {
		PCISPH::Vec3 acceleration = (particleSet->forces[i] + particleSet->pressureForce[i]) / particleSet->particleMass;
		particleSet->predictVelocity[i] = particleSet->velocity[i] + timeStep * acceleration;
		particleSet->predictPosition[i] = particleSet->position[i] + particleSet->predictVelocity[i] * timeStep;
	}
}

__global__
void updatePressure(Grid* grid, DeviceParticleSetPointer* particleSet, size_t* offset, Kernel* kernel, float referenceDensity, float densityVarianceScale) {

	for (int i = 0; i < particleSet->count; i++) {
		float fDensity = 0.f;

		glm::u64vec3 boundBoxMin = grid->getGridPos(particleSet->position[i] - PCISPH::Vec3(grid->cellSize));
		glm::u64vec3 boundBoxMax = grid->getGridPos(particleSet->position[i] + PCISPH::Vec3(grid->cellSize));

		for (auto z = boundBoxMin.z; z <= boundBoxMax.z; z++) {
			for (auto y = boundBoxMin.y; y <= boundBoxMax.y; y++) {
				for (auto x = boundBoxMin.x; x <= boundBoxMax.x; x++) {
					size_t cellIndex = grid->linearIndex(PCISPH::uVec3(x, y, z));
					for (size_t j = offset[cellIndex]; j < offset[cellIndex + 1]; j++) {
						PCISPH::Vec3 r = particleSet->predictPosition[i] - particleSet->predictPosition[j];
						float rLen = PCISPH::length(r);
						if (rLen > kernel->H || rLen < EPSILON) return;
						fDensity += kernel->poly6Kernel(r);
					}
				}
			}
		}

		float bDensity = 0.f;

		float density = fDensity * particleSet->particleMass + bDensity;
		float densityVariation = max(0.f, density - referenceDensity);
		particleSet->maxDensityErr = max(particleSet->maxDensityErr, densityVariation);

		particleSet->pressure[i] += densityVarianceScale * densityVariation;
	}
}

__global__
void updatePressureForce(Grid* grid, DeviceParticleSetPointer* particleSet, size_t* offset, Kernel* kernel) {
	for (int i = 0; i < particleSet->count; i++) {
		PCISPH::Vec3 pressureForce(0.f);

		glm::u64vec3 boundBoxMin = grid->getGridPos(particleSet->position[i] - PCISPH::Vec3(grid->cellSize));
		glm::u64vec3 boundBoxMax = grid->getGridPos(particleSet->position[i] + PCISPH::Vec3(grid->cellSize));

		for (auto z = boundBoxMin.z; z <= boundBoxMax.z; z++) {
			for (auto y = boundBoxMin.y; y <= boundBoxMax.y; y++) {
				for (auto x = boundBoxMin.x; x <= boundBoxMax.x; x++) {
					size_t cellIndex = grid->linearIndex(PCISPH::uVec3(x, y, z));
					for (size_t j = offset[cellIndex]; j < offset[cellIndex + 1]; j++) {
						PCISPH::Vec3 r = particleSet->predictPosition[i] - particleSet->predictPosition[j];
						float rLen = PCISPH::length(r);
						if (rLen > kernel->H || rLen < 1e-5) return;
						float term1 = particleSet->pressure[i] / particleSet->density[i] / particleSet->density[i];
						float term2 = particleSet->pressure[j] / particleSet->density[j] / particleSet->density[j];
						pressureForce -= particleSet->particleMass * (term1 + term2) * kernel->spikyKernelGradient(r);
					}
				}
			}
		}

		pressureForce *= particleSet->particleMass;
		particleSet->pressureForce[i] = pressureForce;
	}
}

__global__
void updateVelocityAndPosition(DeviceParticleSetPointer* particleSet, float timeStep) {
	for (int i = 0; i < particleSet->count; i++) {
		PCISPH::Vec3 acceleration = (particleSet->forces[i] + particleSet->pressureForce[i]) / particleSet->particleMass;
		particleSet->velocity[i] += acceleration * timeStep;
		particleSet->position[i] += particleSet->velocity[i] * timeStep;
	}
}

void Simulator::updateScene() {
	// TODO
}

void Simulator::relax() {
	update(10000);

	thrust::fill(particleSet.velocity.begin(), particleSet.velocity.end(), scene->initVelocity);
}

__device__
void collisionFunc(DeviceParticleSetPointer* particleSet, float restitution, size_t i, const PCISPH::Vec3 &n, const float d)
{
	particleSet->position[i] += d * n;
	particleSet->velocity[i] -= (1 + restitution) * PCISPH::dot(particleSet->velocity[i], n) * n;
}


__global__
void handleCollision(DeviceParticleSetPointer* particleSet, float restitution, PCISPH::Vec3 box) {

	for (int i = 0; i < particleSet->count; i++) {
		const auto &p = particleSet->position[i];
		PCISPH::Vec3 n;
		if (p.x < 0) {
			collisionFunc(particleSet, restitution, i, PCISPH::Vec3(1.f, 0.f, 0.f), -p.x);
		}
		if (p.x > box.x) {
			collisionFunc(particleSet, restitution, i, PCISPH::Vec3(-1.f, 0.f, 0.f), p.x - box.x);
		}
		if (p.y < 0) {
			collisionFunc(particleSet, restitution, i, PCISPH::Vec3(0.f, 1.f, 0.f), -p.y);
		}
		if (p.y > box.y) {
			collisionFunc(particleSet, restitution, i, PCISPH::Vec3(0.f, -1.f, 0.f), p.y - box.y);
		}
		if (p.z < 0) {
			collisionFunc(particleSet, restitution, i, PCISPH::Vec3(0.f, 0.f, 1.f), -p.z);
		}
		if (p.z > box.z) {
			collisionFunc(particleSet, restitution, i, PCISPH::Vec3(0.f, 0.f, -1.f), p.z - box.z);
		}
	}
}

void Simulator::initParticlesPositions() {
	float r = scene->particleRadius;
	float d = 2 * r;

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

void Simulator::update(const size_t maxIterations) {


	updateFluidGrid();
	computeDensity<<<1, 1>>>(thrust::raw_pointer_cast(d_fluidGrid), thrust::raw_pointer_cast(particleSetPointer), thrust::raw_pointer_cast(offset.data()), thrust::raw_pointer_cast(d_kernel));
	computeNormal<<<1,1>>>(thrust::raw_pointer_cast(d_fluidGrid), thrust::raw_pointer_cast(particleSetPointer), thrust::raw_pointer_cast(offset.data()), thrust::raw_pointer_cast(d_kernel));
	computeForces<<<1,1>>>(thrust::raw_pointer_cast(d_fluidGrid), thrust::raw_pointer_cast(particleSetPointer), thrust::raw_pointer_cast(offset.data()), thrust::raw_pointer_cast(d_kernel), scene->referenceDensity, scene->viscosityCoefficient, scene->surfaceTensionCoefficient, scene->gravity);
	clearPressureAndPressureForce();
	size_t iterations = 0;
	while (iterations < maxIterations) {
		particleSet.maxDensityErr = 0.f;
		predictVelocityAndPosition<<<1,1>>>(thrust::raw_pointer_cast(particleSetPointer), scene->timeStep);
		updatePressure<<<1,1>>>(thrust::raw_pointer_cast(d_fluidGrid), thrust::raw_pointer_cast(particleSetPointer), thrust::raw_pointer_cast(offset.data()), thrust::raw_pointer_cast(d_kernel), scene->referenceDensity, densityVarianceScale);
		updatePressureForce<<<1,1>>>(thrust::raw_pointer_cast(d_fluidGrid), thrust::raw_pointer_cast(particleSetPointer), thrust::raw_pointer_cast(offset.data()), thrust::raw_pointer_cast(d_kernel));
		static const int MIN_ITERATIONS = 3;
		static const float yita = 0.01;
		static const float TOL = yita * this->scene->referenceDensity;
		if (++iterations >= MIN_ITERATIONS && this->particleSet.maxDensityErr < TOL) {
			break;
		}
	}

	updateVelocityAndPosition<<<1,1>>>(thrust::raw_pointer_cast(particleSetPointer), scene->timeStep);
	handleCollision<<<1,1>>>(thrust::raw_pointer_cast(particleSetPointer), scene->restitution, scene->boxSize);

	static size_t count = 0;
	std::cout << "maxDensityErr = " << this->particleSet.maxDensityErr << std::endl;
	std::cout << "iterations = " << iterations << std::endl;
	std::cout << "count = " << count << std::endl;
	std::cout << std::endl;
	count++;
}

