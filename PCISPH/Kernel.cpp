#include "Kernel.h"
#include <cmath>

const float PI = 3.1415926;

Kernel::Kernel() {
	const float defaultH = 3.0f; // TODO: 3.0f is just a test value
	setSmoothLength(defaultH);
	initFactors();
}

Kernel::Kernel(float h) {
	setSmoothLength(h);
	initFactors();
}

Kernel::~Kernel() {}

void Kernel::setSmoothLength(const float h) {
	this->h = h;
	this->h2 = h * h;
	this->h6 = this->h2 * this->h2 * this->h2;
	this->h9 = h * this->h2 * this->h6;
}

void Kernel::initFactors() {
	this->poly6Factor = 315.f / (64.f * PI * this->h9);
	this->poly6GradFactor = -945.f / (32.f * PI * this->h9);
	this->poly6LapFactor = 945.f / (8.f * PI * this->h9);

	this->spikyGradFactor = -45.f / (PI * this->h6);

	this->viscosityLapFactor = 45.f / (PI * this->h6);
}

inline float Kernel::poly6Kernel(const PCISPH::Vec3 &r) {
	float rNorm = this->normalize(r);
	if (rNorm <= this->h) {
		float val = this->h2 - rNorm * rNorm;
		return this->poly6Factor * val * val * val;
	}
	else {
		return 0;
	}
}

inline PCISPH::Vec3 Kernel::poly6KernelGradient(const PCISPH::Vec3 &r) {
	float rNorm = this->normalize(r);
	if (rNorm <= this->h) {
		float val = this->h2 - rNorm * rNorm;
		return this->poly6GradFactor * val * val * r;
	}
	else {
		return PCISPH::Vec3(0, 0, 0);
	}
}

inline float Kernel::poly6KernelLaplacian(const PCISPH::Vec3 &r) {
	float rNorm = this->normalize(r);
	if (rNorm <= this->h) {
		float val = this->h2 - rNorm * rNorm;
		return this->poly6LapFactor * val * val * (rNorm * rNorm - 0.75 * val);
	}
	else {
		return 0;
	}
}

inline PCISPH::Vec3 Kernel::spikyKernelGradient(const PCISPH::Vec3 &r) {
	float rNorm = this->normalize(r);
	if (rNorm <= this->h) {
		float val = this->h - rNorm;
		return this->spikyGradFactor * val * val / rNorm * r;
	}
	else {
		return PCISPH::Vec3(0, 0, 0);
	}
}

inline float Kernel::viscosityKernelLaplacian(const PCISPH::Vec3 &r) {
	float rNorm = this->normalize(r);
	if (rNorm <= this->h) {
		return this->viscosityLapFactor * (this->h - rNorm);
	}
	else {
		return 0;
	}
}

inline float Kernel::normalize(const PCISPH::Vec3 &r) {
	return sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
}