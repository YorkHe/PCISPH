#pragma once

#include "Types.h"

class Kernel {
public:
	__host__ __device__
	Kernel():h(0), H(this->h) {}
	__host__ __device__
	~Kernel() {}

	__host__ __device__
	void operator = (const Kernel k) {
		this->H = k.H;
		this->init(k.H);
	}

	__host__ __device__
	void init(const float h) {
		this->H = h;
		this->h = h;
		this->h2 = h * h;
		this->h6 = this->h2 * this->h2 * this->h2;
		this->h9 = h * this->h2 * this->h6;

		const float PI = 3.1415926f;

		this->poly6Factor = 315.f / (64.f * PI * this->h9);
		this->poly6GradFactor = -945.f / (32.f * PI * this->h9);
		this->poly6LapFactor = -945.f / (32.f * PI * this->h9);

		this->spikyGradFactor = -45.f / (PI * this->h6);

		this->viscosityLapFactor = 45.f / (PI * this->h6);

		this->cohesionFactor = 32.0f / (PI * this->h9);
	}

	/*
	* Smooth kernels
	* h means smooth length. if h = 0 as default, use the previously setted value.
	*/

	// Poly6 kernel is used to interpolate density
	__device__
	float poly6Kernel(const PCISPH::Vec3 &r) {
		float rNorm = this->normalize(r);
		if (rNorm <= this->h) {
			float val = this->h2 - rNorm * rNorm;
			return this->poly6Factor * val * val * val;
		}
		else {
			return 0;
		}
	}
	__device__
	PCISPH::Vec3 poly6KernelGradient(const PCISPH::Vec3 &r) {
		float rNorm = this->normalize(r);
		if (rNorm <= this->h) {
			float val = this->h2 - rNorm * rNorm;
			return this->poly6GradFactor * val * val * r;
		}
		else {
			return PCISPH::Vec3(0, 0, 0);
		}
	}

	PCISPH::Vec3 h_poly6KernelGradient(const PCISPH::Vec3 &r) {
		float rNorm = this->slow_normalize(r);
		if (rNorm <= this->h) {
			float val = this->h2 - rNorm * rNorm;
			return this->poly6GradFactor * val * val * r;
		}
		else {
			return PCISPH::Vec3(0, 0, 0);
		}
	}



	__device__
	float poly6KernelLaplacian(const PCISPH::Vec3 &r) {
		float rNorm = this->normalize(r);
		if (rNorm <= this->h) {
			float r2 = rNorm * rNorm;
			return this->poly6LapFactor * (this->h2 - r2) * (3 * this->h2 - 7 * r2);
		}
		else {
			return 0;
		}
	}

	// Spiky kernel is used to interpolate pressure
	//inline float spikyKernel(const PCISPH::Vec3 &r);
	__device__
	PCISPH::Vec3 spikyKernelGradient(const PCISPH::Vec3 &r) {
		float rNorm = this->normalize(r);
		if (rNorm <= this->h) {
			float val = this->h - rNorm;
			return this->spikyGradFactor * val * val / rNorm * r;
		}
		else {
			return PCISPH::Vec3(0, 0, 0);
		}
	}

	PCISPH::Vec3 h_spikyKernelGradient(const PCISPH::Vec3 &r) {
		float rNorm = this->slow_normalize(r);
		if (rNorm <= this->h) {
			float val = this->h - rNorm;
			return this->spikyGradFactor * val * val / rNorm * r;
		}
		else {
			return PCISPH::Vec3(0, 0, 0);
		}
	}

	// viscosity kernel
	//inline float viscosityKernel(const PCISPH::Vec3 &r);
	__device__
	float viscosityKernelLaplacian(const PCISPH::Vec3 &r) {
		float rNorm = this->normalize(r);
		if (rNorm <= this->h) {
			return this->viscosityLapFactor * (this->h - rNorm);
		}
		else {
			return 0;
		}
	}

	// cohesion kernel
	__device__
	float cohesionKernel(const float r) {
		if (2 * r > this->h && r <= h) {
			float t = this->h - r;
			return this->cohesionFactor * t * t * t * r * r * r;
		}
		else if (r > 0 && 2 * r <= h) {
			float t = this->h - r;
			return this->cohesionFactor * (2 * t * t * t * r * r * r - this->h6 / 64.0f);
		}
		else {
			return 0;
		}
	}

private:
	float h;
public:
	float H;

private:
	float h2;
	float h6;
	float h9;

	float poly6Factor;
	float poly6GradFactor;
	float poly6LapFactor;

	//float spikyFactor;
	float spikyGradFactor;

	//float viscosityFactor;
	float viscosityLapFactor;

	float cohesionFactor;

	__host__ __device__
	void initFactors();

	float slow_normalize(const PCISPH::Vec3 &r) {
		return sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
	}

	__device__
	float normalize(const PCISPH::Vec3 &r) {
		return sqrtf(r.x * r.x + r.y * r.y + r.z * r.z);
	}
};

