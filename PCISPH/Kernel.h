#pragma once

#include "Types.h"

class Kernel {
public:
	Kernel();
	Kernel(float h);
	~Kernel();

	void setSmoothLength(const float h);

	/*
	* Smooth kernels
	* h means smooth length. if h = 0 as default, use the previously setted value.
	*/

	// Poly6 kernel is used to interpolate density
	inline float poly6Kernel(const PCISPH::Vec3 &r);
	inline PCISPH::Vec3 poly6KernelGradient(const PCISPH::Vec3 &r);
	inline float poly6KernelLaplacian(const PCISPH::Vec3 &r);

	// Spiky kernel is used to interpolate pressure
	//inline float spikyKernel(const PCISPH::Vec3 &r);
	inline PCISPH::Vec3 spikyKernelGradient(const PCISPH::Vec3 &r);

	// viscosity kernel
	//inline float viscosityKernel(const PCISPH::Vec3 &r);
	inline float viscosityKernelLaplacian(const PCISPH::Vec3 &r);

private:
	float h;

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

	void initFactors();
	inline float normalize(const PCISPH::Vec3 &r);
};

