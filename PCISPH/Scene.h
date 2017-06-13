#pragma once

#include <string>
#include <vector>
#include "Types.h"

class Scene
{
public:
	Scene(const std::string &path);
	~Scene();

	enum {STATIC, FLOW} mode;

	float referenceDensity;
	float particleRadius;
	float viscosityCoefficient;
	float surfaceTensionCoefficient;
	float gasStiffness;
	size_t kernelParticles;

	float timeStep;

	float restitution;
	PCISPH::Vec3 gravity;
	PCISPH::Vec3 boxSize;
	PCISPH::Vec3 fluidSize;
	PCISPH::Vec3 fluidPosition;
	PCISPH::Vec3 initVelocity;

private:
	void parseFromFile(const std::string &path);
};

