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

	//long particleNumber;
	//float particleRadius;  // It seems nothing to do with particle's radius
	float referenceDensity;
	//float surfaceTension;  //TODO: What is surface tension?
	float viscosityCoefficient;
	float timeStep;
	float particleRadius;
	PCISPH::Vec3 gravity;
	PCISPH::Vec3 boxSize;
	PCISPH::Vec3 fluidSize;
	PCISPH::Vec3 fluidPosition;

private:
	void parseFromFile(const std::string &path);
};

