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

	long particleNumber;
	//float particleRadius;  // It seems nothing to do with particle's radius
	float referenceDensity;
	//float surfaceTension;  //TODO: What is surface tension?
	float viscosityCoefficient;
	float timeStep;
	PCISPH::Vec3 gravity;

	std::vector<PCISPH::Vec3> boxBounds;
	std::vector<PCISPH::Vec3> fluidBounds;

private:
	void parseFromFile(const std::string &path);
};

