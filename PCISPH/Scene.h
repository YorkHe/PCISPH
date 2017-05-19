#pragma once
#include <string>
#include "Types.h"
#include <vector>

class Scene
{
public:
	Scene();
	~Scene();
	void parseFromFile(std::string path);

	struct _constant
	{
		float particleRadius = 0.0;
		float restDensity = 0.0;
		float surfaceTension = 0.0;
		float viscosity = 0.0;
		float timeStep = 0.0;
		PCISPH::Vec3 gravity;
	};

	struct _scene
	{
		PCISPH::Vec3 cameraPosition;
		PCISPH::Vec3 cameraTarget;
		std::vector<PCISPH::Vec3> worldBounds;
		std::vector<PCISPH::Vec3> boxBounds;
	};


	_constant constant;

	_scene scene;

};



