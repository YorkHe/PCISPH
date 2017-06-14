#pragma once

#include "Renderer.h"
#include "Simulator.h"
#include "Scene.h"

class System
{
public:
	System();
	~System();

	void init(const std::string &configPath);
	void start();
	void terminate();

private:
	Scene *scene;
	Simulator *simulator;
	Window *window;
	Renderer *renderer;
};

