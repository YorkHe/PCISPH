#pragma once
#include "Scene.h"

class Engine
{
public:
	Engine();
	~Engine();

	void init(Scene*);

private:
	Scene* mScene;

};

