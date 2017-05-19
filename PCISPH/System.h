#pragma once
#include "Scene.h"
#include "Window.h"
#include "Engine.h"
#include <boost/program_options.hpp>

class System
{
public:
	System();
	~System();

	void init(boost::program_options::variables_map vm);
	void start();
	void terminate();

private:

	Scene* mScene;

	Window* mWindow;
	Engine* mEngine;


};

