#pragma once
#include "Window.h"
#include "Engine.h"

class System
{
public:
	System();
	~System();

	void init();
	void start();
	void terminate();

private:

	Window* mWindow;
	Engine* mEngine;


};

