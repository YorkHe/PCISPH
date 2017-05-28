#pragma once

#include "Window.h"
#include "Shader.h"
#include "Scene.h"
#include "Types.h"

class Renderer
{
public:
	Renderer();
	~Renderer();

	void init(Window *window, const Scene *scene);

	void draw(const std::vector<PCISPH::Vec3> &points);

private:
	Window *window;
	const Scene *scene;

	// TODO: Shaders
	Shader sceneShader;
	Shader particleShader;
};

