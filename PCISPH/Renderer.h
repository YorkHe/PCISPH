#pragma once

#include "Window.h"
#include "Shader.h"
#include "Scene.h"
#include "Types.h"
#include "Camera.h"
#include "ParticleSet.h"

class Renderer
{
public:
	Renderer();
	~Renderer();

	void init(Window *window, const Scene *scene);

	void draw(const ParticleSet &particleSet);

private:
	Window *window;
	const Scene *scene;
	Camera camera;
	// TODO: Shaders
	Shader sceneShader;
	Shader particleShader;

	// codes below are just for test
	GLuint sceneVAO;
	GLuint particleVAO;
	void initSceneVAO();
	void initParticleVAO(const std::vector<PCISPH::Vec3> &points);
	void drawIndividualParticle(const std::vector<PCISPH::Vec3> &points);
	void drawMesh(const std::vector<PCISPH::Vec3> &points);
};

