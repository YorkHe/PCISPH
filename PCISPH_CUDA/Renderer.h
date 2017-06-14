#pragma once

#include "Window.h"
#include "Shader.h"
#include "Scene.h"
#include "Types.h"
#include "MarchingCube.h"
#include "Camera.h"
#include "ParticleSet.h"

class Renderer
{
public:
	Renderer();
	~Renderer();

	void init(Window *window, const Scene *scene);

	void draw(const ParticleSet &particleSet);

	double lastTime;
	int frames = 0;
private:

	const float gridSize = 0.0002;

	const ParticleSet* mParticleSet;
	Window *window;
	const Scene *scene;
	Camera camera;
	// TODO: Shaders
	Shader sceneShader;
	Shader particleShader;
	MarchingCube* marchingCube;

	glm::mat4 mvp;

	// codes below are just for test
	GLuint sceneVAO;
	GLuint particleVAO;
	GLuint meshVAO;
	void initSceneVAO();
	void updateParticleVAO(const std::vector<PCISPH::Vec3> &points);
	void updateMeshVAO();

	void drawScene();

	void drawParticle();
	void drawMesh();
	void drawFps();
};

