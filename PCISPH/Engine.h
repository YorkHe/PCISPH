#pragma once
#include "Scene.h"
#include "Shader.h"
#include "Camera.h"
#include "Particle.h"

class Engine
{
public:
	Engine();
	~Engine();

	void init(Scene*);
	void draw();
	void loadCamera();
	void initShaders();
	void initParticles();

	void render();
	void updateSimulation();



private:
	Scene* mScene;
	Camera* mCamera;

	Particle* mParticles;
};

