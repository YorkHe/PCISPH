#include "Engine.h"



Engine::Engine()
{
}


Engine::~Engine()
{
}

void Engine::init(Scene* scene)
{
	mScene = scene;
	mParticles = new Particle();

	initShaders();
	loadCamera();
	initParticles();
}

void Engine::initParticles()
{
	PCISPH::Vec3 leftBottom = mScene->scene.boxBounds[0];
	PCISPH::Vec3 rightTop = mScene->scene.boxBounds[1];

	float step = mScene->constant.particleRadius * 1.5f;

	for (float x = leftBottom.x; x < rightTop.x; x+=step)
	{
		for (float y = leftBottom.y; y < rightTop.y; y+=step)
		{
			for (float z = leftBottom.z; z< rightTop.z; z+=step)
			{
				mParticles->appendParticle(x, y, z);
			}
		}
	}
}


void Engine::initShaders()
{
}

void Engine::loadCamera()
{
	mCamera = new Camera;
	mCamera->Position = mScene->scene.cameraPosition;
	mCamera->Up = PCISPH::Vec3(0.f, 1.0f, 0.f);
	mCamera->Front = mScene->scene.cameraTarget - mScene->scene.cameraPosition;

	mCamera->Yaw = 0;
	mCamera->Pitch = 0;
}

void Engine::draw()
{
	render();

	updateSimulation();
}

void Engine::updateSimulation()
{

}

void Engine::render()
{
	mParticles->draw(mCamera);
}



