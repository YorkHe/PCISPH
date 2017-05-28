#include "Renderer.h"



Renderer::Renderer():window(nullptr), scene(nullptr)
{
}


Renderer::~Renderer()
{

}

void Renderer::init(Window *window, const Scene *scene) {
	this->window = window;
	this->scene = scene;
	this->particleShader = Shader("shaders/particle.vert", "shader/particle.frag");
	this->sceneShader = Shader("shaders/scene.vert", "shader/scene.frag");
}

void Renderer::draw(const std::vector<PCISPH::Vec3> &points) {
	// TODO: Draw the scene and particles
}