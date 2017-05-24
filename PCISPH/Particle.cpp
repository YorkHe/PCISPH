#include "Particle.h"

#include <glm/gtc/type_ptr.hpp>

Particle::Particle()
{

	radius = 2.0;
	mShader = Shader("shaders/particle.vert", "shaders/particle.frag");
}

Particle::Particle(float r)
{
	Particle();
	radius = r;
}


Particle::~Particle()
{

}

void Particle::initBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(position), &position[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);

	glBindVertexArray(0);
}

void Particle::appendParticle(float x, float y, float z)
{
	vec3 p = { x, y, z };
	vec3 zero = { 0, 0, 0 };
	position.push_back(p);

	velocity.push_back(zero);
	density.push_back(zero);
	pressure.push_back(zero);
	extForce.push_back(zero);
	pressureForce.push_back(zero);
}


void Particle::draw(Camera* camera)
{
	mShader.use();
	glUniformMatrix4fv(glGetUniformLocation(mShader.program, "MVP"), 1, GL_FALSE, glm::value_ptr(camera->getViewMatrix()));
	glBindVertexArray(VAO);
	glDrawArrays(GL_POINTS, 0, position.size());
	glBindVertexArray(0);
	mShader.unUse();
}
