#pragma once
#include <vector>
#include <gl/glew.h>
#include "Shader.h"
#include "Camera.h"

class Particle
{
public:
	Particle();
	Particle(float r);
	~Particle();
	void draw(Camera* camera);
	void initBuffer();

	void appendParticle(float x, float y, float z);

	float radius;

	typedef struct _vec3
	{
		float x;
		float y;
		float z;
	}vec3;

	std::vector<vec3> position;
	std::vector<vec3> velocity;
	std::vector<vec3> density;
	std::vector<vec3> pressure;
	std::vector<vec3> extForce;
	std::vector<vec3> pressureForce;

	GLuint VAO;
	GLuint VBO;

private:
	Shader mShader;
};

