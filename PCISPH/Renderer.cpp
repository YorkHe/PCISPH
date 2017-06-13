#include "Renderer.h"

#include <glm/gtc/type_ptr.hpp>

Renderer::Renderer():window(nullptr), scene(nullptr), camera()
{
}


Renderer::~Renderer()
{

}

void Renderer::init(Window *window, const Scene *scene) {
	this->window = window;
	this->scene = scene;
	this->particleShader = Shader("shaders/particle.vert", "shaders/particle.frag");
	this->sceneShader = Shader("shaders/scene.vert", "shaders/scene.frag");

	// codes below are just for test
	glGenVertexArrays(1, &(this->sceneVAO));
	glGenVertexArrays(1, &(this->particleVAO));
	this->initSceneVAO();
}

void Renderer::initSceneVAO() {
	PCISPH::Vec3 worldP1(0, 0, 0);
	PCISPH::Vec3 worldP2 = this->scene->boxSize;

	GLfloat vertices[]{
		worldP1.x, worldP1.y, worldP1.z,
		worldP1.x, worldP1.y, worldP2.z,
		worldP2.x, worldP1.y, worldP2.z,
		worldP2.x, worldP1.y, worldP1.z,
		worldP1.x, worldP2.y, worldP1.z,
		worldP1.x, worldP2.y, worldP2.z,
		worldP2.x, worldP2.y, worldP2.z,
		worldP2.x, worldP2.y, worldP1.z
	};

	GLuint indices[] = {
		0, 1,
		1, 2,
		2, 3,
		3, 0,
		0, 4,
		1, 5,
		2, 6,
		3, 7,
		4, 5,
		5, 6,
		6, 7,
		7, 4
	};

	GLuint VBO, EBO;
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(sceneVAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
}

void Renderer::initParticleVAO(const std::vector<PCISPH::Vec3> &points) {
	GLuint VBO;
	glGenBuffers(1, &VBO);

	glBindVertexArray(this->particleVAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(points[0]) * points.size(), &points[0], GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(points[0]), (void*)0);
	glEnableVertexAttribArray(0);
	glBindVertexArray(0);

	glDeleteBuffers(1, &VBO);
}

void Renderer::draw(const ParticleSet &particleSet) {
	// TODO: Draw the scene and particles

	// The codes below are just for test
	glm::mat4 model;
	model = glm::translate(model, glm::vec3(-0.25f, -0.25f, -0.25f));
	model = glm::scale(model, glm::vec3(1.0f, 1.0f, 1.0f));

	glm::mat4 mvp = camera.getProjViewMatrix() * model;

	sceneShader.use();
	glUniformMatrix4fv(glGetUniformLocation(sceneShader.program, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
	glBindVertexArray(this->sceneVAO);
	glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, (GLvoid*)0);
	glBindVertexArray(0);
	sceneShader.unUse();

	initParticleVAO(particleSet.position);
	particleShader.use();
	glPointSize(2);
	glUniformMatrix4fv(glGetUniformLocation(particleShader.program, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
	glBindVertexArray(this->particleVAO);
	glDrawArrays(GL_POINTS, 0, (GLsizei)particleSet.count);
	glBindVertexArray(0);
	particleShader.unUse();
}