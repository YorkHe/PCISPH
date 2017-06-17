#pragma once
#define GLM_FORCE_CUDA
#include <cuda.h>
#include <glm/glm.hpp>
#include <glm/detail/precision.hpp>

namespace PCISPH
{
	typedef glm::vec2 Vec2;
	typedef glm::vec3 Vec3;
	typedef glm::ivec3 iVec3;

	typedef struct _tri {
		int a;
		int b;
		int c;
	}Triangle;
	typedef struct _vertex {
		Vec3 position;
		double value;
	} Vertex;

	typedef struct _grid {
		Vertex* v[8];
	} Grid;

	typedef glm::u64vec3 uVec3;

	typedef glm::u64 uint;

}
