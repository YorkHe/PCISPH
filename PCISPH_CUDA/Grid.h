#pragma once

#include "Types.h"
#include <vector>
#include <functional>

#include <thrust/device_vector.h>

class Grid
{
public:
	__host__ __device__
	Grid();
	__host__ __device__
	~Grid();

	__host__ __device__
	void init(const PCISPH::Vec3 boxSize, float cellSize);
	void update(const thrust::host_vector<PCISPH::Vec3> &positions, const thrust::host_vector<size_t> offset, std::function<void(size_t, size_t)> swap);
	
	//template<typename Func>
	float cellSize;
	PCISPH::Vec3 boxSize;
	PCISPH::uVec3 gridSize;
	size_t cellNumber;

	__host__ __device__
	PCISPH::uVec3 getGridPos(const PCISPH::Vec3 &pos) const;


	__host__ __device__
	size_t linearIndex(const PCISPH::Vec3 &pos) const;
	__host__ __device__
	size_t linearIndex(const PCISPH::uVec3 &gridPos) const;
};
