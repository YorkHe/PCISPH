#pragma once

#include "Types.h"
#include <vector>
#include <functional>

class Grid
{
public:
	Grid();
	~Grid();

	void init(const PCISPH::Vec3 boxSize, float cellSize);
	void update(const std::vector<PCISPH::Vec3> &positions, std::function<void(size_t, size_t)> swap);
	
	//template<typename Func>
	void query(const PCISPH::Vec3 &pos, std::function<void(size_t)> func) const;

private:
	float cellSize;
	PCISPH::Vec3 boxSize;
	PCISPH::uVec3 gridSize;
	size_t cellNumber;
	std::vector<size_t> offset;

	PCISPH::uVec3 getGridPos(const PCISPH::Vec3 &pos) const;
	size_t linearIndex(const PCISPH::Vec3 &pos) const;
	size_t linearIndex(const PCISPH::uVec3 &gridPos) const;
};

