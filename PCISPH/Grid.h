#pragma once

#include "Types.h"
#include <vector>

class Grid
{
public:
	Grid();
	~Grid();

	typedef std::vector<size_t> GridUnit;
	typedef std::vector<GridUnit> GridUnitSet;

	void init(const PCISPH::Vec3 boxSize, float unitSize);
	const PCISPH::iVec3 getSize() const { return size; }
	const float getUnitSize() const { return unitSize; }
	GridUnit& getGridUnit(const PCISPH::iVec3 gridPos);

	void insert(const PCISPH::Vec3 pos, size_t index);
	void getNeighborGridUnits(const PCISPH::iVec3 gridPos, GridUnitSet &neighbors) const;
	//void getNeighbors(const PCISPH::Vec3 pos, std::vector<size_t> &neighbors) const;
	void move(const PCISPH::Vec3 oldPos, const PCISPH::Vec3 newPos, const size_t index);

private:
	PCISPH::iVec3 size;
	float unitSize;
	PCISPH::Vec3 boxSize;
	GridUnitSet grid;

	PCISPH::iVec3 getGridPos(const PCISPH::Vec3 &pos) const;
	size_t getGridIndex(const PCISPH::iVec3 &pos) const;
};

