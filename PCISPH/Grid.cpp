#include "Grid.h"



Grid::Grid()
{
}


Grid::~Grid()
{
}

void Grid::init(const PCISPH::Vec3 boxSize, float unitSize) {
	this->boxSize = boxSize;
	this->unitSize = unitSize;

	grid.clear();

	int gridNumx = ceil(boxSize.x / unitSize);
	int gridNumy = ceil(boxSize.y / unitSize);
	int gridNumz = ceil(boxSize.z / unitSize);
	this->size = PCISPH::iVec3(gridNumx, gridNumy, gridNumz);

	grid.resize(gridNumx * gridNumy * gridNumz);
}

Grid::GridUnit& Grid::getGridUnit(const PCISPH::iVec3 gridPos) {
	return this->grid[this->getGridIndex(gridPos)];
}

void Grid::insert(const PCISPH::Vec3 pos, size_t index) {
	PCISPH::iVec3 gridPos = this->getGridPos(pos);
	size_t gridIndex = this->getGridIndex(gridPos);
	this->grid[gridIndex].push_back(index);
}

void Grid::getNeighborGridUnits(const PCISPH::iVec3 gridPos, Grid::GridUnitSet &neighbors) const {
	neighbors.clear();
	for (int dx = -1; dx <= 1; dx++) {
		if (gridPos.x + dx < 0 || gridPos.x + dx >= size.x) {
			continue;
		}
		for (int dy = -1; dy <= 1; dy++) {
			if (gridPos.y + dy < 0 || gridPos.y + dy >= size.y) {
				continue;
			}
			for (int dz = -1; dz <= 1; dz++) {
				if (gridPos.z + dz < 0 || gridPos.z + dz >= size.z) {
					continue;
				}
				PCISPH::iVec3 pos = gridPos + PCISPH::iVec3(dx, dy, dz);
				size_t posIndex = this->getGridIndex(pos);
				neighbors.push_back(this->grid[posIndex]);
			}
		}
	}
}

//void Grid::getNeighbors(const PCISPH::Vec3 pos, std::vector<size_t> &neighbors) const {
//	PCISPH::iVec3 gridPos = this->getGridPos(pos);
//	size_t gridIndex = this->getGridIndex(gridPos);
//
//}

inline PCISPH::iVec3 Grid::getGridPos(const PCISPH::Vec3 &pos) const {
	return PCISPH::iVec3(pos / this->unitSize);
}

inline size_t Grid::getGridIndex(const PCISPH::iVec3 &pos) const {
	return pos.x * size.y * size.z + pos.y * size.z + pos.z;
}
