#include "MarchingCube.h"

MarchingCube::MarchingCube(const Scene* scene, const HostParticleSet* particle)
{
	mScene = scene;
	mParticleSet = particle;

	mKernel = new Kernel();
	mKernel->init(8 * scene->particleRadius);
}


int MarchingCube::Polygonise(PCISPH::Grid grid)
{
	int cubeIndex = 0;

	if (grid.v[0]->value < isoLevel) cubeIndex |= 1;
	if (grid.v[1]->value < isoLevel) cubeIndex |= 2;
	if (grid.v[2]->value < isoLevel) cubeIndex |= 4;
	if (grid.v[3]->value < isoLevel) cubeIndex |= 8;
	if (grid.v[4]->value < isoLevel) cubeIndex |= 16;
	if (grid.v[5]->value < isoLevel) cubeIndex |= 32;
	if (grid.v[6]->value < isoLevel) cubeIndex |= 64;
	if (grid.v[7]->value < isoLevel) cubeIndex |= 128;

	if (edgeTable[cubeIndex] == 0) return 0;

	PCISPH::Vec3 vertList[12];

	if (edgeTable[cubeIndex] & 1)
		vertList[0] = vertexInterpolation(isoLevel, grid.v[0], grid.v[1]);
	if (edgeTable[cubeIndex] & 2)
		vertList[1] = vertexInterpolation(isoLevel, grid.v[1], grid.v[2]);
	if (edgeTable[cubeIndex] & 4)
		vertList[2] = vertexInterpolation(isoLevel, grid.v[2], grid.v[3]);
	if (edgeTable[cubeIndex] & 8)
		vertList[3] = vertexInterpolation(isoLevel, grid.v[3], grid.v[0]);
	if (edgeTable[cubeIndex] & 16)
		vertList[4] = vertexInterpolation(isoLevel, grid.v[4], grid.v[5]);
	if (edgeTable[cubeIndex] & 32)
		vertList[5] = vertexInterpolation(isoLevel, grid.v[5], grid.v[6]);
	if (edgeTable[cubeIndex] & 64)
		vertList[6] = vertexInterpolation(isoLevel, grid.v[6], grid.v[7]);
	if (edgeTable[cubeIndex] & 128)
		vertList[7] = vertexInterpolation(isoLevel, grid.v[7], grid.v[4]);
	if (edgeTable[cubeIndex] & 256)
		vertList[8] = vertexInterpolation(isoLevel, grid.v[0], grid.v[4]);
	if (edgeTable[cubeIndex] & 512)
		vertList[9] = vertexInterpolation(isoLevel, grid.v[1], grid.v[5]);
	if (edgeTable[cubeIndex] & 1024)
		vertList[10] = vertexInterpolation(isoLevel, grid.v[2], grid.v[6]);
	if (edgeTable[cubeIndex] & 2048)
		vertList[11] = vertexInterpolation(isoLevel, grid.v[3], grid.v[7]);

	for (int i = 0; triTable[cubeIndex][i] != -1; i += 3)
	{
		VerticePosition.push_back(vertList[triTable[cubeIndex][i]]);
		VerticePosition.push_back(vertList[triTable[cubeIndex][i + 1]]);
		VerticePosition.push_back(vertList[triTable[cubeIndex][i + 2]]);
	}

	return 1;
}

PCISPH::Vec3 MarchingCube::vertexInterpolation(double isolevel, PCISPH::Vertex* p1, PCISPH::Vertex* p2)
{
	double mu;
	PCISPH::Vec3 p;
	if (abs(isolevel - p1->value) < 0.00001)
		return p1->position;
	if (abs(isolevel - p2->value) < 0.00001)
		return p2->position;
	if (abs(p1->value - p2->value) < 0.00001)
		return p1->position;

	mu = (isolevel - p1->value) / (p2->value - p1->value);
	p.x = p1->position.x + mu * (p2->position.x - p1->position.x);
	p.y = p1->position.y + mu * (p2->position.y - p1->position.y);
	p.z = p1->position.z + mu * (p2->position.z - p1->position.z);

	return p;
}

void MarchingCube::updateParticles(const HostParticleSet* particle)
{
	mParticleSet = particle;
}


void MarchingCube::buildMesh()
{
	for (float x = 0; x < mScene->boxSize[0]; x += gridSize)
	{
		for (float y = 0; y < mScene->boxSize[1]; y += gridSize)
		{
			for (float z = 0; z < mScene->boxSize[2]; z += gridSize)
			{
				PCISPH::Vertex* v = new PCISPH::Vertex();

				v->position.x = x;
				v->position.y = y;
				v->position.z = z;

				v->value = calculateScalar(v->position);

				mGridVertices.push_back(v);
			}
		}
	}

	xLength = trunc(mScene->boxSize[0] / gridSize) + 1;
	yLength = trunc(mScene->boxSize[1] / gridSize) + 1;
	zLength = trunc(mScene->boxSize[2] / gridSize) + 1;

	for (int i = 0; i < mGridVertices.size(); i++)
	{
		PCISPH::Vertex* v = mGridVertices[i];
		if (v->position.x + gridSize >= mScene->boxSize[0] ||
			v->position.y + gridSize >= mScene->boxSize[1] ||
			v->position.z + gridSize >= mScene->boxSize[2])
		{
			continue;
		}

		PCISPH::Grid *g = new PCISPH::Grid();
		g->v[1] = mGridVertices[i];
		g->v[0] = mGridVertices[i + 1];
		g->v[2] = mGridVertices[i + zLength];
		g->v[3] = mGridVertices[i + zLength + 1];

		g->v[5] = mGridVertices[i + zLength * yLength];
		g->v[4] = mGridVertices[i + zLength * yLength + 1];
		g->v[6] = mGridVertices[i + zLength * yLength + zLength];
		g->v[7] = mGridVertices[i + zLength * yLength + zLength + 1];

		mGrids.push_back(g);
	}

	triCount = 0;
	for (auto g : mGrids)
	{
		triCount += Polygonise(*g);
	}
}




double MarchingCube::calculateScalar(PCISPH::Vec3 x)
{

	double result = 0.0;

	for (int i = 0; i < mParticleSet->count; i++)
	{
		PCISPH::Vec3 p = mParticleSet->position[i];

		result += mParticleSet->particleMass * mKernel->poly6Kernel(p - x);
	}

	return result;
}