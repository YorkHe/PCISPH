#include "Scene.h"
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "Engine.h"
#include <iostream>


using boost::property_tree::ptree;
using boost::property_tree::read_json;

template<typename T>
std::vector<T> as_vector(ptree const& pt, ptree::key_type const& key)
{
	std::vector<T> r;
	for (auto& item : pt.get_child(key))
		r.push_back(item.second.get_value<T>());

	return r;
}

Scene::Scene()
{
}


Scene::~Scene()
{
}

void Scene::parseFromFile(std::string path)
{
	ptree pt;
	std::fstream fs;
	fs.open(path, std::fstream::in);
	read_json(fs, pt);
	fs.close();

	ptree ptConstant = pt.get_child("constant");
	ptree ptScene = pt.get_child("scene");

	constant.particleNumbers = ptConstant.get<long>("particleNumbers");
	constant.particleRadius = ptConstant.get<float>("particleRadius");
	constant.restDensity = ptConstant.get<float>("restDensity");
	constant.surfaceTension = ptConstant.get<float>("surfaceTension");
	constant.viscosity = ptConstant.get<float>("viscosity");
	constant.timeStep = ptConstant.get<float>("timeStep");

	std::vector<float> gravity = as_vector<float>(ptConstant, "gravity");
	constant.gravity[0] = gravity[0];
	constant.gravity[1] = gravity[1];
	constant.gravity[2] = gravity[2];

	std::vector<float> cameraPosition = as_vector<float>(ptScene, "cameraPosition");
	std::vector<float> cameraTarget = as_vector<float>(ptScene, "cameraTarget");

	scene.cameraPosition[0] = cameraPosition[0];
	scene.cameraPosition[1] = cameraPosition[1];
	scene.cameraPosition[2] = cameraPosition[2];

	scene.cameraTarget[0] = cameraTarget[0];
	scene.cameraTarget[1] = cameraTarget[1];
	scene.cameraTarget[2] = cameraTarget[2];

	for (auto &item : ptScene.get_child("worldBounds"))
	{
		std::vector<float> f = as_vector<float>(item.second, "");
		scene.worldBounds.push_back(PCISPH::Vec3(f[0], f[1], f[2]));
	}

	for (auto &item : ptScene.get_child("boxBounds"))
	{
		std::vector<float> f = as_vector<float>(item.second, "");
		scene.boxBounds.push_back(PCISPH::Vec3(f[0], f[1], f[2]));
	}

}
