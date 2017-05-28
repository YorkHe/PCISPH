#include "Scene.h"
#include <fstream>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using boost::property_tree::ptree;
using boost::property_tree::read_json;

Scene::Scene(const std::string &path)
{
	this->parseFromFile(path);
}


Scene::~Scene()
{
}

template<typename T>
std::vector<T> as_vector(ptree const& pt, ptree::key_type const& key)
{
	std::vector<T> r;
	for (auto& item : pt.get_child(key))
		r.push_back(item.second.get_value<T>());

	return r;
}

void Scene::parseFromFile(const std::string &path) {
	ptree pt;
	std::fstream fs;

	try {
		fs.open(path, std::fstream::in);
		read_json(fs, pt);
		fs.close();
	}
	catch (std::ifstream::failure e) {
		std::cerr << "Can not open configuration file." << std::endl;
		throw e;
		return;
	}

	std::string modeStr = pt.get<std::string>("mode");
	if (modeStr == "FLOW") {
		this->mode = FLOW;
	}
	else {
		this->mode = STATIC;
	}
	this->particleNumber = pt.get<long>("particleNumber");
	this->referenceDensity = pt.get<float>("referenceDensity");
	this->viscosityCoefficient = pt.get<float>("viscosityCoefficient");
	this->timeStep = pt.get<float>("timeStep");
	std::vector<float> g = as_vector<float>(pt, "gravity");
	this->gravity.x = g[0];
	this->gravity.y = g[1];
	this->gravity.z = g[2];

	for (auto &item : pt.get_child("boxBounds"))
	{
		std::vector<float> f = as_vector<float>(item.second, "");
		this->boxBounds.push_back(PCISPH::Vec3(f[0], f[1], f[2]));
	}

	for (auto &item : pt.get_child("fluidBounds"))
	{
		std::vector<float> f = as_vector<float>(item.second, "");
		this->fluidBounds.push_back(PCISPH::Vec3(f[0], f[1], f[2]));
	}

	std::cout << "Read Scene" << std::endl;
	std::cout << "\tmode = " << modeStr << std::endl;
	std::cout << "\tparticleNumber = " << this->particleNumber << std::endl;
	std::cout << "\trefernceDensity = " << this->referenceDensity << std::endl;
	std::cout << "\tviscosityCoefficient = " << this->viscosityCoefficient << std::endl;
	std::cout << "\ttimeStep = " << this->timeStep << std::endl;
	std::cout << "\tgravity = [ " << this->gravity.x << ", " << gravity.y << ", " << gravity.z << " ]" << std::endl;
	std::cout << "\tboxBounds = [ " << this->boxBounds[0].x << ", " << boxBounds[0].y << ", " << boxBounds[0].z << " ] - [ ";
	std::cout << this->boxBounds[1].x << ", " << boxBounds[1].y << ", " << boxBounds[1].z << " ]" << std::endl;
	std::cout << "\tfluidBounds = [ " << this->fluidBounds[0].x << ", " << fluidBounds[0].y << ", " << fluidBounds[0].z << " ] - [ ";
	std::cout << this->fluidBounds[1].x << ", " << fluidBounds[1].y << ", " << fluidBounds[1].z << " ]" << std::endl;
}