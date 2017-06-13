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

	this->referenceDensity = pt.get<float>("referenceDensity");
	this->particleRadius = pt.get<float>("particleRadius");
	this->viscosityCoefficient = pt.get<float>("viscosityCoefficient");
	this->surfaceTensionCoefficient = pt.get<float>("surfaceTensionCoefficient");
	this->gasStiffness = pt.get<float>("gasStiffness");
	this->kernelParticles = pt.get<size_t>("kernelParticles");

	this->timeStep = pt.get<float>("timeStep");

	this->restitution = pt.get<float>("restitution");
	std::vector<float> g = as_vector<float>(pt, "gravity");
	this->gravity.x = g[0];
	this->gravity.y = g[1];
	this->gravity.z = g[2];

	std::vector<float> bSize = as_vector<float>(pt, "boxSize");
	this->boxSize.x = bSize[0];
	this->boxSize.y = bSize[1];
	this->boxSize.z = bSize[2];

	std::vector<float> fSize = as_vector<float>(pt, "fluidSize");
	this->fluidSize.x = fSize[0];
	this->fluidSize.y = fSize[1];
	this->fluidSize.z = fSize[2];

	std::vector<float> fPos = as_vector<float>(pt, "fluidPosition");
	this->fluidPosition.x = fPos[0];
	this->fluidPosition.y = fPos[1];
	this->fluidPosition.z = fPos[2];

	std::vector<float> initV = as_vector<float>(pt, "initVelocity");
	this->initVelocity.x = initV[0];
	this->initVelocity.y = initV[1];
	this->initVelocity.z = initV[2];

	std::cout << "Read Scene" << std::endl;
	std::cout << "\t mode = " << modeStr << std::endl;
	std::cout << "\t refernceDensity = " << this->referenceDensity << std::endl;
	std::cout << "\t particleRadius = " << this->particleRadius << std::endl;
	std::cout << "\t viscosityCoefficient = " << this->viscosityCoefficient << std::endl;
	std::cout << "\t surfaceTensionCoefficient = " << this->surfaceTensionCoefficient << std::endl;
	std::cout << "\t gasStiffness = " << this->gasStiffness << std::endl;
	std::cout << "\t kernelParticles = " << this->kernelParticles << std::endl;

	std::cout << "\t timeStep = " << this->timeStep << std::endl;

	std::cout << "\t restitution = " << this->restitution << std::endl;
	std::cout << "\t gravity = [ " << this->gravity.x << ", " << gravity.y << ", " << gravity.z << " ]" << std::endl;
	std::cout << "\t boxSize = [ " << this->boxSize.x << ", " << boxSize.y << ", " << boxSize.z << " ]" << std::endl;
	std::cout << "\t fluidSize = [ " << this->fluidSize.x << ", " << fluidSize.y << ", " << fluidSize.z << " ]" << std::endl;
	std::cout << "\t fluidPosition = [ " << this->fluidPosition.x << ", " << fluidPosition.y << ", " << fluidPosition.z << " ]" << std::endl;
	std::cout << "\t initVelocity = [ " << this->initVelocity.x << ", " << initVelocity.y << ", " << initVelocity.z << " ]" << std::endl;
}