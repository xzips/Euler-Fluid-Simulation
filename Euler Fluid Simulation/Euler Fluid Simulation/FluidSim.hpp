#pragma once

#include <iostream>
#include "SFML/Graphics.hpp"
#include <vector>

enum class FieldType {
	U_FIELD,
	V_FIELD,
	S_FIELD,
	P_FIELD,
	M_FIELD
};


enum class ObstacleType {
	SQUARE,
	CIRCLE,
	IMAGE
};


class FluidSim;

class Obstacle
{
public:
	Obstacle(float x, float y, float radius) {
		this->x = x;
		this->y = y;
		this->radius = radius;
		this->type = ObstacleType::CIRCLE;
	}

	Obstacle(float x, float y, float width, float height) {
		this->x = x;
		this->y = y;
		this->width = width;
		this->height = height;
		this->type = ObstacleType::SQUARE;
	}

	Obstacle(sf::Image image, float x, float y, float scale) {
		this->modelImage = image;
		this->type = ObstacleType::IMAGE;
		this->scale = scale;
		this->x = x;
		this->y = y;

		modelTexture.loadFromImage(modelImage);
		sprite.setTexture(modelTexture);
	}

	Obstacle() {}


	void CreateFromRect(float x, float y, float width, float height) {
		this->x = x;
		this->y = y;
		this->width = width;
		this->height = height;
		this->type = ObstacleType::SQUARE;
	}

	void CreateFromCircle(float x, float y, float radius) {
		this->x = x;
		this->y = y;
		this->radius = radius;
		this->type = ObstacleType::CIRCLE;
	}

	void CreateFromImage(sf::Image image, float x, float y, float scale) {
		this->modelImage = image;
		this->type = ObstacleType::IMAGE;
		
		modelTexture.loadFromImage(modelImage);
		sprite.setTexture(modelTexture);

	}
	

	float x = 0.5;
	float y = 0.5;
	float radius = 0.1;
	float width = 0.1;
	float height = 0.1;
	float scale = 1;
	ObstacleType type = ObstacleType::CIRCLE;
	sf::Image modelImage;
	sf::Texture modelTexture;
	sf::Sprite sprite;

	void SetObstacleSField(FluidSim* fluidSim);

	void DrawObstaclePretty(sf::RenderWindow& window);
	

	
};


class RectangularDyeSource
{
public:
	RectangularDyeSource(float x, float y, float width, float height, float density = 0.f) {
		this->x = x;
		this->y = y;
		this->width = width;
		this->height = height;
		this->density = density;
	}

	RectangularDyeSource() {}
	float x = 0.5;
	float y = 0.5;
	float width = 0.1;
	float height = 0.1;
	float density = 0.f;
};

class SimParameters
{
public:
	float density = 1000.f;
	size_t n_x = 300 + 2;
	size_t n_y = 200 + 2;
	size_t n_cells = n_x * n_y;
	
	float gridSpacingX = 1.f / (float)n_x;
	float gridSpacingY = 1.f / (float)n_y;

	float domainHeight = 1.f;
	float domainWidth = (float)n_x / (float)n_y;
	
	float gravity = 0.f;
	size_t numIterations = 100;
	float overRelaxation = 1.9;

	float windTunnelSpeed = 0.f;

	float obstacleX;
	float obstacleY;
	float obstacleRadius;

	std::vector<RectangularDyeSource> dyeSources;
	std::vector<Obstacle> obstacles;


	void SetGridSize(size_t n_x, size_t n_y) {
		this->n_x = n_x + 2;
		this->n_y = n_y + 2;
		this->n_cells = this->n_x * this->n_y;
		this->gridSpacingX = 1.f / (float)n_x;
		this->gridSpacingY = 1.f / (float)n_y;
	}

};

class DisplayParameters {
public:

	
	size_t windowWidth = 800;
	size_t windowHeight = 800;
	size_t maxFps = 60;
	size_t frameCount = 0;
	
	
	bool paused = false;
	bool showObstacle = true;
	bool showStreamlines = false;
	bool showVelocities = false;
	bool showPressure = false;
	bool showSmoke = true;
};


class FluidSim {

public:
	FluidSim(SimParameters simParams, DisplayParameters displayParams);

	void ApplyGravity(float dt, float gravity);

	void SolveIncompressibility(size_t numIterations, float dt);
	
	void Extrapolate();

	float SampleField(float x, float y, FieldType fieldType);

	void AdvectVel(float dt);

	void AdvectSmoke(float dt);

	void UpdateSField();

	void Simulate(float dt);

	bool CheckFieldsExploded();

	void Reset();

	void SetupWindTunnelBoundaries();

	void SetObstacle(float x, float y, float r);

	void ApplyDyeSources();

	void AddDyeSource(float x, float y, float width, float height, float density);

	float avgU(size_t i, size_t j);
	float avgV(size_t i, size_t j);
	
	SimParameters simParams;
	DisplayParameters displayParams;

	sf::Image image;
	sf::Texture texture;
	sf::Sprite sprite;

	void Render(sf::RenderWindow& window);

	

	float * uField; //horizontal velocity field
	float * vField; //vertical velocity field

	float * sField;

	float * newUField;
	float * newVField;

	float * pField;

	float * mField;
	float * newMField;




};