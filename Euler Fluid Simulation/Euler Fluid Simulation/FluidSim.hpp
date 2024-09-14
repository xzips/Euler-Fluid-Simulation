#include <iostream>


enum class FieldType {
	U_FIELD,
	V_FIELD,
	S_FIELD,
	P_FIELD,
	M_FIELD
};


/*

gravity : -9.81,
		dt : 1.0 / 120.0,
		numIters : 100,
		frameNr : 0,
		overRelaxation : 1.9,
		obstacleX : 0.0,
		obstacleY : 0.0,
		obstacleRadius: 0.15,
		paused: false,
		sceneNr: 0,
		showObstacle: false,
		showStreamlines: false,
		showVelocities: false,
		showPressure: false,
		showSmoke: true,
		fluid: null
*/


class Obstacle
{
	//todo later
};

class SimParameters
{
public:
	float density = 1000.f;
	size_t n_x = 200 + 2;
	size_t n_y = 200 + 2;
	size_t n_cells = n_x * n_y;
	
	float gridSpacing = 1.f / (float)n_x;
	
	float gravity = 0.f;
	size_t numIterations = 50;
	float overRelaxation = 1.9;

	float obstacleX;
	float obstacleY;
	float obstacleRadius;


	void SetGridSize(size_t n_x, size_t n_y) {
		this->n_x = n_x + 2;
		this->n_y = n_y + 2;
		this->n_cells = this->n_x * this->n_y;
		this->gridSpacing = 1.f / (float)n_x;
	}

};

class DisplayParameters {
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

	void Simulate(float dt);

	bool CheckFieldsExploded();

	void SetObstacle(float x, float y, float r);

	float avgU(size_t i, size_t j);
	float avgV(size_t i, size_t j);
	
	SimParameters simParams;
	DisplayParameters displayParams;

	


	float * uField; //horizontal velocity field
	float * vField; //vertical velocity field

	float * sField;

	float * newUField;
	float * newVField;

	float * pField;

	float * mField;
	float * newMField;




};